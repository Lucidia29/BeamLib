#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Element/Timoshenko2D.h>
#include <cmath>
#include <cstdio>

// Element-level Phi -> 0 degeneration: the local stiffness produced by
// Timoshenko2D::computeElement must converge to the EB2D bending block as
// Phi = 12 E I / (kappa G A L^2) -> 0.
//
// Construct geometry / material so that Phi ~ 1e-9, then compare every nonzero
// entry of ke against the corresponding EB2D entry. This catches:
//   - wrong Phi algebra
//   - wrong kappa source (kappa_y vs kappa_z)
//   - wrong sign on the (2 - Phi) entry
//   - misplaced (1 + Phi) denominator
int main()
{
    using namespace beamlib;

    // Pick a slender geometry with high G to drive Phi ~ 1e-9.
    SectionProperties props;
    props.E       = 200e9;
    props.G       = 80e9;
    props.A       = 1.0;                       // unit area
    props.Iz      = 1.0e-3;                    // small I
    props.kappa_z = 5.0 / 6.0;

    // L chosen so that Phi = 12 EI / (kappa G A L^2) ~ 1e-9.
    //   Phi = 12 * 200e9 * 1e-3 / (5/6 * 80e9 * 1 * L^2)
    //       = 2.4e9 / (6.667e10 * L^2)
    //       = 0.036 / L^2
    //   Phi = 1e-9  =>  L^2 = 3.6e7  =>  L = 6000
    const double L = 6000.0;
    const Vec3 xA(0, 0, 0);
    const Vec3 xB(L, 0, 0);

    const double EI  = props.E  * props.Iz;
    const double GAk = props.kappa_z * props.G * props.A;
    const double Phi = 12.0 * EI / (GAk * L * L);

    if (Phi > 1e-7) {
        std::printf("FAIL setup: Phi=%g is not small enough for this test\n", Phi);
        return 1;
    }

    VecN<6> disp = VecN<6>::Zero();    // ke does not depend on dispVec for linear element
    auto rT = Timoshenko2D::computeElement(xA, xB, disp, props);
    auto rE = EulerBernoulli2D::computeElement(xA, xB, disp, props);

    // Entry-by-entry comparison. Tolerance has to accommodate the residual
    // Phi correction itself (~1e-9), so we use 1e-7 relative.
    const double tol = 1e-7;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            const double t = rT.ke(i, j);
            const double e = rE.ke(i, j);
            const double denom = std::fabs(e) > 1e-20 ? std::fabs(e) : 1.0;
            const double err = std::fabs(t - e) / denom;
            if (err > tol) {
                std::printf("FAIL: ke(%d,%d) Timo=%g EB=%g rel err %g (tol %g, Phi=%g)\n",
                            i, j, t, e, err, tol, Phi);
                return 1;
            }
        }
    }

    // Spot-check: the (2-Phi) entry signs are correct. At Phi -> 0 the (2,5)
    // entry equals +2 EI/L (EB's 2 EI/L). A wrong (2 + Phi) implementation
    // would still give 2 EI/L at Phi=0 but the sign would be wrong for
    // larger Phi; check that explicitly at a moderate Phi too.
    {
        SectionProperties propsM = props;
        const double L2 = 100.0;       // shorter beam -> moderate Phi
        const Vec3 xB2(L2, 0, 0);
        auto r2 = Timoshenko2D::computeElement(xA, xB2, disp, propsM);
        const double Phi2 = 12.0 * EI / (GAk * L2 * L2);
        const double expected_25 = (2.0 - Phi2) * L2 * L2 * EI / (L2 * L2 * L2 * (1.0 + Phi2));
        const double got_25 = r2.ke(2, 5);
        const double err25 = std::fabs(got_25 - expected_25) / std::fabs(expected_25);
        if (err25 > 1e-12) {
            std::printf("FAIL moderate-Phi: ke(2,5)=%g, expected %g (Phi=%g, rel %g)\n",
                        got_25, expected_25, Phi2, err25);
            return 1;
        }
        // Make sure the sign reflects (2 - Phi), not (2 + Phi): for our setup,
        // Phi2 ~ 0.0036 < 2, so the entry is still positive, but it must be
        // smaller than the EB equivalent 2 EI/L.
        const double eb_25 = 2.0 * EI / L2;
        if (got_25 >= eb_25) {
            std::printf("FAIL: ke(2,5)=%g should be smaller than EB's 2EI/L=%g for Phi>0\n",
                        got_25, eb_25);
            return 1;
        }
    }

    std::printf("PASS test_timo2d_phi_limit: Phi=%.3e, all 36 ke entries agree with EB2D "
                "to rel < %g; (2,5) sign behaves as (2 - Phi).\n", Phi, tol);
    return 0;
}
