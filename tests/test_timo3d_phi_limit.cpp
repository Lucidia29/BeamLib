#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Element/Timoshenko3D.h>
#include <cmath>
#include <cstdio>

// Element-level Phi -> 0 degeneration: as Phi_y, Phi_z -> 0 the Timoshenko 3D
// local stiffness must converge to the EB3D local stiffness entry-by-entry.
//
// This is the 3D counterpart to test_timo2d_phi_limit. Choose geometry such
// that both Phi_y and Phi_z are very small (~1e-9), then compare all 144
// entries of ke against the EB3D value to relative 1e-7.
//
// Important: use Iy != Iz and kappa_y != kappa_z so a wrong plane-mapping
// (e.g. kappa_y where kappa_z should be) would not coincidentally satisfy
// the test. Two moderate-Phi spot-checks then verify the sign of the
// (2 - Phi) entry in EACH plane independently, catching plane swaps.
int main()
{
    using namespace beamlib;

    SectionProperties props;
    props.E       = 200e9;
    props.G       = 80e9;
    props.A       = 1.0;
    props.Iy      = 1.0e-3;
    props.Iz      = 3.0e-3;        // !! different from Iy on purpose
    props.Ix      = props.Iy + props.Iz;
    props.kappa_y = 5.0 / 6.0;
    props.kappa_z = 0.9;           // !! different from kappa_y on purpose

    // Choose L large so both Phi values are ~ 1e-9 or below.
    //   Phi_y = 12 EIz / (kappa_y GA L^2)
    //   Phi_z = 12 EIy / (kappa_z GA L^2)
    // Pick L = 7000 m -> Phi_y ~ 1.6e-9, Phi_z ~ 5.7e-10, both tiny.
    const double L = 7000.0;
    const Vec3 xA(0, 0, 0);
    const Vec3 xB(L, 0, 0);

    const double EIy  = props.E * props.Iy;
    const double EIz  = props.E * props.Iz;
    const double GAky = props.kappa_y * props.G * props.A;
    const double GAkz = props.kappa_z * props.G * props.A;
    const double Phi_y = 12.0 * EIz / (GAky * L * L);
    const double Phi_z = 12.0 * EIy / (GAkz * L * L);
    if (Phi_y > 1e-7 || Phi_z > 1e-7) {
        std::printf("FAIL setup: Phi_y=%g Phi_z=%g not small enough\n", Phi_y, Phi_z);
        return 1;
    }

    VecN<12> disp = VecN<12>::Zero();   // linear element: ke independent of disp
    auto rT = Timoshenko3D::computeElement(xA, xB, disp, props);
    auto rE = EulerBernoulli3D::computeElement(xA, xB, disp, props);

    const double tol = 1e-7;
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            const double t = rT.ke(i, j);
            const double e = rE.ke(i, j);
            const double denom = std::fabs(e) > 1e-20 ? std::fabs(e) : 1.0;
            const double err = std::fabs(t - e) / denom;
            if (err > tol) {
                std::printf("FAIL: ke(%d,%d) Timo=%g EB=%g rel err %g (tol %g)\n",
                            i, j, t, e, err, tol);
                std::printf("  Phi_y=%g, Phi_z=%g\n", Phi_y, Phi_z);
                return 1;
            }
        }
    }

    // Moderate-Phi spot-check on the (2 - Phi) entries, one per plane.
    //
    // xz plane uses indices (4, 10) for the theta_y - theta_y off-diagonal,
    // value (2 - Phi_z) L^2 * EIy / (L^3 (1 + Phi_z)). Decreasing L makes
    // Phi_z larger; using L = 50, Phi_z ~ 4.4e-7 (tiny but nonzero).
    {
        const double L2 = 50.0;
        const Vec3 xB2(L2, 0, 0);
        auto r2 = Timoshenko3D::computeElement(xA, xB2, disp, props);
        const double Phi_z_2 = 12.0 * EIy / (GAkz * L2 * L2);
        const double expected_4_10 =
            (2.0 - Phi_z_2) * L2 * L2 * EIy / (L2 * L2 * L2 * (1.0 + Phi_z_2));
        if (std::fabs(r2.ke(4, 10) - expected_4_10) / std::fabs(expected_4_10) > 1e-12) {
            std::printf("FAIL moderate-Phi xz: ke(4,10)=%g, expected %g\n",
                        r2.ke(4, 10), expected_4_10);
            return 1;
        }
        // Must be strictly less than the EB equivalent 2 EIy / L for any
        // nonzero Phi.
        if (r2.ke(4, 10) >= 2.0 * EIy / L2) {
            std::printf("FAIL xz: ke(4,10)=%g should be < 2 EIy/L=%g for Phi>0\n",
                        r2.ke(4, 10), 2.0 * EIy / L2);
            return 1;
        }
    }

    // xy plane uses indices (5, 11) for theta_z - theta_z off-diagonal,
    // value (2 - Phi_y) L^2 * EIz / (L^3 (1 + Phi_y)).
    {
        const double L2 = 50.0;
        const Vec3 xB2(L2, 0, 0);
        auto r2 = Timoshenko3D::computeElement(xA, xB2, disp, props);
        const double Phi_y_2 = 12.0 * EIz / (GAky * L2 * L2);
        const double expected_5_11 =
            (2.0 - Phi_y_2) * L2 * L2 * EIz / (L2 * L2 * L2 * (1.0 + Phi_y_2));
        if (std::fabs(r2.ke(5, 11) - expected_5_11) / std::fabs(expected_5_11) > 1e-12) {
            std::printf("FAIL moderate-Phi xy: ke(5,11)=%g, expected %g\n",
                        r2.ke(5, 11), expected_5_11);
            return 1;
        }
        if (r2.ke(5, 11) >= 2.0 * EIz / L2) {
            std::printf("FAIL xy: ke(5,11)=%g should be < 2 EIz/L=%g for Phi>0\n",
                        r2.ke(5, 11), 2.0 * EIz / L2);
            return 1;
        }
    }

    std::printf("PASS test_timo3d_phi_limit: Phi_y=%.3e Phi_z=%.3e; all 144 ke "
                "entries match EB3D to rel<%g; (2-Phi) signs verified in both planes.\n",
                Phi_y, Phi_z, tol);
    return 0;
}
