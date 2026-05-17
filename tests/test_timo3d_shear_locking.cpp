#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Element/Timoshenko3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// 3D shear-locking check: very slender beam (L/h = 1000) with only 4
// elements, loaded simultaneously in BOTH bending planes (F_y and F_z). A
// locking element would over-stiffen and give tip deflections orders of
// magnitude smaller than the analytical EB solution.
//
// Both bending planes must remain unlocked; this is the 3D counterpart to
// the Timo2D shear-locking test.
int main()
{
    using namespace beamlib;

    const double E       = 200e9;
    const double nu      = 0.3;
    const double G       = E / (2.0 * (1.0 + nu));
    const double rho     = 7800.0;

    const double Lbeam   = 1.0;
    const double h       = Lbeam / 1000.0;
    const double b       = 0.01;
    const double A       = b * h;
    const double Iy      = b * h * h * h / 12.0;
    const double Iz      = h * b * b * b / 12.0;
    const double Ix      = Iy + Iz;
    const double kappa_y = 5.0 / 6.0;
    const double kappa_z = 5.0 / 6.0;

    const int    nElems = 4;     // intentionally coarse
    const double dx     = Lbeam / nElems;
    const double Fy     = -1e-3;
    const double Fz     = -1e-3;

    BeamModel<Timoshenko3D> m;
    m.props.E = E; m.props.G = G; m.props.rho = rho;
    m.props.A = A; m.props.Iy = Iy; m.props.Iz = Iz; m.props.Ix = Ix;
    m.props.kappa_y = kappa_y; m.props.kappa_z = kappa_z;
    m.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) m.nodes[i].x0 = Vec3(i * dx, 0, 0);
    m.nodes[0].fixAll();
    m.nodes[nElems].load[1] = Fy;
    m.nodes[nElems].load[2] = Fz;
    for (int e = 0; e < nElems; ++e) {
        ElementConn c; c.nodeA = e; c.nodeB = e + 1;
        m.elements.push_back(c);
    }
    m.buildDofMap();
    NRResult res =
        NewtonRaphsonSolver<Timoshenko3D>::solveOneStep(m, m.getExternalForceVector());
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL NR: converged=%d iters=%d\n",
                    int(res.converged), res.iterations);
        return 1;
    }

    // Analytical Timoshenko cantilever solutions for each plane.
    const double EIy = E * Iy;
    const double EIz = E * Iz;
    const double GAky = kappa_y * G * A;
    const double GAkz = kappa_z * G * A;

    const double uz_an = Fz * Lbeam * Lbeam * Lbeam / (3.0 * EIy) + Fz * Lbeam / GAkz;
    const double uy_an = Fy * Lbeam * Lbeam * Lbeam / (3.0 * EIz) + Fy * Lbeam / GAky;
    const double thy_an = Fz * Lbeam * Lbeam / (2.0 * EIy);
    const double thz_an = Fy * Lbeam * Lbeam / (2.0 * EIz);

    const auto& tip = m.nodes[nElems].dof;

    auto check_rel = [](double a, double b, double tol, const char* tag) {
        const double denom = std::fabs(b) > 1e-30 ? std::fabs(b) : 1.0;
        const double err   = std::fabs(a - b) / denom;
        if (err > tol) {
            std::printf("FAIL [LOCKING SUSPECTED %s]: got %g, analytical %g (rel %g, tol %g)\n",
                        tag, a, b, err, tol);
            std::printf("   Ratio FEM/analytical = %.6f (1.0 = no locking; <<1 = locking)\n",
                        a / b);
            return 1;
        }
        return 0;
    };

    const double tol = 1e-6;
    if (check_rel(tip[1], uy_an,  tol, "u_y"))      return 1;
    if (check_rel(tip[2], uz_an,  tol, "u_z"))      return 1;
    if (check_rel(tip[4], thy_an, tol, "theta_y"))  return 1;
    if (check_rel(tip[5], thz_an, tol, "theta_z"))  return 1;

    std::printf("PASS test_timo3d_shear_locking: L/h=%.0f, %d elements; "
                "u_y rel=%.2e, u_z rel=%.2e. Both planes locking-free.\n",
                Lbeam / h, nElems,
                std::fabs(tip[1] - uy_an) / std::fabs(uy_an),
                std::fabs(tip[2] - uz_an) / std::fabs(uz_an));
    return 0;
}
