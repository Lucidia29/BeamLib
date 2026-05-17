#include <BeamLib/Element/Timoshenko3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Two independent bending/shear planes verified separately on the same
// cantilever geometry:
//   Subtest 1 (F_z at tip, xz plane):
//     u_z(L)    = F_z L^3 / (3 E I_y) + F_z L / (kappa_z G A)
//     theta_y(L)= F_z L^2 / (2 E I_y)
//   Subtest 2 (F_y at tip, xy plane):
//     u_y(L)    = F_y L^3 / (3 E I_z) + F_y L / (kappa_y G A)
//     theta_z(L)= F_y L^2 / (2 E I_z)
//
// Iy != Iz and kappa_y != kappa_z so the two planes are independently
// distinguishable. This catches any swap between (I_y <-> I_z) or
// (kappa_y <-> kappa_z) at the element level.

namespace {

using namespace beamlib;

struct Sec {
    double E       = 200e9;
    double G       = 80e9;
    double rho     = 7800.0;
    double A       = 0.01;
    double Iy      = 8.0e-6;       // intentionally different from Iz
    double Iz      = 1.6e-5;
    double Ix      = 2.4e-5;
    double kappa_y = 5.0 / 6.0;
    double kappa_z = 0.85;          // intentionally != kappa_y
};

BeamModel<Timoshenko3D> buildCantilever(const Sec& s, double Lbeam, int nElems)
{
    BeamModel<Timoshenko3D> model;
    model.props.E       = s.E;
    model.props.G       = s.G;
    model.props.rho     = s.rho;
    model.props.A       = s.A;
    model.props.Iy      = s.Iy;
    model.props.Iz      = s.Iz;
    model.props.Ix      = s.Ix;
    model.props.kappa_y = s.kappa_y;
    model.props.kappa_z = s.kappa_z;

    const double dx = Lbeam / nElems;
    model.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        model.nodes[i].x0 = Vec3(i * dx, 0.0, 0.0);
    }
    model.nodes[0].fixAll();
    for (int e = 0; e < nElems; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        model.elements.push_back(c);
    }
    return model;
}

int check_rel(double got, double expected, double tol, const char* tag)
{
    const double denom = std::fabs(expected) > 1e-30 ? std::fabs(expected) : 1.0;
    const double err = std::fabs(got - expected) / denom;
    if (err > tol) {
        std::printf("FAIL [%s]: got %g expected %g (rel %g, tol %g)\n",
                    tag, got, expected, err, tol);
        return 1;
    }
    return 0;
}

int check_near_zero(double got, double tol, const char* tag)
{
    if (std::fabs(got) > tol) {
        std::printf("FAIL [%s]: expected ~0, got %g (tol %g)\n", tag, got, tol);
        return 1;
    }
    return 0;
}

int subtest_Fz()
{
    Sec s;
    const double Lbeam = 1.0;
    const int    nElems = 10;
    auto model = buildCantilever(s, Lbeam, nElems);

    const double Fz = -1000.0;
    model.nodes[nElems].load[2] = Fz;

    model.buildDofMap();
    NRConfig cfg; cfg.tol = 1e-7; cfg.maxIter = 5;
    NRResult res =
        NewtonRaphsonSolver<Timoshenko3D>::solveOneStep(model, model.getExternalForceVector(), cfg);
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL Fz: NR converged=%d iters=%d\n",
                    int(res.converged), res.iterations);
        return 1;
    }

    const double EIy = s.E * s.Iy;
    const double GAk = s.kappa_z * s.G * s.A;
    const double uz_expect = Fz * Lbeam * Lbeam * Lbeam / (3.0 * EIy) + Fz * Lbeam / GAk;
    const double th_expect = Fz * Lbeam * Lbeam / (2.0 * EIy);

    const auto& tip = model.nodes[nElems].dof;
    if (check_rel(tip[2], uz_expect, 1e-9, "Fz tip u_z")) return 1;
    if (check_rel(tip[4], th_expect, 1e-9, "Fz tip theta_y")) return 1;
    // Unrelated DOFs must stay zero.
    if (check_near_zero(tip[0], 1e-14, "Fz tip u_x")) return 1;
    if (check_near_zero(tip[1], 1e-14, "Fz tip u_y")) return 1;
    if (check_near_zero(tip[3], 1e-14, "Fz tip theta_x")) return 1;
    if (check_near_zero(tip[5], 1e-14, "Fz tip theta_z")) return 1;

    std::printf("PASS subtest_Fz (xz plane, EIy=%g, kappa_z=%g): u_z=%g, theta_y=%g\n",
                EIy, s.kappa_z, tip[2], tip[4]);
    return 0;
}

int subtest_Fy()
{
    Sec s;
    const double Lbeam = 1.0;
    const int    nElems = 10;
    auto model = buildCantilever(s, Lbeam, nElems);

    const double Fy = -1000.0;
    model.nodes[nElems].load[1] = Fy;

    model.buildDofMap();
    NRConfig cfg; cfg.tol = 1e-7; cfg.maxIter = 5;
    NRResult res =
        NewtonRaphsonSolver<Timoshenko3D>::solveOneStep(model, model.getExternalForceVector(), cfg);
    if (!res.converged || res.iterations != 1) {
        std::printf("FAIL Fy: NR converged=%d iters=%d\n",
                    int(res.converged), res.iterations);
        return 1;
    }

    const double EIz = s.E * s.Iz;
    const double GAk = s.kappa_y * s.G * s.A;
    const double uy_expect = Fy * Lbeam * Lbeam * Lbeam / (3.0 * EIz) + Fy * Lbeam / GAk;
    const double th_expect = Fy * Lbeam * Lbeam / (2.0 * EIz);

    const auto& tip = model.nodes[nElems].dof;
    if (check_rel(tip[1], uy_expect, 1e-9, "Fy tip u_y")) return 1;
    if (check_rel(tip[5], th_expect, 1e-9, "Fy tip theta_z")) return 1;
    if (check_near_zero(tip[0], 1e-14, "Fy tip u_x")) return 1;
    if (check_near_zero(tip[2], 1e-14, "Fy tip u_z")) return 1;
    if (check_near_zero(tip[3], 1e-14, "Fy tip theta_x")) return 1;
    if (check_near_zero(tip[4], 1e-14, "Fy tip theta_y")) return 1;

    std::printf("PASS subtest_Fy (xy plane, EIz=%g, kappa_y=%g): u_y=%g, theta_z=%g\n",
                EIz, s.kappa_y, tip[1], tip[5]);
    return 0;
}

} // namespace

int main()
{
    if (subtest_Fz()) return 1;
    if (subtest_Fy()) return 1;
    std::printf("PASS test_timo3d_cantilever (both planes)\n");
    return 0;
}
