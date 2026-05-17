#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

namespace {

using namespace beamlib;

struct CantileverProps {
    double E   = 200e9;
    double G   = 80e9;
    double A   = 0.01;
    double Iy  = 8.333e-6;
    double Iz  = 8.333e-6;
    double Ix  = 16.667e-6;
};

struct CantileverModel {
    BeamModel<EulerBernoulli3D> model;
    int                          nElems;
};

CantileverModel buildCantilever(const CantileverProps& p, double Lbeam,
                                int nElems)
{
    CantileverModel cm;
    cm.nElems = nElems;
    auto& model = cm.model;

    model.props.E  = p.E;
    model.props.G  = p.G;
    model.props.A  = p.A;
    model.props.Iy = p.Iy;
    model.props.Iz = p.Iz;
    model.props.Ix = p.Ix;

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
    return cm;
}

int check_relative(double got, double expected, double tol,
                   const char* tag)
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

// Case 1: cantilever along global x, tip force F_z = -1000 (downward).
// Expected analytical (EB cantilever in xz plane, stiffness E*Iy):
//   u_z(L)    = F_z * L^3 / (3 E Iy)   (negative for F_z < 0)
//   theta_y(L)= F_z * L^2 / (2 E Iy)   (negative; BeamLib structural convention
//                                        theta_y = d u_z / d x, so slope negative)
// All other DOFs at tip should be ~0.
int case1_tip_Fz()
{
    CantileverProps p;
    const double Lbeam = 1.0;
    const int nElems = 10;
    CantileverModel cm = buildCantilever(p, Lbeam, nElems);
    auto& model = cm.model;

    const double Fz = -1000.0;
    model.nodes[nElems].load[2] = Fz;   // index 2 = u_z

    model.buildDofMap();
    VecX Fext = model.getExternalForceVector();
    NRConfig cfg;
    cfg.tol     = 1e-7;
    cfg.maxIter = 5;
    NRResult res =
        NewtonRaphsonSolver<EulerBernoulli3D>::solveOneStep(model, Fext, cfg);

    if (!res.converged) {
        std::printf("FAIL case1: NR did not converge (iters=%d res=%g)\n",
                    res.iterations, res.finalResidual);
        return 1;
    }
    if (res.iterations != 1) {
        std::printf("FAIL case1: expected 1 NR correction step, got %d\n",
                    res.iterations);
        return 1;
    }

    const double EIy = p.E * p.Iy;
    const double uz_expected  = Fz * Lbeam * Lbeam * Lbeam / (3.0 * EIy);
    const double thy_expected = Fz * Lbeam * Lbeam / (2.0 * EIy);

    const auto& tip = model.nodes[nElems].dof;
    const double ux = tip[0], uy = tip[1], uz = tip[2];
    const double tx = tip[3], ty = tip[4], tz = tip[5];

    if (check_relative(uz, uz_expected, 1e-9, "tip u_z")) return 1;
    if (check_relative(ty, thy_expected, 1e-9, "tip theta_y")) return 1;
    if (uz >= 0.0) {
        std::printf("FAIL case1: tip u_z must be negative, got %g\n", uz);
        return 1;
    }
    if (ty >= 0.0) {
        std::printf("FAIL case1: tip theta_y must be negative, got %g\n", ty);
        return 1;
    }
    // Unrelated DOFs: u_x (no axial load), u_y, theta_x, theta_z all ~0.
    if (check_near_zero(ux, 1e-14, "tip u_x"))     return 1;
    if (check_near_zero(uy, 1e-14, "tip u_y"))     return 1;
    if (check_near_zero(tx, 1e-14, "tip theta_x")) return 1;
    if (check_near_zero(tz, 1e-14, "tip theta_z")) return 1;

    std::printf("PASS case1 tip F_z: u_z=%g, theta_y=%g, iters=%d\n",
                uz, ty, res.iterations);
    return 0;
}

// Case 2: same cantilever, tip force F_y = -1000.
// Expected (xy-plane bending, stiffness E*Iz):
//   u_y(L)    = F_y * L^3 / (3 E Iz)             (negative)
//   theta_z(L)= F_y * L^2 / (2 E Iz)             (negative; structural
//                                                  theta_z = d u_y / d x)
// All other DOFs at tip ~0.
int case2_tip_Fy()
{
    CantileverProps p;
    const double Lbeam = 1.0;
    const int nElems = 10;
    CantileverModel cm = buildCantilever(p, Lbeam, nElems);
    auto& model = cm.model;

    const double Fy = -1000.0;
    model.nodes[nElems].load[1] = Fy;   // index 1 = u_y

    model.buildDofMap();
    VecX Fext = model.getExternalForceVector();
    NRConfig cfg;
    cfg.tol     = 1e-7;
    cfg.maxIter = 5;
    NRResult res =
        NewtonRaphsonSolver<EulerBernoulli3D>::solveOneStep(model, Fext, cfg);

    if (!res.converged) {
        std::printf("FAIL case2: NR did not converge (iters=%d res=%g)\n",
                    res.iterations, res.finalResidual);
        return 1;
    }
    if (res.iterations != 1) {
        std::printf("FAIL case2: expected 1 NR correction step, got %d\n",
                    res.iterations);
        return 1;
    }

    const double EIz = p.E * p.Iz;
    const double uy_expected  = Fy * Lbeam * Lbeam * Lbeam / (3.0 * EIz);
    const double thz_expected = Fy * Lbeam * Lbeam / (2.0 * EIz);

    const auto& tip = model.nodes[nElems].dof;
    const double ux = tip[0], uy = tip[1], uz = tip[2];
    const double tx = tip[3], ty = tip[4], tz = tip[5];

    if (check_relative(uy, uy_expected, 1e-9, "tip u_y")) return 1;
    if (check_relative(tz, thz_expected, 1e-9, "tip theta_z")) return 1;
    if (uy >= 0.0) {
        std::printf("FAIL case2: tip u_y must be negative, got %g\n", uy);
        return 1;
    }
    if (tz >= 0.0) {
        std::printf("FAIL case2: tip theta_z must be negative, got %g\n", tz);
        return 1;
    }
    if (check_near_zero(ux, 1e-14, "tip u_x"))     return 1;
    if (check_near_zero(uz, 1e-14, "tip u_z"))     return 1;
    if (check_near_zero(tx, 1e-14, "tip theta_x")) return 1;
    if (check_near_zero(ty, 1e-14, "tip theta_y")) return 1;

    std::printf("PASS case2 tip F_y: u_y=%g, theta_z=%g, iters=%d\n",
                uy, tz, res.iterations);
    return 0;
}

} // namespace

int main()
{
    if (case1_tip_Fz()) return 1;
    if (case2_tip_Fy()) return 1;
    std::printf("PASS test_eb3d_cantilever (case1 + case2)\n");
    return 0;
}
