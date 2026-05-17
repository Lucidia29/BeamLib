#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Element/Timoshenko2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <cmath>
#include <cstdio>

namespace {

using namespace beamlib;

int check_symmetry(const MatX& M, const char* tag)
{
    const int n = static_cast<int>(M.rows());
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            const double d = std::fabs(M(i, j) - M(j, i));
            const double s = 1.0 + std::fabs(M(i, j));
            if (d > 1e-12 * s) {
                std::printf("FAIL [%s] symmetry M(%d,%d)=%g vs M(%d,%d)=%g\n",
                            tag, i, j, M(i, j), j, i, M(j, i));
                return 1;
            }
        }
    }
    return 0;
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

// Sub-test 1: horizontal single element, lambda = I so M_g = M_l directly.
// Verify all Timoshenko local mass entries and confirm the Timo-EB delta
// equals exactly rho*I*L/6 * [[2,1],[1,2]] on the (2,5) sub-block.
int test_single_horizontal()
{
    SectionProperties props;
    props.E       = 200e9;
    props.G       = 80e9;
    props.rho     = 7800.0;
    props.A       = 0.01;
    props.Iz      = 8.333e-6;
    props.kappa_z = 5.0 / 6.0;

    const double L = 1.0;
    const Vec3 xA(0, 0, 0);
    const Vec3 xB(L, 0, 0);

    auto mr = Timoshenko2D::computeMass(xA, xB, props);
    MatMN<6, 6> M = mr.me;

    MatX Md = MatX(M);
    if (check_symmetry(Md, "single-horizontal")) return 1;

    const double rhoAL = props.rho * props.A  * L;
    const double rhoIL = props.rho * props.Iz * L;
    const double bend  = rhoAL / 420.0;

    // Axial block (0, 3): same as EB
    if (check_rel(M(0, 0), 2.0 * rhoAL / 6.0, 1e-12, "M(0,0) axial")) return 1;
    if (check_rel(M(3, 3), 2.0 * rhoAL / 6.0, 1e-12, "M(3,3) axial")) return 1;
    if (check_rel(M(0, 3), 1.0 * rhoAL / 6.0, 1e-12, "M(0,3) axial off")) return 1;

    // Translational Hermite block on (1, 4) and u-theta cross terms
    if (check_rel(M(1, 1), 156.0 * bend, 1e-12, "M(1,1)")) return 1;
    if (check_rel(M(1, 4),  54.0 * bend, 1e-12, "M(1,4)")) return 1;
    if (check_rel(M(1, 2),  22.0 * L * bend, 1e-12, "M(1,2)")) return 1;
    if (check_rel(M(1, 5), -13.0 * L * bend, 1e-12, "M(1,5)")) return 1;
    if (check_rel(M(2, 4),  13.0 * L * bend, 1e-12, "M(2,4)")) return 1;
    if (check_rel(M(4, 5), -22.0 * L * bend, 1e-12, "M(4,5)")) return 1;

    // Theta-theta block: Hermite translational + rotary inertia
    // M(2,2) = 4 L^2 rhoAL / 420   +   rho*I*L / 3
    // M(2,5) = -3 L^2 rhoAL / 420  +   rho*I*L / 6
    const double m22_expect = 4.0 * L * L * bend + rhoIL / 3.0;
    const double m25_expect = -3.0 * L * L * bend + rhoIL / 6.0;
    const double m55_expect = 4.0 * L * L * bend + rhoIL / 3.0;
    if (check_rel(M(2, 2), m22_expect, 1e-12, "M(2,2) = Hermite + rotary")) return 1;
    if (check_rel(M(5, 5), m55_expect, 1e-12, "M(5,5) = Hermite + rotary")) return 1;
    if (check_rel(M(2, 5), m25_expect, 1e-12, "M(2,5) = Hermite + rotary")) return 1;

    // Cross-check: Timoshenko mass delta against EB2D mass must equal exactly
    // rho*I*L/6 * [[2,1],[1,2]] on the (2,5) sub-block, zero elsewhere.
    auto eb = EulerBernoulli2D::computeMass(xA, xB, props);
    MatMN<6, 6> delta = mr.me - eb.me;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            double d = delta(i, j);
            double exp = 0.0;
            if (i == 2 && j == 2) exp = 2.0 * rhoIL / 6.0;
            else if (i == 5 && j == 5) exp = 2.0 * rhoIL / 6.0;
            else if ((i == 2 && j == 5) || (i == 5 && j == 2)) exp = 1.0 * rhoIL / 6.0;
            if (std::fabs(d - exp) > 1e-12 * (1.0 + std::fabs(exp))) {
                std::printf("FAIL: Timo-EB mass delta at (%d,%d): got %g expected %g\n",
                            i, j, d, exp);
                return 1;
            }
        }
    }

    // Rigid axial translation: u^T M u with u_x = 1 at both nodes -> total mass rho*A*L
    VecN<6> u = VecN<6>::Zero();
    u[0] = 1.0; u[3] = 1.0;
    const double rigid_axial = u.dot(mr.me * u);
    if (check_rel(rigid_axial, rhoAL, 1e-12, "rigid axial translation")) return 1;

    return 0;
}

// Sub-test 2: transformed mass. Vertical (+z) beam: local axial mass should
// appear on the global u_z diagonal, local bending translational mass on
// global u_x diagonal, rotary inertia term unaffected (theta_y invariant).
int test_vertical_transformed()
{
    BeamModel<Timoshenko2D> model;
    model.props.E       = 200e9;
    model.props.G       = 80e9;
    model.props.rho     = 7800.0;
    model.props.A       = 0.01;
    model.props.Iz      = 8.333e-6;
    model.props.kappa_z = 5.0 / 6.0;

    const double L = 1.0;
    model.nodes.resize(2);
    model.nodes[0].x0 = Vec3(0, 0, 0);
    model.nodes[1].x0 = Vec3(0, 0, L);   // vertical, +z direction
    ElementConn c;
    c.nodeA = 0;
    c.nodeB = 1;
    model.elements.push_back(c);

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);
    MatX Md = MatX(M);
    if (check_symmetry(Md, "vertical")) return 1;

    const double rhoAL = model.props.rho * model.props.A  * L;
    const double rhoIL = model.props.rho * model.props.Iz * L;
    const double bend  = rhoAL / 420.0;

    // T_n for vertical = [[0,-1,0],[1,0,0],[0,0,1]]. After T^T M_l T:
    //   global u_x  <- local u_z (Hermite translational): 156 bend
    //   global u_z  <- local u_x (axial):                  2 rhoAL/6
    //   global theta_y unchanged: 4 L^2 bend + rhoIL/3 (full Timoshenko value)
    if (check_rel(Md(0, 0), 156.0 * bend, 1e-12, "vert: local bend -> global u_x"))
        return 1;
    if (check_rel(Md(1, 1), 2.0 * rhoAL / 6.0, 1e-12, "vert: local axial -> global u_z"))
        return 1;
    const double th_expect = 4.0 * L * L * bend + rhoIL / 3.0;
    if (check_rel(Md(2, 2), th_expect, 1e-12, "vert: theta_y invariant (Timoshenko)"))
        return 1;

    return 0;
}

} // namespace

int main()
{
    if (test_single_horizontal())     return 1;
    if (test_vertical_transformed())  return 1;
    std::printf("PASS test_timo2d_mass (all sub-tests)\n");
    return 0;
}
