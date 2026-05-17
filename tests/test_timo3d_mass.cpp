#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Element/Timoshenko3D.h>
#include <BeamLib/Math/BeamMath3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <cmath>
#include <cstdio>

// Direct Timoshenko 3D mass coverage:
//   6.1 single horizontal element: every nonzero local-mass entry checked
//       against closed form; Timo-EB delta verified to be exactly the two
//       rotary inertia 2x2 blocks (one per theta DOF).
//   6.2 generic-orientation transformed mass: re-derive M_g = T^T M_l T with
//       buildTransformation3D (uses S*lambda*S) and computeMass, compare
//       entry-by-entry with BeamModel<Timoshenko3D>::assembleMass output.
//   6.3 two-element model: free-DOF reduction + shared-node accumulation.

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

struct Sec {
    double E = 200e9, G = 80e9, rho = 7800.0;
    double A = 0.01;
    double Iy = 8.0e-6, Iz = 1.6e-5;
    double Ix = 2.4e-5;
    double kappa_y = 5.0 / 6.0;
    double kappa_z = 0.85;
};

void setProps(BeamModel<Timoshenko3D>& m, const Sec& s)
{
    m.props.E = s.E; m.props.G = s.G; m.props.rho = s.rho;
    m.props.A = s.A; m.props.Iy = s.Iy; m.props.Iz = s.Iz; m.props.Ix = s.Ix;
    m.props.kappa_y = s.kappa_y; m.props.kappa_z = s.kappa_z;
}

int test_single_horizontal()
{
    Sec s;
    BeamModel<Timoshenko3D> model;
    setProps(model, s);

    const double L = 1.0;
    model.nodes.resize(2);
    model.nodes[0].x0 = Vec3(0, 0, 0);
    model.nodes[1].x0 = Vec3(L, 0, 0);
    ElementConn c; c.nodeA = 0; c.nodeB = 1; c.refVector = Vec3(0, 1, 0);
    model.elements.push_back(c);

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);
    if (M.rows() != 12) {
        std::printf("FAIL horizontal: size %lld\n", (long long)M.rows());
        return 1;
    }
    MatX Md = MatX(M);
    if (check_symmetry(Md, "horizontal")) return 1;

    const double rhoAL  = s.rho * s.A  * L;
    const double rhoIxL = s.rho * s.Ix * L;
    const double rhoIyL = s.rho * s.Iy * L;
    const double rhoIzL = s.rho * s.Iz * L;
    const double L2     = L * L;
    const double bend   = rhoAL / 420.0;

    // Axial
    if (check_rel(Md(0, 0), 2.0 * rhoAL / 6.0, 1e-12, "axial M(0,0)")) return 1;
    if (check_rel(Md(0, 6), 1.0 * rhoAL / 6.0, 1e-12, "axial M(0,6)")) return 1;

    // Torsion
    if (check_rel(Md(3, 3), 2.0 * rhoIxL / 6.0, 1e-12, "torsion M(3,3)")) return 1;
    if (check_rel(Md(3, 9), 1.0 * rhoIxL / 6.0, 1e-12, "torsion M(3,9)")) return 1;

    // xy bending Hermite translational + theta_z rotary inertia on (5,11)
    if (check_rel(Md(1, 1),  156.0 * bend,           1e-12, "xy M(1,1)")) return 1;
    if (check_rel(Md(1, 7),   54.0 * bend,           1e-12, "xy M(1,7)")) return 1;
    if (check_rel(Md(1, 5),   22.0 * L * bend,       1e-12, "xy M(1,5)")) return 1;
    if (check_rel(Md(1, 11), -13.0 * L * bend,       1e-12, "xy M(1,11)")) return 1;
    // theta_z entries get rotary inertia delta
    const double m55_exp  = 4.0 * L2 * bend + 2.0 * rhoIzL / 6.0;
    const double m511_exp = -3.0 * L2 * bend + 1.0 * rhoIzL / 6.0;
    if (check_rel(Md(5, 5),   m55_exp,  1e-12, "xy M(5,5)=Hermite+rot_z")) return 1;
    if (check_rel(Md(11, 11), m55_exp,  1e-12, "xy M(11,11)=Hermite+rot_z")) return 1;
    if (check_rel(Md(5, 11),  m511_exp, 1e-12, "xy M(5,11)=Hermite+rot_z")) return 1;

    // xz bending Hermite translational + theta_y rotary inertia on (4,10)
    if (check_rel(Md(2, 2),  156.0 * bend,           1e-12, "xz M(2,2)")) return 1;
    if (check_rel(Md(2, 8),   54.0 * bend,           1e-12, "xz M(2,8)")) return 1;
    if (check_rel(Md(2, 4),   22.0 * L * bend,       1e-12, "xz M(2,4)")) return 1;
    if (check_rel(Md(2, 10), -13.0 * L * bend,       1e-12, "xz M(2,10)")) return 1;
    const double m44_exp  = 4.0 * L2 * bend + 2.0 * rhoIyL / 6.0;
    const double m410_exp = -3.0 * L2 * bend + 1.0 * rhoIyL / 6.0;
    if (check_rel(Md(4, 4),   m44_exp,  1e-12, "xz M(4,4)=Hermite+rot_y")) return 1;
    if (check_rel(Md(10, 10), m44_exp,  1e-12, "xz M(10,10)=Hermite+rot_y")) return 1;
    if (check_rel(Md(4, 10),  m410_exp, 1e-12, "xz M(4,10)=Hermite+rot_y")) return 1;

    // Delta against EB3D mass: exactly two rotary inertia 2x2 blocks.
    auto rE = EulerBernoulli3D::computeMass(Vec3::Zero(), Vec3(L, 0, 0), model.props);
    auto rT = Timoshenko3D::computeMass(Vec3::Zero(), Vec3(L, 0, 0), model.props);
    MatMN<12, 12> delta = rT.me - rE.me;
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            double d = delta(i, j);
            double exp = 0.0;
            // theta_y rotary inertia (indices 4, 10)
            if ((i == 4 && j == 4) || (i == 10 && j == 10)) exp = 2.0 * rhoIyL / 6.0;
            else if ((i == 4 && j == 10) || (i == 10 && j == 4)) exp = 1.0 * rhoIyL / 6.0;
            // theta_z rotary inertia (indices 5, 11)
            else if ((i == 5 && j == 5) || (i == 11 && j == 11)) exp = 2.0 * rhoIzL / 6.0;
            else if ((i == 5 && j == 11) || (i == 11 && j == 5)) exp = 1.0 * rhoIzL / 6.0;
            if (std::fabs(d - exp) > 1e-12 * (1.0 + std::fabs(exp))) {
                std::printf("FAIL Timo-EB delta at (%d,%d): got %g expected %g\n",
                            i, j, d, exp);
                return 1;
            }
        }
    }

    // Rigid translation along each axis -> u^T M u = rho*A*L
    auto rigid_check = [&](int d, double exp, const char* tag) {
        VecN<12> u = VecN<12>::Zero();
        u[d] = 1.0; u[d + 6] = 1.0;
        const double m = u.dot(rT.me * u);
        return check_rel(m, exp, 1e-12, tag);
    };
    if (rigid_check(0, rhoAL, "rigid u_x")) return 1;
    if (rigid_check(1, rhoAL, "rigid u_y")) return 1;
    if (rigid_check(2, rhoAL, "rigid u_z")) return 1;
    // Rigid theta_x -> rho * Ix * L
    if (rigid_check(3, rhoIxL, "rigid theta_x")) return 1;

    return 0;
}

int test_generic_rotated_consistency()
{
    Sec s;
    BeamModel<Timoshenko3D> model;
    setProps(model, s);

    const double L = 1.5;
    const Vec3 xA(0.4, -0.3, 0.6);
    const Vec3 dir = Vec3(2.0, 1.0, -1.0).normalized();
    const Vec3 xB  = xA + L * dir;
    const Vec3 refVec(0.1, 1.0, 0.3);

    model.nodes.resize(2);
    model.nodes[0].x0 = xA;
    model.nodes[1].x0 = xB;
    ElementConn c; c.nodeA = 0; c.nodeB = 1; c.refVector = refVec;
    model.elements.push_back(c);

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);
    MatX Md = MatX(M);
    if (Md.rows() != 12) return 1;
    if (check_symmetry(Md, "generic-rotated")) return 1;

    MatMN<12, 12> T = buildTransformation3D(xA, xB, refVec);
    auto mr = Timoshenko3D::computeMass(Vec3::Zero(), Vec3(L, 0, 0), model.props);
    MatMN<12, 12> Mg_expected = T.transpose() * mr.me * T;

    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            const double got = Md(i, j);
            const double exp = Mg_expected(i, j);
            const double tol = 1e-10 * (1.0 + std::fabs(exp));
            if (std::fabs(got - exp) > tol) {
                std::printf("FAIL generic-rotated M(%d,%d)=%g vs expected %g\n",
                            i, j, got, exp);
                return 1;
            }
        }
    }
    return 0;
}

int test_two_elements_fixed_root()
{
    Sec s;
    BeamModel<Timoshenko3D> model;
    setProps(model, s);

    const double L = 1.0;
    model.nodes.resize(3);
    for (int i = 0; i < 3; ++i) model.nodes[i].x0 = Vec3(i * L, 0, 0);
    model.nodes[0].fixAll();
    for (int e = 0; e < 2; ++e) {
        ElementConn c; c.nodeA = e; c.nodeB = e + 1; c.refVector = Vec3(0, 1, 0);
        model.elements.push_back(c);
    }

    model.buildDofMap();
    SpMat M;
    model.assembleMass(M);
    if (M.rows() != 12) {
        std::printf("FAIL two-element size: %lld\n", (long long)M.rows());
        return 1;
    }
    MatX Md = MatX(M);
    if (check_symmetry(Md, "two-element")) return 1;

    const double rhoAL  = s.rho * s.A  * L;
    const double rhoIxL = s.rho * s.Ix * L;
    const double rhoIyL = s.rho * s.Iy * L;
    const double rhoIzL = s.rho * s.Iz * L;
    const double L2     = L * L;
    const double bend   = rhoAL / 420.0;

    // Free-DOF layout: node 1 -> free DOFs 0..5, node 2 -> 6..11.
    // Interior node receives both elements' contributions.
    if (check_rel(Md(0, 0), 4.0 * rhoAL / 6.0,  1e-12, "n1 axial")) return 1;
    if (check_rel(Md(6, 6), 2.0 * rhoAL / 6.0,  1e-12, "n2 axial")) return 1;
    if (check_rel(Md(3, 3), 4.0 * rhoIxL / 6.0, 1e-12, "n1 torsion")) return 1;
    if (check_rel(Md(9, 9), 2.0 * rhoIxL / 6.0, 1e-12, "n2 torsion")) return 1;
    // Theta_y interior gets 2 x (Hermite + rotary):
    const double m44_one = 4.0 * L2 * bend + 2.0 * rhoIyL / 6.0;
    if (check_rel(Md(4, 4),  2.0 * m44_one, 1e-12, "n1 theta_y diag")) return 1;
    if (check_rel(Md(10, 10), m44_one,      1e-12, "n2 theta_y diag")) return 1;
    // Theta_z interior gets 2 x (Hermite + rotary):
    const double m55_one = 4.0 * L2 * bend + 2.0 * rhoIzL / 6.0;
    if (check_rel(Md(5, 5),  2.0 * m55_one, 1e-12, "n1 theta_z diag")) return 1;
    if (check_rel(Md(11, 11), m55_one,      1e-12, "n2 theta_z diag")) return 1;

    return 0;
}

} // namespace

int main()
{
    if (test_single_horizontal())           return 1;
    if (test_generic_rotated_consistency()) return 1;
    if (test_two_elements_fixed_root())     return 1;
    std::printf("PASS test_timo3d_mass (all sub-tests)\n");
    return 0;
}
