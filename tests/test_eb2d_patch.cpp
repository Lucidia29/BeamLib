#include <BeamLib/Element/EulerBernoulli2D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

namespace {

int test_stiffness_matrix_entries() {
    using namespace beamlib;
    const double L = 2.5;
    SectionProperties props;
    props.E  = 200e9;
    props.A  = 0.02;
    props.Iz = 1.0e-5;

    const Vec3 xA(0, 0, 0), xB(L, 0, 0);
    VecN<6> disp = VecN<6>::Zero();
    ElementResult<3> r = EulerBernoulli2D::computeElement(xA, xB, disp, props);

    const double tol = 1e-12;
    const double EA  = props.E * props.A;
    const double EI  = props.E * props.Iz;
    const double L2  = L * L;
    const double L3  = L2 * L;

    auto checkRel = [&](double got, double expected, const char* tag) {
        const double err = std::fabs(got - expected) / std::fabs(expected);
        if (err > tol) {
            std::printf("FAIL [%s]: got %g expected %g (rel %g)\n",
                        tag, got, expected, err);
            return 1;
        }
        return 0;
    };

    if (checkRel(r.ke(0, 0), EA / L,            "K(0,0)=EA/L"))           return 1;
    if (checkRel(r.ke(3, 3), EA / L,            "K(3,3)=EA/L"))           return 1;
    if (checkRel(r.ke(0, 3), -EA / L,           "K(0,3)=-EA/L"))          return 1;
    if (checkRel(r.ke(1, 1), 12.0 * EI / L3,    "K(1,1)=12EI/L^3"))       return 1;
    if (checkRel(r.ke(4, 4), 12.0 * EI / L3,    "K(4,4)=12EI/L^3"))       return 1;
    if (checkRel(r.ke(1, 4), -12.0 * EI / L3,   "K(1,4)=-12EI/L^3"))      return 1;
    if (checkRel(r.ke(2, 2), 4.0 * EI / L,      "K(2,2)=4EI/L"))          return 1;
    if (checkRel(r.ke(5, 5), 4.0 * EI / L,      "K(5,5)=4EI/L"))          return 1;
    if (checkRel(r.ke(2, 5), 2.0 * EI / L,      "K(2,5)=2EI/L"))          return 1;
    if (checkRel(r.ke(1, 2), 6.0 * EI / L2,     "K(1,2)=6EI/L^2"))        return 1;
    if (checkRel(r.ke(1, 5), 6.0 * EI / L2,     "K(1,5)=6EI/L^2"))        return 1;
    if (checkRel(r.ke(2, 4), -6.0 * EI / L2,    "K(2,4)=-6EI/L^2"))       return 1;
    if (checkRel(r.ke(4, 5), -6.0 * EI / L2,    "K(4,5)=-6EI/L^2"))       return 1;

    // Symmetry
    for (int i = 0; i < 6; ++i)
        for (int j = i + 1; j < 6; ++j) {
            const double d = std::fabs(r.ke(i, j) - r.ke(j, i));
            if (d > 1e-14 * (1.0 + std::fabs(r.ke(i, j)))) {
                std::printf("FAIL symmetry K(%d,%d) vs K(%d,%d): %g vs %g\n",
                            i, j, j, i, r.ke(i, j), r.ke(j, i));
                return 1;
            }
        }

    // Zero-displacement residual
    if (r.re.norm() > 1e-12) {
        std::printf("FAIL: zero-disp residual nonzero, norm=%g\n", r.re.norm());
        return 1;
    }
    return 0;
}

int test_rigid_translation_zero_residual() {
    using namespace beamlib;
    const double L = 1.5;
    SectionProperties props;
    props.E  = 200e9;
    props.A  = 0.02;
    props.Iz = 1.0e-5;
    const Vec3 xA(0, 0, 0), xB(L, 0, 0);

    VecN<6> disp;
    disp << 0.123, -0.456, 0.0, 0.123, -0.456, 0.0;
    ElementResult<3> r = EulerBernoulli2D::computeElement(xA, xB, disp, props);
    if (r.re.norm() > 1e-12) {
        std::printf("FAIL rigid translation: re norm %g\n", r.re.norm());
        return 1;
    }
    return 0;
}

int test_constant_axial_strain() {
    using namespace beamlib;
    const int    nElems = 5;
    const double Lbeam  = 2.0;
    const double dx     = Lbeam / nElems;

    BeamModel<EulerBernoulli2D> model;
    model.props.E  = 200e9;
    model.props.A  = 0.01;
    model.props.Iz = 8.333e-6;
    model.nodes.resize(nElems + 1);
    for (int i = 0; i <= nElems; ++i) {
        model.nodes[i].x0 = Vec3(i * dx, 0, 0);
        model.nodes[i].fixed[1] = true;
        model.nodes[i].fixed[2] = true;
    }
    model.nodes[0].fixed[0] = true;
    const double Fx = 5000.0;
    model.nodes[nElems].load[0] = Fx;

    for (int e = 0; e < nElems; ++e) {
        ElementConn c;
        c.nodeA = e;
        c.nodeB = e + 1;
        model.elements.push_back(c);
    }
    model.buildDofMap();
    VecX Fext = model.getExternalForceVector();
    NRResult res = NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(model, Fext);
    if (!res.converged) {
        std::printf("FAIL axial: NR did not converge\n");
        return 1;
    }

    const double EA              = model.props.E * model.props.A;
    const double ux_end_expected = Fx * Lbeam / EA;

    for (int i = 0; i <= nElems; ++i) {
        const double xi             = i * dx;
        const double ux_expected_i  = ux_end_expected * (xi / Lbeam);
        const double ux_i           = model.nodes[i].dof[0];
        const double err            = std::fabs(ux_i - ux_expected_i) /
                                       (std::fabs(ux_end_expected) + 1e-30);
        if (err > 1e-10) {
            std::printf("FAIL axial linearity at node %d: got %g expected %g\n",
                        i, ux_i, ux_expected_i);
            return 1;
        }
    }
    return 0;
}

int test_theta_y_sign_under_Fz() {
    using namespace beamlib;
    const double L = 1.0;
    BeamModel<EulerBernoulli2D> model;
    model.props.E  = 200e9;
    model.props.A  = 0.01;
    model.props.Iz = 8.333e-6;
    model.nodes.resize(2);
    model.nodes[0].x0 = Vec3(0, 0, 0);
    model.nodes[1].x0 = Vec3(L, 0, 0);
    model.nodes[0].fixAll();
    model.nodes[1].load[1] = -1000.0;

    ElementConn c;
    c.nodeA = 0;
    c.nodeB = 1;
    model.elements.push_back(c);

    model.buildDofMap();
    VecX Fext = model.getExternalForceVector();
    NRResult res = NewtonRaphsonSolver<EulerBernoulli2D>::solveOneStep(model, Fext);
    if (!res.converged) {
        std::printf("FAIL theta_y sign: not converged\n");
        return 1;
    }
    const double uz  = model.nodes[1].dof[1];
    const double thy = model.nodes[1].dof[2];
    if (uz >= 0.0) {
        std::printf("FAIL theta_y sign: u_z should be negative, got %g\n", uz);
        return 1;
    }
    if (thy >= 0.0) {
        std::printf("FAIL theta_y sign: theta_y should be negative for downward F_z, got %g\n", thy);
        return 1;
    }
    return 0;
}

} // namespace

int main() {
    if (test_stiffness_matrix_entries())      return 1;
    if (test_rigid_translation_zero_residual()) return 1;
    if (test_constant_axial_strain())         return 1;
    if (test_theta_y_sign_under_Fz())         return 1;
    std::printf("PASS test_eb2d_patch (all sub-tests)\n");
    return 0;
}
