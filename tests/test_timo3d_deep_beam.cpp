#include <BeamLib/Element/EulerBernoulli3D.h>
#include <BeamLib/Element/Timoshenko3D.h>
#include <BeamLib/Model/BeamModel.h>
#include <BeamLib/Solver/NewtonRaphson.h>
#include <cmath>
#include <cstdio>

// Deep-beam test in BOTH bending planes. For a tip-loaded cantilever the
// Timoshenko-to-EB deflection ratio is plane-specific:
//   xz plane (under F_z): ratio_z = 1 + 3 EIy / (kappa_z GA L^2)
//   xy plane (under F_y): ratio_y = 1 + 3 EIz / (kappa_y GA L^2)
//
// Section: rectangular b x h with h = L/4 (L/h = 4, deep). Pick b != h and
// kappa_y != kappa_z so the two planes have visibly different ratios. Both
// numerical ratios must match the analytical ones to relative 1e-9, which
// pins the constitutive constants of each plane.

namespace {

using namespace beamlib;

struct DeepSec {
    double E       = 200e9;
    double nu      = 0.3;
    double G;
    double rho     = 7800.0;
    double L       = 1.0;
    // h along z, b along y. Deep along z (h = L/4); slimmer along y (b = L/6).
    double h, b;
    double A;
    double Iy;   // ∫ z^2 dA = b h^3 / 12
    double Iz;   // ∫ y^2 dA = h b^3 / 12
    double Ix;
    double kappa_y = 5.0 / 6.0;
    double kappa_z = 5.0 / 6.0;
    DeepSec() {
        G = E / (2.0 * (1.0 + nu));
        h = L / 4.0;          // deep about z (L/h = 4)
        b = L / 6.0;          // less deep about y (L/b = 6)
        A  = b * h;
        Iy = b * h * h * h / 12.0;
        Iz = h * b * b * b / 12.0;
        Ix = Iy + Iz;
    }
};

template <typename Elem>
double solveTipDeflection(const DeepSec& s, int nElems, double Fy, double Fz, int comp)
{
    BeamModel<Elem> m;
    m.props.E = s.E; m.props.G = s.G; m.props.rho = s.rho;
    m.props.A = s.A; m.props.Iy = s.Iy; m.props.Iz = s.Iz; m.props.Ix = s.Ix;
    m.props.kappa_y = s.kappa_y; m.props.kappa_z = s.kappa_z;

    const double dx = s.L / nElems;
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
        NewtonRaphsonSolver<Elem>::solveOneStep(m, m.getExternalForceVector());
    if (!res.converged) {
        std::printf("FAIL: NR did not converge for one of the deep-beam solves\n");
        std::exit(1);
    }
    return m.nodes[nElems].dof[comp];
}

} // namespace

int main()
{
    DeepSec s;
    const int nElems = 10;
    const double F = -1000.0;

    // ---- xz plane (F_z only) ----
    const double uz_EB = solveTipDeflection<EulerBernoulli3D>(s, nElems, 0.0, F, /*u_z*/2);
    const double uz_T  = solveTipDeflection<Timoshenko3D>   (s, nElems, 0.0, F, /*u_z*/2);
    const double ratio_z    = uz_T / uz_EB;
    const double ratio_z_an = 1.0 + 3.0 * s.E * s.Iy / (s.kappa_z * s.G * s.A * s.L * s.L);
    const double err_z = std::fabs(ratio_z - ratio_z_an) / ratio_z_an;

    if (std::fabs(uz_T) <= std::fabs(uz_EB)) {
        std::printf("FAIL xz plane: |u_z_Timo|=%g must exceed |u_z_EB|=%g\n",
                    uz_T, uz_EB);
        return 1;
    }
    if (err_z > 1e-9) {
        std::printf("FAIL xz plane: numerical ratio %g vs analytical %g (rel %g)\n",
                    ratio_z, ratio_z_an, err_z);
        return 1;
    }

    // ---- xy plane (F_y only) ----
    const double uy_EB = solveTipDeflection<EulerBernoulli3D>(s, nElems, F, 0.0, /*u_y*/1);
    const double uy_T  = solveTipDeflection<Timoshenko3D>   (s, nElems, F, 0.0, /*u_y*/1);
    const double ratio_y    = uy_T / uy_EB;
    const double ratio_y_an = 1.0 + 3.0 * s.E * s.Iz / (s.kappa_y * s.G * s.A * s.L * s.L);
    const double err_y = std::fabs(ratio_y - ratio_y_an) / ratio_y_an;

    if (std::fabs(uy_T) <= std::fabs(uy_EB)) {
        std::printf("FAIL xy plane: |u_y_Timo|=%g must exceed |u_y_EB|=%g\n",
                    uy_T, uy_EB);
        return 1;
    }
    if (err_y > 1e-9) {
        std::printf("FAIL xy plane: numerical ratio %g vs analytical %g (rel %g)\n",
                    ratio_y, ratio_y_an, err_y);
        return 1;
    }

    // Plane independence: with b != h the two ratios MUST differ. If they
    // happened to be equal it would suggest a single shear term controlling
    // both planes (i.e., a bug).
    if (std::fabs(ratio_z - ratio_y) < 1e-9) {
        std::printf("FAIL: deep-beam ratios coincide (%g, %g) -- expected "
                    "plane-specific values for asymmetric section\n",
                    ratio_z, ratio_y);
        return 1;
    }

    std::printf("PASS test_timo3d_deep_beam:\n"
                "  xz: u_z_EB=%g u_z_Timo=%g ratio=%.6f (analytical %.6f, %.2f%% shear)\n"
                "  xy: u_y_EB=%g u_y_Timo=%g ratio=%.6f (analytical %.6f, %.2f%% shear)\n",
                uz_EB, uz_T, ratio_z, ratio_z_an, 100.0 * (ratio_z - 1.0) / ratio_z,
                uy_EB, uy_T, ratio_y, ratio_y_an, 100.0 * (ratio_y - 1.0) / ratio_y);
    return 0;
}
