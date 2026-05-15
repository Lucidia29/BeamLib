#include <BeamLib/Core/Types.h>
#include <BeamLib/Core/Node.h>
#include <BeamLib/Core/SectionProperties.h>
#include <BeamLib/Element/ElementBase.h>
#include <cstdio>

int main() {
    // Types: fixed-size aliases
    beamlib::Vec3 v = beamlib::Vec3::Zero();
    beamlib::Mat3 m = beamlib::Mat3::Identity();
    beamlib::VecX vx(10);
    vx.setZero();
    beamlib::SpMat sp(3, 3);
    beamlib::VecN<6> v6 = beamlib::VecN<6>::Zero();
    beamlib::MatMN<6, 6> m66 = beamlib::MatMN<6, 6>::Zero();

    // Node2D
    beamlib::Node2D n2d(beamlib::Vec3(1.0, 0.0, 2.0));
    n2d.fixAll();
    for (int i = 0; i < 3; ++i) {
        if (!n2d.fixed[i]) {
            std::printf("FAIL: Node2D.fixAll() did not set fixed[%d]\n", i);
            return 1;
        }
    }
    n2d.load[2] = -1000.0;

    // Node3D
    beamlib::Node3D n3d;
    n3d.x0 = beamlib::Vec3(0.0, 0.0, 0.0);
    n3d.fixed[0] = true;
    n3d.fixed[1] = true;
    n3d.fixed[2] = true;

    // SectionProperties
    beamlib::SectionProperties sec;
    sec.E = 200e9;
    sec.G = 80e9;
    sec.A = 0.01;
    sec.Iy = 8.333e-6;
    sec.Iz = 8.333e-6;
    sec.Ix = 1.0e-5;
    sec.rho = 7800.0;

    beamlib::Mat3 C = sec.C_axialShear();
    if (C(0, 0) != sec.E * sec.A) {
        std::printf("FAIL: C_axialShear(0,0) mismatch\n");
        return 1;
    }
    if (C(1, 1) != sec.kappa_y * sec.G * sec.A) {
        std::printf("FAIL: C_axialShear(1,1) mismatch\n");
        return 1;
    }

    beamlib::Mat3 D = sec.D_bendTorsion();
    if (D(0, 0) != sec.G * sec.Ix) {
        std::printf("FAIL: D_bendTorsion(0,0) mismatch\n");
        return 1;
    }
    if (D(1, 1) != sec.E * sec.Iy) {
        std::printf("FAIL: D_bendTorsion(1,1) mismatch\n");
        return 1;
    }

    if (sec.kappa() != sec.kappa_z) {
        std::printf("FAIL: kappa() != kappa_z\n");
        return 1;
    }

    // ElementResult
    beamlib::ElementResult<3> res2d;
    if (res2d.re.size() != 6 || res2d.ke.rows() != 6 || res2d.ke.cols() != 6) {
        std::printf("FAIL: ElementResult<3> dimensions\n");
        return 1;
    }

    beamlib::ElementResult<6> res3d;
    if (res3d.re.size() != 12 || res3d.ke.rows() != 12 || res3d.ke.cols() != 12) {
        std::printf("FAIL: ElementResult<6> dimensions\n");
        return 1;
    }

    // ElementMassResult
    beamlib::ElementMassResult<3> mass2d;
    if (mass2d.me.rows() != 6 || mass2d.me.cols() != 6) {
        std::printf("FAIL: ElementMassResult<3> dimensions\n");
        return 1;
    }

    // ElementConn
    beamlib::ElementConn conn;
    conn.nodeA = 0;
    conn.nodeB = 1;
    if (conn.refVector != beamlib::Vec3(0, 1, 0)) {
        std::printf("FAIL: ElementConn default refVector\n");
        return 1;
    }

    std::printf("All compile checks passed.\n");
    return 0;
}
