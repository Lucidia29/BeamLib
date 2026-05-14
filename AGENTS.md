# BeamLib — 梁单元有限元库

## 项目概述
C++17 + Eigen3 的独立梁单元有限元库，覆盖 2D/3D 小变形（Euler-Bernoulli）和大变形/大转动（几何精确梁）。

## 硬约束
- **namespace**: 所有代码在 `beamlib` 命名空间内
- **C++ 标准**: C++17，不使用 C++20 特性
- **依赖**: 仅 Eigen3（header-only），无其他第三方库
- **构建**: CMake ≥ 3.16
- **测试**: 纯 main() 函数 + ctest，不用 gtest/catch2

## 目录约定
```
include/BeamLib/Core/       # 基础类型 (Types.h, Node.h, SectionProperties.h)
include/BeamLib/Math/       # 数学工具 (BeamMath3D.h, Rotation2D.h)
include/BeamLib/Element/    # 单元 (ElementBase.h, 各元素类)
include/BeamLib/Model/      # 模型和装配 (DofMap.h, BeamModel.h)
include/BeamLib/Solver/     # 求解器 (LinearSolver.h, NewtonRaphson.h, 时间积分)
include/BeamLib/Load/       # 载荷 (LoadTypes.h, LoadManager.h)
src/                        # 仅 GeomExact3D.cpp（其余 header-only）
tests/                      # 测试可执行文件
```

## 元素类型接口约定
每个元素类型是一个 struct，提供 static constexpr 和 static 方法：
```cpp
struct SomeElement {
    static constexpr int nDofsPerNode = N;       // 必须
    static constexpr bool hasTransformation = true/false; // 必须
    
    static ElementResult<N> computeElement(       // 必须
        const Vec3& xA, const Vec3& xB,
        const VecN<2*N>& dispVec,
        const SectionProperties& props);
    
    static MatMN<2*N, 2*N> computeTransformation( // 当 hasTransformation=true 时必须
        const Vec3& xA, const Vec3& xB);
    
    static ElementMassResult<N> computeMass(      // 可选，动力学需要
        const Vec3& xA, const Vec3& xB,
        const SectionProperties& props);
};
```

## DOF 顺序约定
- **2D (nDofsPerNode=3)**: [u_x, u_z, θ_y] — xz 平面内
- **3D (nDofsPerNode=6)**: [u_x, u_y, u_z, θ_x, θ_y, θ_z]
- **单元 DOF**: [nodeA DOFs, nodeB DOFs]

## 坐标约定
- 梁的参考方向为局部 x 轴
- 2D 问题在 xz 平面内（x 水平，z 竖直）
- 3D 变换矩阵默认用参考点 (0,1,0) 构建局部 y 轴，当梁平行于全局 y 轴时切换到 (0,0,1)

## 数值精度要求
- 线性问题（EB 梁）: 1 次 NR 迭代收敛，与解析解相对误差 < 1e-10
- 非线性问题（GE 梁小变形）: 与解析解相对误差 < 1e-4
- 非线性问题（GE 梁大变形）: NR 收敛，无 NaN/Inf
- 显式动力学: 无阻尼能量守恒误差 < 1%

## 参考代码（只读，不修改）
- 已有 C++: `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\`
- MATLAB: `C:\_ZW\DEM\01_res_and_dev\FSI_dev\JaphyBeam\`

## 代码风格
- 不写注释，除非解释 WHY（非 WHAT）
- 不加 error handling / fallback，除非在系统边界
- 优先 header-only，只有复杂实现才用 .cpp
- 变量命名遵循物理约定（E, G, rho, Iy, Iz 等）
