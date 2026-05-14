# BeamLib — Codex 执行方案 (Architect-Builder 模式)

## 协作模式

```
每个 Batch 的循环：
  Claude Code (Architect)  →  定义接口 spec + 验收标准 + 数值基准
  Codex (Builder)          →  实现 + 自测 + 迭代直到验收标准通过
  Claude Code (Reviewer)   →  Review 代码正确性，重点检查数值精度和边界条件
```

**协议层**: Git — Codex 在 feature branch 上工作，完成后 Claude review。
**项目规范**: 见 `CLAUDE.md`（Codex 必读）。

---

## Batch 1: 项目脚手架 + Core 类型

### 你的任务
创建完整的项目骨架，让后续代码可以 include 和编译。

### 1.1 创建目录和 CMakeLists.txt

**输出**:
```
BeamLib/
├── CMakeLists.txt
├── include/BeamLib/{Core,Math,Element,Model,Solver,Load}/
├── src/       (空，后续放 GeomExact3D.cpp)
└── tests/     (空，后续放测试)
```

**CMakeLists.txt 规格**:
- `cmake_minimum_required(VERSION 3.16)`, `project(BeamLib LANGUAGES CXX)`, C++17
- Eigen3 查找：支持 `-DEIGEN3_INCLUDE_DIR=<path>` 和 `find_package(Eigen3)`
- 参考 `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\CMakeLists.txt` 的 Eigen 查找逻辑，照搬即可
- 定义 `beamlib` 静态库 target（源文件列表先留空或放一个占位 .cpp）
- `enable_testing()`

### 1.2 Core 类型文件

创建以下三个文件，内容**严格按照下面给出的代码**，不做修改：

**`include/BeamLib/Core/Types.h`** — Eigen 别名 + 稀疏矩阵类型:
```cpp
#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>

namespace beamlib {

using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Mat2 = Eigen::Matrix2d;
using Mat3 = Eigen::Matrix3d;
using VecX = Eigen::VectorXd;
using MatX = Eigen::MatrixXd;
using SpMat = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

template <int N>
using VecN = Eigen::Matrix<double, N, 1>;

template <int M, int N>
using MatMN = Eigen::Matrix<double, M, N>;

using Tensor4_3x3 = std::array<std::array<Mat3, 3>, 3>;

} // namespace beamlib
```

**`include/BeamLib/Core/Node.h`** — 模板化节点:
```cpp
#pragma once
#include "Types.h"

namespace beamlib {

template <int NDofsPerNode>
struct Node {
    static constexpr int nDofs = NDofsPerNode;
    Vec3 x0 = Vec3::Zero();
    VecN<NDofsPerNode> dof = VecN<NDofsPerNode>::Zero();
    std::array<bool, NDofsPerNode> fixed = {};
    VecN<NDofsPerNode> load = VecN<NDofsPerNode>::Zero();

    Node() = default;
    explicit Node(const Vec3& x0_) : x0(x0_) {}
    void fixAll() { fixed.fill(true); }
};

using Node2D = Node<3>;
using Node3D = Node<6>;

} // namespace beamlib
```

**`include/BeamLib/Core/SectionProperties.h`** — 截面属性:
```cpp
#pragma once
#include "Types.h"

namespace beamlib {

struct SectionProperties {
    double E = 0.0;
    double G = 0.0;
    double rho = 0.0;
    double A = 0.0;
    double Iy = 0.0;
    double Iz = 0.0;
    double Ix = 0.0;
    double kappa = 5.0 / 6.0;

    Mat3 C_axialShear() const {
        Mat3 C = Mat3::Zero();
        C(0, 0) = E * A;
        C(1, 1) = kappa * G * A;
        C(2, 2) = kappa * G * A;
        return C;
    }

    Mat3 D_bendTorsion() const {
        Mat3 D = Mat3::Zero();
        D(0, 0) = G * Ix;
        D(1, 1) = E * Iy;
        D(2, 2) = E * Iz;
        return D;
    }
};

} // namespace beamlib
```

### 1.3 ElementBase + BeamMath3D 迁移

**`include/BeamLib/Element/ElementBase.h`** — 元素结果模板:
```cpp
#pragma once
#include "../Core/Types.h"

namespace beamlib {

template <int NDofsPerNode>
struct ElementResult {
    static constexpr int elemDofs = 2 * NDofsPerNode;
    VecN<elemDofs> re = VecN<elemDofs>::Zero();
    MatMN<elemDofs, elemDofs> ke = MatMN<elemDofs, elemDofs>::Zero();
};

template <int NDofsPerNode>
struct ElementMassResult {
    static constexpr int elemDofs = 2 * NDofsPerNode;
    MatMN<elemDofs, elemDofs> me = MatMN<elemDofs, elemDofs>::Zero();
};

struct ElementConn {
    int nodeA = 0;
    int nodeB = 0;
};

} // namespace beamlib
```

**`include/BeamLib/Math/BeamMath3D.h`** — 从现有代码迁移:
- 源文件: `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\include\BeamMath.h`
- 操作: 复制全部内容，仅将 `namespace beam_fem` 改为 `namespace beamlib`
- **不修改任何算法代码** — 这是经过 14 个测试验证的代码

### 验收标准
- [ ] `cmake -B build -DEIGEN3_INCLUDE_DIR=<path>` 成功
- [ ] 写一个 `tests/test_compile_check.cpp` include 所有头文件，编译通过
- [ ] `Node<3>` 和 `Node<6>` 可实例化，`fixAll()` 工作正常

### Claude Review 重点
- BeamMath3D.h 与原文件逐行 diff，确认仅 namespace 变更
- Types.h 中的模板别名是否正确支持固定大小矩阵

---

## Batch 2: 核心 Pipeline — GeomExact3D + 装配 + 求解

### 你的任务
实现 3D 几何精确梁元素、稀疏装配模型、Newton-Raphson 求解器，并用一个悬臂梁测试验证端到端正确性。

### 2.1 GeomExact3D 元素

**输出文件**: `include/BeamLib/Element/GeomExact3D.h` + `src/GeomExact3D.cpp`

**接口** (GeomExact3D.h):
```cpp
namespace beamlib {
struct GeomExact3D {
    static constexpr int nDofsPerNode = 6;
    static constexpr bool hasTransformation = true;
    using Result = ElementResult<6>;

    static Result computeElement(
        const Vec3& xA, const Vec3& xB,
        const VecN<12>& dispVec,
        const SectionProperties& props);

    static MatMN<12, 12> computeTransformation(
        const Vec3& xA, const Vec3& xB);
};
}
```

**实现** (GeomExact3D.cpp):
- `computeElement`: 从 `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\src\BeamElement3D.cpp` 复制 `computeElement` 和 `computeResidual`、`computeStiffBlock` 两个 static 函数
- 将 `(const Mat3& C, const Mat3& D)` 参数改为 `(const SectionProperties& props)`，在函数开头 `Mat3 C = props.C_axialShear(); Mat3 D = props.D_bendTorsion();`
- 其余算法零改动
- `computeTransformation`: 从 `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\src\BeamModel.cpp` 的 `BeamModel::calcTransformationMatrix` 复制

### 2.2 DofMap

**输出文件**: `include/BeamLib/Model/DofMap.h`

**接口**: header-only 模板类 `DofMap<NDofsPerNode>`
- `void build(const std::vector<Node<NDofsPerNode>>& nodes)` — 遍历节点，为非固定 DOF 分配连续编号
- `int numFree() const`
- `int freeDofIndex(int nodeIdx, int localDof) const` — 固定 DOF 返回 -1

参考: `BeamModel::buildDofMap()` 和 `BeamModel::freeDofIndex()` 的逻辑。

### 2.3 BeamModel 稀疏装配

**输出文件**: `include/BeamLib/Model/BeamModel.h`

**接口**: 模板类 `BeamModel<ElemType>`

```cpp
template <typename ElemType>
class BeamModel {
public:
    static constexpr int nDPN = ElemType::nDofsPerNode;
    static constexpr int elemDofs = 2 * nDPN;

    std::vector<Node<nDPN>> nodes;
    std::vector<ElementConn> elements;
    SectionProperties props;

    void buildDofMap();
    int numFreeDofs() const;
    int freeDofIndex(int nodeIdx, int localDof) const;

    void assemble(VecX& R, SpMat& K) const;
    void assembleResidualOnly(VecX& R) const;
    VecX getExternalForceVector() const;
    void scatterFreeDofs(const VecX& delta);
    VecX gatherFreeDofs() const;
};
```

**装配逻辑** (assemble 方法):
1. 预分配 `std::vector<Triplet>`，reserve `elements.size() * elemDofs * elemDofs`
2. 对每个 element:
   - 收集两节点 DOF → `dispVec`
   - 用 `if constexpr (ElemType::hasTransformation)` 判断是否需要坐标变换
   - 若需要变换: `T = ElemType::computeTransformation(xA, xB)`
     - `dispLocal = T * dispGlobal`
     - `xA_local = lambda * xA`, `xB_local = lambda * xB` (lambda = T 左上 3×3)
     - `result = ElemType::computeElement(xA_local, xB_local, dispLocal, props)`
     - `re_g = T^T * re_l`, `ke_g = T^T * ke_l * T`
   - 否则直接调用
   - 构建 location map，用 Triplet 收集
3. `K.setFromTriplets(...)`

参考: `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\src\BeamModel.cpp`，但改用稀疏矩阵。

### 2.4 NewtonRaphson 求解器

**输出文件**: `include/BeamLib/Solver/LinearSolver.h` + `include/BeamLib/Solver/NewtonRaphson.h`

LinearSolver: `Eigen::SparseLU` 包装，static 方法 `VecX solve(const SpMat& K, const VecX& rhs)`。

NewtonRaphson: 模板类 `NewtonRaphsonSolver<ElemType>`
- `NRResult solveOneStep(BeamModel<ElemType>& model, const VecX& Fext)`
  - 循环: assemble → R -= Fext → 检查 ||R|| < tol → solve → scatter
  - 与参考 `BeamSolver::solveOneStep` 逻辑一致
- `std::vector<StepResult> solveStatic(BeamModel<ElemType>& model, int nSteps, int tipNodeIdx)`
  - 载荷步分步：Fn = n/nSteps * f_total
  - 与参考 `BeamSolver::solveStatic` 逻辑一致

### 2.5 端到端测试: GE3D 悬臂梁

**输出文件**: `tests/test_ge3d_cantilever.cpp`

**设置**:
```
E=1e6, G=5e5, kappa=1.0, A=1, Iy=Iz=Ix=1
L=1.0, 10 等长单元, 11 节点沿 x 轴
节点 0 全固定, 节点 10 施加 F_z = 1.0
1 步 Newton-Raphson, tol=1e-7
```

**数值基准**:
- 末端 u_z 解析解: `FL³/(3EI) = 1.0/(3×1e6) = 3.3333e-7`
- 末端 θ_y 解析解: `FL²/(2EI) = 1.0/(2×1e6) = 5.0e-7`
- 容差: 相对误差 < 1e-4（GE 梁在小力下接近线性解）

**程序退出**: return 0 全通过, return 1 有失败。打印每个检查点的 computed vs reference 值。

**CMakeLists.txt 更新**: 添加 `test_ge3d_cantilever` target + `add_test`。

### 验收标准
- [ ] `cmake --build build` 编译通过（零 warning 最好）
- [ ] Newton-Raphson ≤ 5 次迭代收敛
- [ ] u_z 相对误差 < 1e-4
- [ ] θ_y 相对误差 < 1e-4
- [ ] `ctest --test-dir build` 全部 PASS

### Claude Review 重点
- `computeElement` 算法代码与原文件 diff — **必须零算法改动**
- `assemble` 中坐标变换逻辑: `dispLocal = T * dispGlobal`，`re_g = T^T * re_l`，`ke_g = T^T * ke_l * T` 顺序是否正确
- 稀疏 vs 稠密结果是否一致
- Newton-Raphson 残差计算: `R_int - F_ext`，符号约定是否与参考一致（参考中是 `R -= Fn`，即内力残差减去外力）

---

## Batch 3: 大变形验证 + 显式动力学

### 你的任务
验证 GE3D 的大变形能力，并实现显式时间积分器。

### 3.1 大变形测试

**输出文件**: `tests/test_ge3d_large_deformation.cpp`

**设置**:
```
E=1e6, G=5e5, kappa=1.0, A=1, Iy=Iz=Ix=1
L=1.0, 10 单元
节点 0 固定, 节点 10 施加 F_z = 5e5 (大力，使挠度≈0.5L)
20 步载荷步, tol=1e-7
```

**验证**:
1. 所有 20 个载荷步收敛
2. 无 NaN/Inf
3. 末端位移合理（u_z 大于 0，u_x 因大变形应非零）
4. 最终步残差 < 1e-6

### 3.2 显式中心差分积分器

**输出文件**: `include/BeamLib/Solver/TimeIntegrator.h` + `include/BeamLib/Solver/ExplicitCentral.h`

**接口**: 模板类 `ExplicitCentralDifference<ElemType>`

```cpp
template <typename ElemType>
class ExplicitCentralDifference {
public:
    double alpha = 0.0;  // Rayleigh 阻尼

    ExplicitCentralDifference(BeamModel<ElemType>& model);
    void step(double dt);
    BeamEnergy computeEnergy() const;
    double estimateCriticalTimestep() const;

    std::vector<VecN<ElemType::nDofsPerNode>> velocity;
    std::vector<VecN<ElemType::nDofsPerNode>> acceleration;
    int stepCount() const;
};
```

**实现参考**: `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\src\ExplicitBeamIntegrator.cpp`

关键算法:
- **集中质量**: `m_trans = props.rho * props.A * h * 0.5`, `J_rot = props.rho * (props.Iy + props.Iz) * h * 0.5`
- **velocity-Verlet**: 半步速度 → 全步位移 → 计算加速度 → 完成速度
- **加速度**: `a = (F_ext - R_int - alpha * M * v) / M`（逐 DOF）
- **临界时间步**: `dt_crit = h_min / sqrt(E / rho)`
- **应变能**: 参照原代码，需要在局部坐标系下调用 BeamMath3D 的 strain 函数

### 3.3 显式动力学测试

**输出文件**: `tests/test_explicit_free_vibration.cpp`

**设置**:
```
E=1e6, G=5e5, rho=1000, kappa=1.0
A=0.01, Iy=Iz=1e-4, Ix=2e-4
L=1.0, 10 单元
节点 0 固定
初始条件: 节点 10 的 dof(2) = 0.001 (初始 u_z 位移)
dt = 0.5 * dt_critical, 运行 1000 步
```

**验证**:
- 1000 步无 NaN/Inf
- 无阻尼: 总能量(动能+应变能)波动 < 1%

### 验收标准
- [ ] 大变形: 20 步全收敛，无 NaN/Inf
- [ ] 显式: 1000 步无 NaN/Inf
- [ ] 能量守恒: `(E_max - E_min) / E_mean < 0.01`
- [ ] `ctest` 全部 PASS

### Claude Review 重点
- 大变形结果是否物理合理（大力下 u_x 应因缩短效应为负）
- 显式积分器的 velocity-Verlet 半步逻辑是否正确
- 集中质量计算是否用了 props 而非硬编码
- 应变能计算中的坐标变换是否与 assemble 一致

---

## Batch 4: 2D Euler-Bernoulli 梁

### 你的任务
实现 2D EB 元素并用解析解验证。

### 4.1 Rotation2D + EulerBernoulli2D 元素

**输出文件**: `include/BeamLib/Math/Rotation2D.h` + `include/BeamLib/Element/EulerBernoulli2D.h`

**接口**: struct `EulerBernoulli2D` with `nDofsPerNode = 3`, `hasTransformation = true`

DOF 顺序 (局部坐标): `[u_x, u_z, θ_y]` per node

**刚度矩阵** (来自 MATLAB `EulerBernoulliBeam2DElement.m`):
- 轴向: `K(0,0)=K(3,3)=EA/L`, `K(0,3)=K(3,0)=-EA/L`
- 弯曲:
  ```
  K(1,1)=K(4,4) = 12EI/L³
  K(1,2)=K(2,1) = 6EI/L²
  K(1,4)=K(4,1) = -12EI/L³
  K(1,5)=K(5,1) = 6EI/L²
  K(2,2)=K(5,5) = 4EI/L
  K(2,4)=K(4,2) = -6EI/L²
  K(2,5)=K(5,2) = 2EI/L
  K(4,5)=K(5,4) = -6EI/L²
  ```
  其中 `I = props.Iz`

**残差**: 线性问题 `re = -K * dispVec`

**质量矩阵**: 一致质量，系数 `rho*A/420`:
```
M(0,0)=M(3,3)=140L, M(0,3)=M(3,0)=70L
M(1,1)=M(4,4)=156L, M(1,2)=M(2,1)=22L², M(1,4)=M(4,1)=54L
M(1,5)=M(5,1)=-13L², M(2,2)=M(5,5)=4L³, M(2,4)=M(4,2)=13L²
M(2,5)=M(5,2)=-3L³, M(4,5)=M(5,4)=-22L²
```

**2D 变换矩阵** (`Rotation2D.h`):
```
T = [c  -s  0  0   0  0]
    [s   c  0  0   0  0]
    [0   0  1  0   0  0]
    [0   0  0  c  -s  0]
    [0   0  0  s   c  0]
    [0   0  0  0   0  1]
```
其中 `c = dx_x/L`, `s = dx_z/L` (`dx_z` 对应 xz 平面)

### 4.2 测试: 2D EB 悬臂梁

**输出文件**: `tests/test_eb2d_cantilever.cpp`

**设置**:
```
E=200e9, A=0.01, Iz=8.333e-6, rho=7800
L=1.0, 10 单元
节点 0 全固定 (3 DOFs), 节点 10 施加 F_z = -1000
1 步 Newton-Raphson, tol=1e-7
```

**数值基准** (解析解精确匹配):
- `u_z = -FL³/(3EI) = -1000/(3×200e9×8.333e-6) = -2.0e-4`
- `θ_y = -FL²/(2EI) = -1000/(2×200e9×8.333e-6) = -3.0e-4`
- 容差: 相对误差 < **1e-10**（线性问题必须精确）

### 验收标准
- [ ] 编译通过
- [ ] 沿 x 轴的梁变换矩阵为单位矩阵
- [ ] 单元素 K(0,0)=EA/L, K(1,1)=12EI/L³ 正确
- [ ] **1 次** NR 迭代收敛（线性问题）
- [ ] u_z, θ_y 与解析解相对误差 < 1e-10
- [ ] `ctest` PASS

### Claude Review 重点
- DOF 顺序: MATLAB 用 (x,y) 平面，我们用 (x,z) 平面 — 检查变换矩阵 sin/cos 取的是 z 分量而非 y
- 刚度矩阵对称性
- 线性问题必须 1 次收敛 — 如果不是，说明残差符号或装配有 bug

---

## Batch 5: 3D Euler-Bernoulli 梁

### 你的任务
实现 3D EB 元素并用解析解验证。

### 5.1 EulerBernoulli3D 元素

**输出文件**: `include/BeamLib/Element/EulerBernoulli3D.h`

**接口**: struct `EulerBernoulli3D` with `nDofsPerNode = 6`, `hasTransformation = true`

DOF 顺序 (局部, 0-based): `[u_x(0), u_y(1), u_z(2), θ_x(3), θ_y(4), θ_z(5)]` per node

**刚度矩阵** (来自 MATLAB `EulerBernoulliBeam3DElement.m`, 注意 MATLAB 1-based → C++ 0-based):
- 轴向 (DOF 0,6): `EA/L`
- 扭转 (DOF 3,9): `GIx/L`
- xy 平面弯曲 (DOF 1,5,7,11): 用 `EIz`
  ```
  K(1,1)=K(7,7)=12EIz/L³,  K(1,7)=K(7,1)=-12EIz/L³
  K(1,5)=K(5,1)=6EIz/L²,   K(1,11)=K(11,1)=6EIz/L²
  K(5,7)=K(7,5)=-6EIz/L²,  K(7,11)=K(11,7)=-6EIz/L²
  K(5,5)=K(11,11)=4EIz/L,  K(5,11)=K(11,5)=2EIz/L
  ```
- xz 平面弯曲 (DOF 2,4,8,10): 用 `EIy`，同样的模式

**变换矩阵**: 复用 `GeomExact3D::computeTransformation`

**残差**: `re = -K * dispVec`

**质量矩阵**: 一致质量，直接从 MATLAB 代码转写（包括轴向、扭转、两个弯曲平面）

### 5.2 测试: 3D EB 悬臂梁

**输出文件**: `tests/test_eb3d_cantilever.cpp`

**设置**:
```
E=200e9, G=80e9, rho=7800
A=0.01, Iy=Iz=8.333e-6, Ix=16.667e-6
L=1.0, 10 单元沿 x 轴
节点 0 全固定, 节点 10 施加 F_z = -1000
```

**数值基准**:
- `u_z = -FL³/(3EIy) = -2.0e-4`, `θ_y = -FL²/(2EIy) = -3.0e-4`
- 其余 DOF ≈ 0 (绝对值 < 1e-15)
- 容差: 相对误差 < **1e-10**

### 验收标准
- [ ] 编译通过
- [ ] 单元素: K(0,0)=EA/L, K(3,3)=GIx/L, K(1,1)=12EIz/L³, K(2,2)=12EIy/L³
- [ ] 1 次 NR 收敛
- [ ] u_z, θ_y 精度 < 1e-10
- [ ] 非加载方向 DOF ≈ 0

### Claude Review 重点
- DOF 索引 0-based 与 MATLAB 1-based 的映射是否正确（这是最常见的 off-by-one bug 来源）
- 3D EB 刚度矩阵必须对称
- 质量矩阵中 bending-in-xz 的耦合 DOF 索引

---

## Batch 6: 载荷管理 + 最终集成

### 你的任务
实现 LoadManager，确保所有测试通过。

### 6.1 LoadTypes + LoadManager

**输出文件**: `include/BeamLib/Load/LoadTypes.h` + `include/BeamLib/Load/LoadManager.h`

**接口**: `LoadManager<ElemType>` 模板类
- `addPointLoad(nodeId, dofIndex, value)`
- `addPrescribedDisplacement(nodeId, dofIndex, value)` — 设置 `node.fixed[dofIndex] = true` 并记录约束值
- `applyToModel(model, loadFactor)` — 将点载荷乘以 loadFactor 写入 node.load
- `clear()`

### 6.2 最终验证

确保所有已有测试仍然 PASS: `ctest --test-dir build`

### 验收标准
- [ ] 编译通过
- [ ] `applyToModel` 正确设置 node.load（loadFactor 正确缩放）
- [ ] 全部已有测试 PASS
- [ ] CMakeLists.txt 中所有 test target 均已添加

---

## 执行顺序和依赖

```
Batch 1 (脚手架)    → Codex 可直接开始
Batch 2 (核心管线)  → 依赖 Batch 1
Batch 3 (大变形+动力学) → 依赖 Batch 2
Batch 4 (2D EB)     → 依赖 Batch 2
Batch 5 (3D EB)     → 依赖 Batch 2
Batch 6 (载荷管理)  → 依赖 Batch 2
```

Batch 3/4/5/6 之间无相互依赖，可按任意顺序或并行执行。

## Git 分支策略

```
main                    ← 仅 Claude review 通过后合并
├── feat/batch-1-scaffold
├── feat/batch-2-core-pipeline
├── feat/batch-3-dynamics
├── feat/batch-4-eb2d
├── feat/batch-5-eb3d
└── feat/batch-6-load-manager
```

每个 Batch 完成后提交到对应 branch，Claude review 后合入 main。
