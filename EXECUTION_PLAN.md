# BeamLib — 执行总计划 (v3)

## 文档角色

- `PROJECT_SPEC.md`: 项目总规范和硬约束（以此为准）
- `EXECUTION_PLAN.md`: 执行总计划、架构决策、实现步骤、验收标准
- `REVIEW_NOTES.md`: 审阅意见、风险点、计划修订建议

## 总体策略

### 理论先行

每种梁理论的实现流程:
```
理论推导 (LaTeX PDF) → Codex/GPT Review → 编码 → 独立验证 → Codex Review
```

不跳过理论文档直接编码。MATLAB/C++ 参考代码辅助理解，但公式以第一性原理推导为准。

### 从简到繁

```
EB2D (最简单，建立核心 pipeline)
  → EB2D 后处理补全
    → EB3D (扩展到 3D)
      → Timo2D/3D (加入剪切)
        → GE3D (非线性，迁移现有代码)
          → GE2D (独立 2D 推导)
            → 动力学 + 模态
              → 最终集成 + API 封装
```

### 协作循环

```
每个 Batch:
  CC: 理论文档 + 接口定义 + 实现 + 自测
  GPT/Codex: Review 理论、数值、边界条件
  若有修订 → 写入 REVIEW_NOTES.md → CC 修正 → re-review
```

---

## Batch 1: 项目脚手架 + Core 类型

### 目标
创建完整项目骨架，所有头文件可 include 和编译。

### 1.1 目录和 CMakeLists.txt

**输出**:
```
BeamLib/
├── CMakeLists.txt
├── include/BeamLib/{Core,Math,Element,Model,Solver,Load,PostProcess,API}/
├── src/
├── tests/
└── docs/{theory,tutorials}/
```

**CMakeLists.txt 规格**:
- `cmake_minimum_required(VERSION 3.16)`, `project(BeamLib LANGUAGES CXX)`, C++17
- Eigen3 查找：支持 `-DEIGEN3_INCLUDE_DIR=<path>` 和 `find_package(Eigen3)`
- 参考 `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\CMakeLists.txt` 的 Eigen 查找逻辑
- 定义 `beamlib` 为 INTERFACE library target（纯 header-only 阶段）
- Batch 6 引入 GeomExact3D.cpp 后改为 static library
- `enable_testing()`

### 1.2 Core 类型文件

**`include/BeamLib/Core/Types.h`**:
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

**`include/BeamLib/Core/Node.h`**:
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

**`include/BeamLib/Core/SectionProperties.h`**:
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
    double kappa_y = 5.0 / 6.0;
    double kappa_z = 5.0 / 6.0;

    Mat3 C_axialShear() const {
        Mat3 C = Mat3::Zero();
        C(0, 0) = E * A;
        C(1, 1) = kappa_y * G * A;
        C(2, 2) = kappa_z * G * A;
        return C;
    }

    Mat3 D_bendTorsion() const {
        Mat3 D = Mat3::Zero();
        D(0, 0) = G * Ix;
        D(1, 1) = E * Iy;
        D(2, 2) = E * Iz;
        return D;
    }

    double kappa() const { return kappa_z; }
};

} // namespace beamlib
```

### 1.3 ElementBase

**`include/BeamLib/Element/ElementBase.h`**:
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
    Vec3 refVector = Vec3(0, 1, 0);
};

} // namespace beamlib
```

### 验收标准
- [ ] `cmake -B build -DEIGEN3_INCLUDE_DIR=<path>` 成功
- [ ] `tests/test_compile_check.cpp` include 所有头文件，编译通过
- [ ] `Node<3>` 和 `Node<6>` 可实例化，`fixAll()` 正确
- [ ] INTERFACE library target 正确传递 include path

### Review 重点
- Types.h 模板别名是否正确支持固定大小矩阵运算
- INTERFACE vs static library 的 CMake 行为
- `kappa_y`/`kappa_z` 与 `C_axialShear()` 的一致性

---

## Batch 2A: EB2D 核心 Pipeline

### 目标
实现 2D Euler-Bernoulli 梁 + DofMap + BeamModel 装配 + 线性/NR 求解，验证位移正确性。**只做位移求解，不做后处理。**

### 前置: EB2D 理论文档

**输出**: `docs/theory/01_euler_bernoulli_2d.tex` → PDF

**内容大纲**:
1. 运动学假设: 平截面、法线不变
2. 坐标系: xz 平面，DOF = [u_x, u_z, θ_y]
3. 应变-位移: ε_x = du/dx - z·d²w/dx²
4. 虚功原理 → 弱形式
5. 形函数: 线性 (axial) + Hermite 三次 (bending)
6. 单元刚度矩阵 K^e 完整推导 (6×6)
7. 一致质量矩阵 M^e 推导 (6×6)
8. 2D 坐标变换 T (全局 xz → 局部)
9. 全局刚度: K_g = T^T K_l T
10. 验证设计: 悬臂梁、简支梁解析解

**必须包含附表**: 符号约定表、公式→C++ 索引映射表、验证容差表。

### 2A.1 Rotation2D

**输出**: `include/BeamLib/Math/Rotation2D.h`

2D xz 平面变换矩阵:
```
T = [c  -s  0 | 0   0  0]     c = dx_x / L
    [s   c  0 | 0   0  0]     s = dx_z / L
    [0   0  1 | 0   0  0]     L = ||xB - xA||
    [-----------+---------]
    [0   0  0 | c  -s  0]
    [0   0  0 | s   c  0]
    [0   0  0 | 0   0  1]
```

提供 static 函数 `MatMN<6,6> compute(const Vec3& xA, const Vec3& xB)`。

### 2A.2 EulerBernoulli2D 单元

**输出**: `include/BeamLib/Element/EulerBernoulli2D.h`

```cpp
struct EulerBernoulli2D {
    static constexpr int nDofsPerNode = 3;
    static constexpr bool hasTransformation = true;
    static constexpr bool isLinear = true;

    static ElementResult<3> computeElement(
        const Vec3& xA, const Vec3& xB,
        const VecN<6>& dispVec,
        const SectionProperties& props);

    static MatMN<6, 6> computeTransformation(
        const Vec3& xA, const Vec3& xB,
        const Vec3& refVector);   // 2D 忽略 refVector

    static ElementMassResult<3> computeMass(
        const Vec3& xA, const Vec3& xB,
        const SectionProperties& props);
};
```

**刚度矩阵** (局部坐标, I = props.Iz):
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

**残差**: `re = ke * dispVec` (内力向量)

**一致质量矩阵**: 系数 `ρA/420`。需从理论文档独立推导验证。

### 2A.3 DofMap

**输出**: `include/BeamLib/Model/DofMap.h`

Header-only 模板类 `DofMap<NDofsPerNode>`:
- `void build(const std::vector<Node<NDofsPerNode>>& nodes)` — 为非固定 DOF 分配连续编号
- `int numFree() const`
- `int freeDofIndex(int nodeIdx, int localDof) const` — 固定 DOF 返回 -1

### 2A.4 BeamModel 稀疏装配

**输出**: `include/BeamLib/Model/BeamModel.h`

模板类 `BeamModel<ElemType>`:
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
    void assembleStiffnessOnly(SpMat& K) const;
    void assembleMass(SpMat& M) const;       // 接口预留，Batch 2B 实现
    VecX getExternalForceVector() const;
    void scatterFreeDofs(const VecX& delta);
    VecX gatherFreeDofs() const;
};
```

**装配逻辑** (assemble):
1. `std::vector<Triplet>` reserve `elements.size() * elemDofs * elemDofs`
2. 对每个 element:
   - 收集两节点 DOF → `dispVec`
   - `if constexpr (ElemType::hasTransformation)`:
     - `T = ElemType::computeTransformation(xA, xB, elem.refVector)`
     - `dispLocal = T * dispGlobal`
     - `xA_local = (0,0,0)`, `xB_local = (L,0,0)`
     - `result = ElemType::computeElement(xA_local, xB_local, dispLocal, props)`
     - `re_g = T^T * re_l`, `ke_g = T^T * ke_l * T`
   - 否则直接调用
   - 构建 location map → Triplet 收集
3. `K.setFromTriplets(...)`

**注意**: `computeTransformation` 由 ElemType 提供，BeamModel 只调用统一接口。2D 元素内部忽略 refVector；3D 元素使用 refVector 构建局部坐标系。BeamModel 不嵌入 2D/3D 坐标转换策略。

### 2A.5 LinearSolver

**输出**: `include/BeamLib/Solver/LinearSolver.h`

```cpp
struct SolveResult {
    VecX x;
    bool success = false;
};

struct LinearSolver {
    static SolveResult solve(const SpMat& K, const VecX& rhs);
};
```

`SparseLU` 需报告奇异分解。

### 2A.6 NewtonRaphson 求解器

**输出**: `include/BeamLib/Solver/NewtonRaphson.h`

```cpp
struct NRConfig {
    double tol = 1e-7;
    int maxIter = 20;
};

struct NRResult {
    bool converged = false;
    int iterations = 0;
    double finalResidual = 0.0;
};

template <typename ElemType>
struct NewtonRaphsonSolver {
    static NRResult solveOneStep(
        BeamModel<ElemType>& model,
        const VecX& Fext,
        const NRConfig& config = {});
};
```

**NR 循环**:
```
repeat:
    assemble(R_int, K)
    R = R_int - F_ext
    if ||R|| < tol: converged
    delta = LinearSolver::solve(K, -R)
    model.scatterFreeDofs(delta)
```

### 2A.7 测试

**测试 1**: `tests/test_eb2d_cantilever.cpp`
```
E=200e9, A=0.01, Iz=8.333e-6
L=1.0, 10 单元沿 x 轴
节点 0 全固定, 节点 10 施加 F_z = -1000
```
- u_z = -FL³/(3EI) = -2.0e-4, θ_y = -FL²/(2EI) = -3.0e-4
- 多单元位移: 相对误差 < **1e-9**
- **1 次 NR 迭代**

**测试 2**: Patch tests (合并到一个文件 `tests/test_eb2d_patch.cpp`)
- 纯轴向: 两端不同轴向位移 → 均匀应变 → 精确
- 刚体平移: 所有节点同位移 → 零内力
- 单元素刚度检查: K(0,0)=EA/L, K(1,1)=12EI/L³ → 相对误差 < 1e-12
- **单元素 θ_y 符号测试**: F_z < 0 时验证 θ_y 的预期符号

**测试 3**: 纯端弯矩 `tests/test_eb2d_moment.cpp`
```
节点 10 施加 M_y (纯弯矩)
```
- 线性弯矩分布，位移为二次函数
- 验证弯曲刚度

### 验收标准
- [ ] 理论文档 PDF 完成，含三个附表
- [ ] 编译零 warning
- [ ] 悬臂梁 u_z, θ_y 相对误差 < 1e-9
- [ ] 1 次 NR 迭代收敛
- [ ] Patch tests 全部通过 (含 θ_y 符号测试)
- [ ] `ctest` 全部 PASS

### Review 重点
- 理论文档中刚度矩阵推导与代码实现是否一致
- 残差符号: `re = ke * dispVec` (不是 `-ke * dispVec`)
- 2D 变换矩阵: sin/cos 取 dx_z 而非 dx_y（xz 平面约定）
- θ_y 正方向定义是否一致

---

## Batch 2B: EB2D 后处理 + 扩展测试

### 目标
在 Batch 2A 的基础上补全: 质量矩阵装配、支反力提取、单元内力、简支梁测试、旋转梁测试。

### 2B.1 质量矩阵装配

实现 `BeamModel::assembleMass(SpMat& M)` — 与刚度装配类似，调用 `ElemType::computeMass`，经坐标变换后装入稀疏矩阵。

### 2B.2 支反力提取

**策略**: 全系统重装配（含固定 DOF 位置），取固定 DOF 处的 `R_int - F_ext` 残差分量。注意约束点可能有外力，不能只取 `R_int`。

实现为 `BeamModel` 的方法或独立函数:
```cpp
std::vector<ReactionForce> computeReactions(const BeamModel<ElemType>& model);
```

### 2B.3 单元内力提取

**输出**: `include/BeamLib/PostProcess/InternalForces.h`

EB2D 内力从局部位移通过 Hermite 曲率 B 矩阵在单元两端计算。完整的端点力/弯矩公式在 EB2D 理论文档附表中给出，代码必须按理论文档的公式→索引映射表实现。

元素 static method 提供原始两端内力值。PostProcess 层此阶段只做简单封装。

### 2B.4 Nonzero Prescribed Displacement (设计 + 实现 + 测试)

**硬性要求**: Batch 2B **必须**完成 elimination + RHS correction 的设计、实现和测试。这是 Batch 4 的前置门槛。

实现内容:
- 当 `node.fixed[i] = true` 且 prescribed value ≠ 0 时，修正右端项
- 实现在 `BeamModel::assemble` 或 `DofMap` 层
- Node 结构增加 prescribed value 存储（或通过单独的 map）

测试: 两端固支梁，一端施加非零沉降位移 → 验证位移场和支反力与解析解一致。

### 2B.5 测试

**测试 1**: `tests/test_eb2d_postprocess.cpp` — 悬臂梁后处理
```
同 Batch 2A 悬臂梁设置
```
- 支反力: R_z(node0) = 1000, M_y(node0) = 1000 (相对误差 < 1e-9)
- 单元内力: 根部 V_z = 1000, M_y = 1000; 末端 V_z = 1000, M_y = 0

**测试 2**: `tests/test_eb2d_simply_supported.cpp`
```
E=200e9, A=0.01, Iz=8.333e-6
L=2.0, 20 单元
节点 0: u_x=0, u_z=0 (fixed), θ_y free
节点 20: u_z=0 (fixed), θ_y free, u_x free
均布载荷 q=1000 N/m
```
- 等效节点力: **测试代码内手动组装** (qL_e/2 集中力 + qL_e²/12 弯矩)，不依赖 LoadManager
- 跨中挠度: w_max = 5qL⁴/(384EI), 相对误差 < 1e-8
- 支反力: R_A = R_B = qL/2

**测试 3**: `tests/test_eb2d_rotated_beam.cpp` — 旋转梁 (两个子测试)

**子测试 3a**: 局部横向力
```
45° 梁: xA=(0,0,0), xB=(L/√2, 0, L/√2)
施加局部横向力（沿局部 z 方向），变换到全局坐标系后施加
将全局位移变换回局部坐标，与水平梁（同参数、同局部力）对比
```
- 局部位移相对误差 < 1e-9

**子测试 3b**: 全局力
```
同一 45° 梁，施加全局竖直力 F_z
```
- 验证全局位移中 u_x 和 u_z 均非零（局部轴向/横向耦合）
- 与手算投影值对比

### 验收标准
- [ ] 质量矩阵装配编译通过
- [ ] 支反力正确（悬臂梁、简支梁）
- [ ] 单元内力正确（悬臂梁剪力/弯矩分布）
- [ ] 简支梁跨中挠度精度 < 1e-8
- [ ] 旋转梁两个子测试均通过
- [ ] **Nonzero prescribed displacement 实现并通过测试** (Batch 2B 必须完成项)
- [ ] `ctest` 全部 PASS

### Review 重点
- 支反力提取逻辑: 全系统重装配 vs 装配时记录
- 简支梁等效节点力: 手动组装的公式是否正确
- 旋转梁测试的分解是否严谨
- Elimination + RHS correction 实现的正确性

---

## Batch 3: EB3D

### 目标
扩展到 3D Euler-Bernoulli，引入 3D 坐标变换和扭转。

### 前置: EB3D 理论文档

**输出**: `docs/theory/02_euler_bernoulli_3d.tex` → PDF

新增内容 (相对 EB2D):
- 双平面弯曲 (xy + xz)
- 扭转自由度 (GIx/L)
- 3D 坐标变换矩阵推导 (含 refVector 构建局部 y 轴)

**必须包含附表**: 符号约定表、公式→C++ 索引映射表 (12×12 完整)、验证容差表。

### 3.1 BeamMath3D 迁移

**输出**: `include/BeamLib/Math/BeamMath3D.h`

从 `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\include\BeamMath.h` 复制，仅改 namespace。此文件同时服务 EB3D 变换和后续 GE3D。

### 3.2 EulerBernoulli3D 单元

**输出**: `include/BeamLib/Element/EulerBernoulli3D.h`

```cpp
struct EulerBernoulli3D {
    static constexpr int nDofsPerNode = 6;
    static constexpr bool hasTransformation = true;
    static constexpr bool isLinear = true;
    // ...
};
```

**刚度矩阵** (局部, 12×12):
- 轴向 (DOF 0,6): EA/L
- 扭转 (DOF 3,9): GIx/L
- xy 平面弯曲 (DOF 1,5,7,11): 用 EIz
- xz 平面弯曲 (DOF 2,4,8,10): 用 EIy

**变换矩阵**: 使用 `ElementConn::refVector` 构建局部坐标系。算法提取为共用函数（后续 Timo3D、GE3D 复用）。默认 refVector=(0,1,0)，梁平行于此向量时 fallback (0,0,1)。

**质量矩阵**: 一致质量，从理论文档推导。

### 3.3 测试

**测试 1**: `tests/test_eb3d_cantilever.cpp`
```
E=200e9, G=80e9, rho=7800
A=0.01, Iy=Iz=8.333e-6, Ix=16.667e-6
L=1.0, 10 单元沿 x 轴
节点 0 全固定, 节点 10 施加 F_z = -1000
```
- u_z = -2.0e-4, θ_y = -3.0e-4, 其余 DOF ≈ 0
- 多单元位移: 相对误差 < 1e-9
- 近零 DOF: 绝对值 < 1e-14
- 1 次 NR 迭代

**测试 2**: 扭转测试
```
节点 10 施加 T_x (扭矩)
```
- θ_x(tip) = TL/(GIx), 相对误差 < 1e-9

**测试 3**: 双向弯曲
```
节点 10 同时施加 F_y 和 F_z
```
- 两方向位移独立（因 Iy=Iz 可验证对称性）

**测试 4**: 单元素刚度检查
- K(0,0)=EA/L, K(3,3)=GIx/L, K(1,1)=12EIz/L³, K(2,2)=12EIy/L³
- 容差: 相对 1e-12

**测试 5**: 纯端弯矩
- 同 EB2D 的端弯矩测试，3D 版本

### 验收标准
- [ ] 理论文档完成，含三个附表
- [ ] 12×12 刚度矩阵对称
- [ ] 单元素矩阵检查通过 (1e-12)
- [ ] 1 次 NR 迭代，u_z 精度 < 1e-9
- [ ] 非加载方向 DOF 绝对值 < 1e-14
- [ ] `ctest` PASS

### Review 重点
- DOF 索引: MATLAB 1-based → C++ 0-based 映射 (用理论文档附表验证)
- xz 平面弯曲的耦合 DOF (2,4,8,10) 是否正确
- refVector 处理逻辑与 2D Rotation2D 的隔离

---

## Batch 4: Timoshenko 2D

### 目标
实现含剪切变形的 2D 梁单元，验证与 EB 在细长梁极限下的一致性。

### 前置: Timoshenko 理论文档

**输出**: `docs/theory/03_timoshenko_2d.tex` → PDF

关键内容:
1. Timoshenko 假设: 截面仍为平面，但不再垂直于中性轴
2. 剪切应变 γ = dw/dx - θ
3. 弱形式包含弯曲和剪切两部分
4. **Shear locking 问题**: 全积分下薄梁锁定
5. **Locking-free 方案: Φ-corrected 闭合形式刚度矩阵**
   - `Φ = 12EI/(κGAL²)`
   - 弯曲刚度项含 `1/(1+Φ)` 修正
   - Φ→0 自然退化为 EB
   - 不使用 EAS，不使用选择性缩减积分
6. 质量矩阵: 含旋转惯量 ρI，不含剪切惯量修正
7. 与 EB 的关系: Φ → 0 时精确退化

**必须包含附表**: 符号约定表、完整 6×6 刚度矩阵（含 Φ 项）、验证容差表。

### 4.1 Timoshenko2D 单元

**输出**: `include/BeamLib/Element/Timoshenko2D.h`

```cpp
struct Timoshenko2D {
    static constexpr int nDofsPerNode = 3;
    static constexpr bool hasTransformation = true;
    static constexpr bool isLinear = true;
    // ...
};
```

DOF 顺序同 EB2D: `[u_x, u_z, θ_y]`

**刚度矩阵**: 2-node Φ-corrected 闭合形式。轴向部分与 EB 相同，弯曲部分的分母含 `(1+Φ)`:
```
K_bend(1,1) = 12EI / (L³(1+Φ))
K_bend(1,2) = 6EI / (L²(1+Φ))
K_bend(2,2) = (4+Φ)EI / (L(1+Φ))
K_bend(2,5) = (2-Φ)EI / (L(1+Φ))
...
```

**变换矩阵**: 复用 `Rotation2D::compute`。

**质量矩阵**: 一致 + 含旋转惯量 ρI 项。

### 4.2 Nonzero Prescribed Displacement

Batch 2B 已完成 elimination + RHS correction 的设计、实现和测试（硬性前置）。本 Batch 直接复用，不再提供替代路径。

### 4.3 测试

**测试 1**: 悬臂梁
```
E=200e9, G=80e9, kappa=5/6, A=0.01, Iz=8.333e-6
L=1.0, 10 单元, F_z = -1000
```
- 解析解: `u_z = -FL³/(3EI) - FL/(κGA)` (含剪切项)
- 相对误差 < 1e-9
- 1 次 NR 迭代

**测试 2**: 与 EB 对比 — 细长梁
```
L/h = 100 (很细长)
```
- Timoshenko 结果 ≈ EB 结果 (差异 < 0.1%)

**测试 3**: 与 EB 对比 — 深梁
```
L/h = 5 (粗短)
```
- Timoshenko 挠度明显大于 EB（剪切贡献显著）
- 与含剪切解析解一致

**测试 4**: Shear locking 检查
```
L/h = 1000 (极细长), 2~5 单元
```
- 结果不应出现 locking (挠度不应远小于 EB 解)
- 与 EB 解差异 < 0.1%

**测试 5**: 单元素矩阵检查
- 验证 Φ→0 时退化为 EB 矩阵
- 容差: 1e-12

### 验收标准
- [ ] 理论文档完成，shear locking 方案为 Φ-corrected 闭合形式
- [ ] 悬臂梁精度 < 1e-9
- [ ] 细长梁极限与 EB 一致
- [ ] 深梁剪切效应正确
- [ ] 无 shear locking (极细长粗网格)
- [ ] Φ→0 退化为 EB (单元素检查)
- [ ] `ctest` PASS

### Review 重点
- Φ-corrected 矩阵的每一项是否正确（对照理论文档附表）
- 质量矩阵中旋转惯量项
- kappa 取值: 使用 `props.kappa()` (即 kappa_z)

---

## Batch 5: Timoshenko 3D

### 目标
3D Timoshenko 梁，含双平面剪切和扭转。

### 前置: 理论文档

**输出**: `docs/theory/04_timoshenko_3d.tex` → PDF

相对 2D 新增: 双平面剪切 (kappa_y/kappa_z 分别修正)、扭转。

### 5.1 Timoshenko3D 单元

**输出**: `include/BeamLib/Element/Timoshenko3D.h`

```cpp
struct Timoshenko3D {
    static constexpr int nDofsPerNode = 6;
    static constexpr bool hasTransformation = true;
    static constexpr bool isLinear = true;
    // ...
};
```

DOF 顺序同 EB3D。刚度矩阵 12×12:
- xy 平面: `Φ_z = 12EIz/(kappa_y·GA·L²)`
- xz 平面: `Φ_y = 12EIy/(kappa_z·GA·L²)`

**变换矩阵**: 复用 EB3D 的 3D 变换（含 refVector）。

**质量矩阵**: 一致 + 含旋转惯量 ρI。

### 5.2 测试

- 悬臂梁 (同 EB3D 设置，加入剪切) → 精度 < 1e-9
- 与 EB3D 对比 (细长/深梁)
- 扭转测试 (应与 EB3D 相同)
- Shear locking 检查
- 单元素矩阵: Φ→0 退化为 EB3D

### 验收标准
- [ ] 理论文档完成
- [ ] 12×12 刚度矩阵对称
- [ ] 悬臂梁精度 < 1e-9
- [ ] 与 EB3D 在细长极限一致
- [ ] kappa_y/kappa_z 分别生效
- [ ] 无 shear locking
- [ ] `ctest` PASS

---

## Batch 6: GeomExact3D + 非线性求解

### 目标
迁移已有 3D 几何精确梁实现，建立非线性求解能力。

### 前置: GE 理论文档

**输出**: `docs/theory/05_geometrically_exact.tex` → PDF

**关键内容**:
1. Reissner-Simo 梁理论 (可剪切几何精确梁)
2. 有限转动参数化: Rodrigues 公式 `R(θ) = I + sin|θ|/|θ| · Θ̂ + (1-cos|θ|)/|θ|² · Θ̂²`
3. 应变度量: Γ (axial/shear strain), K (curvature/torsion)
4. 弱形式线性化 → 切线刚度
5. B-矩阵构造
6. 几何刚度贡献
7. Newton-Raphson 载荷步策略
8. 1-point Gauss 积分: Phase 1 保持此选择以与已验证 C++ 代码一致；多点积分评估留后续

**必须包含附表**: 符号约定表、索引映射表、验证容差表。

### 6.1 BeamMath3D

Batch 3 已迁移。

### 6.2 GeomExact3D 单元

**输出**: `include/BeamLib/Element/GeomExact3D.h` + `src/GeomExact3D.cpp`

从 `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\src\BeamElement3D.cpp` 迁移:
- `computeElement`: 1-point Gauss, 12×1 residual + 12×12 tangent
- `computeResidual`, `computeStiffBlock`: 内部 static helpers
- 参数签名: `(const Mat3& C, const Mat3& D)` → `(const SectionProperties& props)`
- **算法零改动** — 仅 namespace 和参数适配

`computeTransformation` 从 `BeamModel::calcTransformationMatrix` 迁移，使用 `refVector` 参数。

**局部坐标传递注意**: BeamModel 装配时传入 `xA_local=(0,0,0)`, `xB_local=(L,0,0)`。参考 C++ 代码使用 `lambda*xA`, `lambda*xB`（保留全局位置信息经旋转）。对 EB/Timo 线性单元两者等价（只依赖长度 L）。对 GE3D 需验证等价性 — 见 Batch 6 测试 4。

**CMakeLists.txt 更新**: INTERFACE → static library，添加 `src/GeomExact3D.cpp`。

### 6.3 Newton-Raphson 载荷步

扩展 NR solver，添加 multi-step:

```cpp
struct LoadStepConfig {
    int nSteps = 1;
    NRConfig nrConfig = {};
};

template <typename ElemType>
struct NewtonRaphsonSolver {
    static NRResult solveOneStep(...);
    static std::vector<NRResult> solveLoadStepping(
        BeamModel<ElemType>& model,
        const VecX& Ftotal,
        const LoadStepConfig& config);
};
```

载荷步: `F_n = (n/nSteps) * F_total`，每步从上一步结果继续。

### 6.4 测试

**测试 1**: 小载荷悬臂梁 `tests/test_ge3d_cantilever.cpp`
```
E=1e6, G=5e5, kappa_y=kappa_z=1.0, A=1, Iy=Iz=Ix=1
L=1.0, 10 单元, F_z = 1.0
1 步 NR, tol=1e-7
```
- 完整 Timoshenko 解析解: **`u_z = FL³/(3EI) + FL/(κGA) = 1/(3×1e6) + 1/(1×5e5×1) = 3.333e-7 + 2.0e-6 = 2.333e-6`**
- 相对误差 < 1e-4

**测试 2**: 大变形 `tests/test_ge3d_large_deformation.cpp`
```
同上参数, F_z = 5e5, 20 步载荷步
```
- 20 步全收敛
- 无 NaN/Inf
- u_z > 0, u_x < 0 (大变形缩短效应)

**测试 3**: 文献定量基准 — **必须通过**

**前置任务**: Batch 6 开始前须从文献中提取完整基准数据（几何、材料、载荷、离散方案、目标位移值、数据来源表格），写入理论文档或独立基准文件 `docs/theory/benchmark_bathe_bolourchi.md`，使测试可直接编写为闭合验收用例。

**Bathe & Bolourchi 45° cantilever bend** (1979):
- 45° 弯曲悬臂梁，末端集中力
- 对比文献中的末端位移数据
- 相对误差 < 5% (考虑不同离散化)

**测试 4**: GE3D 局部坐标 invariance 测试

验证 `xA_local=(0,0,0)`/`xB_local=(L,0,0)` 与参考代码 `lambda*xA`/`lambda*xB` 的等价性:
```
同测试 1 的悬臂梁，但节点初始位置从 x=5.0 开始 (非原点)，
或沿非轴对齐方向放置，验证结果与测试 1 完全一致
```
- 末端位移差异 < 1e-12 (应精确相同)
- 若不一致，说明 GE3D 的 `computeElement` 依赖绝对位置，需调整装配层

### 验收标准
- [ ] 理论文档完成
- [ ] 文献基准数据已提取并文档化
- [ ] `computeElement` 与原代码 diff 仅有 namespace/参数变更
- [ ] 小载荷: ≤ 5 次 NR 迭代，u_z 误差 < 1e-4
- [ ] 大变形: 20 步收敛，无 NaN，u_x < 0
- [ ] **文献基准定量通过** (Bathe & Bolourchi)
- [ ] 局部坐标 invariance 通过
- [ ] `ctest` PASS

### Review 重点
- `computeElement` 算法与原文件逐行 diff
- `assemble` 中坐标变换: `dispLocal = T * dispGlobal`, `re_g = T^T * re_l`, `ke_g = T^T * ke_l * T`
- NR 残差: `R_int - F_ext`，符号与参考一致
- 小载荷基准使用完整 Timoshenko 解析解 (含剪切项)

---

## Batch 7: GeomExact2D

### 目标
独立推导并实现 2D 几何精确梁。

### 实现路径: **独立 2D 推导**

使用标量 θ_y 和 2D 旋转矩阵 (cos/sin)，不从 GE3D 退化。理由:
- 教学清晰: 2D 大转动用 cos/sin 比 Rodrigues 更直观
- 调试简单: 独立实现更容易隔离 bug
- 用 GE3D 面内结果做交叉验证

### 7.1 GeomExact2D 单元

**输出**: `include/BeamLib/Element/GeomExact2D.h`

```cpp
struct GeomExact2D {
    static constexpr int nDofsPerNode = 3;
    static constexpr bool hasTransformation = true;
    static constexpr bool isLinear = false;
    // ...
};
```

2D GE 梁在 xz 平面内。运动学: 位置 (x + u_x, z + u_z)，转角 θ_y (标量)。

### 7.2 测试

- 小载荷: 与 Timo2D 对比 (相对误差 < 1e-4)
- 大变形: 大角度弯曲，NR 收敛
- 与 GE3D 面内结果交叉验证 (相对误差 < 1e-6)

### 验收标准
- [ ] 理论文档 (可在 05_geometrically_exact.tex 中增加 2D 章节)
- [ ] 小载荷与 Timo2D 一致
- [ ] 大变形收敛、无 NaN
- [ ] 与 GE3D 面内结果交叉验证通过
- [ ] `ctest` PASS

---

## Batch 8: 隐式时间积分

### 目标
实现 Newmark-β 和 Generalized-α 隐式时间积分器。

### 前置: 理论文档

**输出**: `docs/theory/06_time_integration.tex` → PDF

**内容**:
1. 半离散运动方程: Mü + Cu̇ + R_int(u) = F_ext(t)
2. Newmark-β 族: β, γ 参数
   - 稳定性条件: `γ ≥ 0.5` 且 `β ≥ 0.25(γ + 0.5)²`
   - 默认: β=0.25, γ=0.5 (常加速度法，无条件稳定，二阶精度，无数值阻尼)
3. Generalized-α (Chung-Hulbert):
   - `α_m = (2ρ_∞ - 1)/(ρ_∞ + 1)`
   - `α_f = ρ_∞/(ρ_∞ + 1)`
   - `γ = 0.5 + α_f - α_m`
   - `β = 0.25(1 + α_f - α_m)²`
   - **必须明确定义评估时刻**: `R_int`, velocity, acceleration 在 α_m 和 α_f 加权时间点的评估位置 (如 `R_int(u_{n+1-α_f})`, `Ma_{n+1-α_m}`)，否则实现可能偏离推导约定
4. 非线性问题: 每个时间步内 NR 迭代
5. 稳定性与精度分析
6. Rayleigh 阻尼: C = α_M · M + β_K · K

### 8.1 质量矩阵确认

所有 6 种单元的 `computeMass` 实现状态:
- EB2D/3D: Batch 2A/3 已实现一致质量
- Timo2D/3D: Batch 4/5 已实现一致 + 旋转惯量
- GE3D: **待定** — 参考 C++ 代码用 lumped mass (为显式设计)；隐式/模态需 consistent mass，须在 GE 理论文档 (Batch 6) 中推导确认
- GE2D: 跟随 GE3D 决定

**确认**: 所有质量矩阵在进入 Batch 8 前已就绪。GE3D/GE2D 的质量矩阵类型须在 Batch 6/7 理论文档中确定，不能拖到 Batch 8。如有缺失，在此 batch 补全。

### 8.2 Newmark-β

**输出**: `include/BeamLib/Solver/ImplicitTimeIntegrator.h`

```cpp
struct NewmarkConfig {
    double beta = 0.25;
    double gamma = 0.5;
    double dt = 0.0;
    int numSteps = 0;
    double alpha_damping = 0.0;
    double beta_damping = 0.0;
};

template <typename ElemType>
class NewmarkIntegrator {
public:
    NewmarkIntegrator(BeamModel<ElemType>& model, const NewmarkConfig& config);
    NRResult stepForward();
    int currentStep() const;
    double currentTime() const;
    double kineticEnergy() const;
    double strainEnergy() const;
};
```

### 8.3 Generalized-α

```cpp
struct GenAlphaConfig {
    double rho_inf = 0.9;
    double dt = 0.0;
    int numSteps = 0;
    double alpha_damping = 0.0;
    double beta_damping = 0.0;
};
```

从 ρ_∞ 推导 α_m, α_f, β, γ，使用 Chung-Hulbert 公式。

### 8.4 测试

**测试 1**: EB2D 悬臂梁自由振动 `tests/test_implicit_dynamics.cpp`
```
初始条件: 节点 10 施加初始 u_z 位移
Newmark β=0.25, γ=0.5, 无阻尼
运行足够步数覆盖多个振动周期
```
- 总能量 (KE + SE) 波动 < 1%

**测试 2**: Generalized-α 数值耗散验证
```
ρ_∞ = 0.9 → 高频有耗散，低频保持
```
- 与 Newmark 对比: 低频响应一致，高频逐渐衰减

**测试 3**: 非线性动力学 (GE3D)
- 大振幅振动，NR 每步收敛

### 验收标准
- [ ] 理论文档完成 (含稳定性条件和 Chung-Hulbert 公式)
- [ ] Newmark: 无阻尼能量守恒 < 1%
- [ ] Gen-α: 数值耗散可控
- [ ] 非线性动力学: 每步 NR 收敛
- [ ] `ctest` PASS

---

## Batch 9: 模态分析

### 目标
实现广义特征值求解，提取固有频率和振型。

### 9.1 模态求解器

**输出**: `include/BeamLib/Solver/ModalSolver.h`

```cpp
struct ModalConfig {
    int numModes = 6;
};

struct ModalResult {
    VecX frequencies;
    VecX angularFrequencies;
    std::vector<VecX> modeShapes;
};

template <typename ElemType>
struct ModalSolver {
    static ModalResult solve(
        BeamModel<ElemType>& model,
        const ModalConfig& config = {});
};
```

**实现要求**:
- 装配 free-DOF dense `Kff` 和 `Mff` (从稀疏矩阵提取)
- 使用 `Eigen::GeneralizedSelfAdjointEigenSolver<MatX>`
- 验证 `Mff` 为 SPD（不是则说明质量矩阵或约束有问题）
- 约束不足 → K 奇异 → 需报告错误
- **Phase 1 限教学/中小规模模型**，大规模稀疏模态分析不在范围内

### 9.2 测试

**测试 1**: EB2D 悬臂梁解析频率 `tests/test_modal_cantilever.cpp`
```
均匀梁, 解析: ω_n = (β_n L)² √(EI/ρAL⁴)
β_1 L = 1.8751, β_2 L = 4.6941, β_3 L = 7.8548
```
- 前 3 阶频率与解析解相对误差 < 1e-6
- 网格加密后误差应减小

**测试 2**: 模态正交性
- 验证 `Φᵀ M Φ ≈ I` (归一化后)
- 验证 `Φᵀ K Φ ≈ diag(ω²)`

**测试 3**: EB3D 模态
- 两个弯曲平面的频率 (Iy=Iz 时应重合)
- 扭转频率

### 验收标准
- [ ] 前 3 阶频率误差 < 1e-6
- [ ] 模态正交性通过
- [ ] Mff SPD 验证
- [ ] `ctest` PASS

---

## Batch 10: 最终集成 + API 封装

### 目标
LoadManager、Type-erased API wrapper、便捷层，最终集成测试。

### 10.1 LoadManager

**输出**: `include/BeamLib/Load/LoadTypes.h` + `include/BeamLib/Load/LoadManager.h`

```cpp
template <typename ElemType>
class LoadManager {
public:
    void addPointLoad(int nodeId, int dofIndex, double value);
    void addPrescribedDisplacement(int nodeId, int dofIndex, double value);
    void addUniformDistributedLoad(int elemId, int direction, double q,
                                    CoordinateSystem cs = CoordinateSystem::Local);
    void applyToModel(BeamModel<ElemType>& model, double loadFactor = 1.0);
    void clear();
};
```

分布载荷: 自动转换为等效节点力（一致等效，含弯矩），`direction` 指定力的方向，`cs` 指定局部/全局坐标系。

### 10.2 Type-erased API (预留)

**输出**: `include/BeamLib/API/AnalysisModel.h`

```cpp
class IAnalysisModel {
public:
    virtual ~IAnalysisModel() = default;
    virtual void solve() = 0;
    virtual AnalysisResult getResult() const = 0;
};

template <typename ElemType>
class TypedAnalysisModel : public IAnalysisModel { ... };
```

Phase 2/3 接口预留。Phase 1 提供接口定义和简单实现。

### 10.3 便捷层 (Convenience API)

```cpp
namespace beamlib::easy {
    auto cantilever_eb2d(double L, int nElems, double E, double A, double Iz);
    auto simply_supported_eb2d(double L, int nElems, ...);
    // ...
}
```

### 10.4 最终集成测试

- 所有已有测试仍 PASS
- 交叉验证矩阵:
  - GE3D 小载荷 ≈ Timo3D (含剪切解析解)
  - Timo (L/h→∞) ≈ EB
  - 2D 结果 ≈ 3D 面内结果 (含 Iy/Iz 映射)
- 综合测试: 多构件框架结构 (如 portal frame)

### 验收标准
- [ ] LoadManager 功能完整 (含分布载荷)
- [ ] 全部已有测试 PASS
- [ ] 交叉验证通过
- [ ] API 接口定义清晰
- [ ] 便捷层可用

---

## 执行顺序和依赖

```
Batch 1 (脚手架)          → 无依赖
Batch 2A (EB2D 核心)      → 依赖 Batch 1
Batch 2B (EB2D 后处理)    → 依赖 Batch 2A
Batch 3 (EB3D)            → 依赖 Batch 2A
Batch 4 (Timo2D)          → 依赖 Batch 2A (+ Batch 2B 的经验)
Batch 5 (Timo3D)          → 依赖 Batch 2A
Batch 6 (GE3D)            → 依赖 Batch 2A
Batch 7 (GE2D)            → 依赖 Batch 6
Batch 8 (隐式动力学)       → 依赖 Batch 2B~7 的质量矩阵
Batch 9 (模态分析)         → 依赖 Batch 2B~7 的质量矩阵
Batch 10 (集成)            → 依赖全部
```

```
                     Batch 1
                        │
                     Batch 2A (EB2D 核心 pipeline)
                   ╱    │    ╲
             Batch 2B  Batch 3  Batch 4   Batch 6
             (后处理)  (EB3D)  (Timo2D)   (GE3D)
                │       │       │           │
                │    Batch 5  Batch 5    Batch 7
                │    (Timo3D)            (GE2D)
                 ╲      │      ╱
                  Batch 8 + 9 (动力学 + 模态)
                        │
                     Batch 10 (集成)
```

建议执行顺序: 1 → 2A → 2B → 3 → 4 → 5 → 6 → 7 → 8 → 9 → 10

Batch 3/4/6 之间无相互依赖，但从简到繁的顺序有助于积累经验。

## 预估工作量

| Batch | 主要内容 | 预估 |
|-------|---------|------|
| 1 | 脚手架 | 0.5 天 |
| 2A | EB2D 核心 pipeline | 2-3 天 (含理论文档) |
| 2B | EB2D 后处理 + 扩展测试 | 1-2 天 |
| 3 | EB3D | 1-2 天 |
| 4 | Timo2D | 2-3 天 (含理论 + locking 研究) |
| 5 | Timo3D | 1-2 天 |
| 6 | GE3D + NR | 3-4 天 (含理论文档 + 迁移验证 + 文献基准) |
| 7 | GE2D | 2-3 天 (独立推导) |
| 8 | 隐式动力学 | 3-4 天 |
| 9 | 模态分析 | 1-2 天 |
| 10 | 集成 + API | 2-3 天 |
| **合计** | | **~20-30 天 (4-6 周)** |

## Git 分支策略

```
main                             ← review 通过后合并
├── feat/batch-1-scaffold
├── feat/batch-2a-eb2d-core
├── feat/batch-2b-eb2d-postprocess
├── feat/batch-3-eb3d
├── feat/batch-4-timo2d
├── feat/batch-5-timo3d
├── feat/batch-6-ge3d
├── feat/batch-7-ge2d
├── feat/batch-8-implicit-dynamics
├── feat/batch-9-modal
└── feat/batch-10-integration
```

每个 Batch 完成后提交到对应 branch，Codex review 后合并到 main。
