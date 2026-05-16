# BeamLib — 项目总规范 (v3)

## 1. 项目概述

C++17 + Eigen3 独立梁单元有限元库。

**双重定位**:
- **(A) 教学工具** (主要): 面向结构力学/有限元课程，学生可通过 GUI 或 AI 自然语言接口建模与分析梁/框架结构
- **(B) 工程分析** (次要): 可用于框架建筑骨架、机械杆件等简单工程结构力学分析

## 2. 产品路线图

```
Phase 1: 核心求解引擎 (C++ library)     ← 当前重点
    6 种梁单元 + 3 种分析类型 + 后处理数据提取
    配套文档: 理论推导手册 + 应用案例教程
    
Phase 2: GUI 可视化 (Teaching frontend)
    交互式前/后处理，面向教学场景的可视化
    
Phase 3: AI 自然语言接口
    Python scripted, LLM-driven 建模与分析
```

Phase 1 的 API 设计必须考虑 Phase 2 (GUI runtime dispatch) 和 Phase 3 (Python bindings) 的需要。

## 3. Phase 1 完整范围

### 3.1 单元类型 (6 种)

| 单元 | 维度 | 理论 | DOFs/node | 适用范围 | 关键特征 |
|------|------|------|-----------|----------|----------|
| EulerBernoulli2D | 2D | Euler-Bernoulli | 3 | 细长构件 (L/h > 20) | 无剪切效应，静力/动力均排除剪切 |
| EulerBernoulli3D | 3D | Euler-Bernoulli | 6 | 细长构件 (L/h > 20) | 无剪切效应，含扭转 |
| Timoshenko2D | 2D | Timoshenko | 3 | 一般构件 (任意 L/h) | 含剪切变形 + 旋转惯量，Φ-corrected locking-free |
| Timoshenko3D | 3D | Timoshenko | 6 | 一般构件 (任意 L/h) | 含剪切变形 + 旋转惯量，含扭转 |
| GeomExact2D | 2D | Geometrically Exact | 3 | 大变形/大转动 | Reissner-Simo 可剪切梁，2D 独立推导 |
| GeomExact3D | 3D | Geometrically Exact | 6 | 大变形/大转动 | Reissner-Simo 可剪切梁，Rodrigues 参数化 |

**理论边界说明**:
- EB 完全排除剪切效应（包括静力刚度和动力学惯量）
- Timoshenko 包含剪切变形和旋转惯量 ρI；Phase 1 不含剪切惯量修正
- GE 是剪切可变形的 Reissner-Simo 梁，包含 axial/shear strain 和 curvature/torsion

### 3.2 分析类型 (3 种)

| 分析类型 | 求解器 | 适用场景 |
|----------|--------|----------|
| 静力分析 | 线性 SparseLU / Newton-Raphson | 线性和非线性静力问题 |
| 动力分析 (隐式) | Newmark-β / Generalized-α | 瞬态响应 |
| 模态分析 | 广义特征值 Kφ = ω²Mφ | 固有频率和振型 |

**注意**: Phase 1 不包含显式动力学。所有时间积分均为隐式格式。分析类型接口须可扩展，以后可增加 buckling、频响分析等。

### 3.3 后处理

Phase 1 输出结构化数据（非可视化），包括：
- **节点位移** (平移 + 转角)
- **支反力** (固定 DOF 的反力)
- **单元内力**: 轴力 N, 剪力 V (2D: V_z; 3D: V_y, V_z), 弯矩 M (2D: M_y; 3D: M_y, M_z), 扭矩 T (3D only)

**反力恢复策略**: 全系统重装配（含固定 DOF），取固定 DOF 处的 `R_int - F_ext` 分量。实现简单，不侵入 assemble 逻辑，对教学规模模型性能可接受。

**内力后处理分层**: 元素 static method 提供原始数据（两端内力值）；PostProcess 层负责 sampling、smoothing、图表约定等后处理策略。Phase 1 的 PostProcess 层只做基本的每单元两端内力。

### 3.4 质量矩阵策略

| 单元 | 质量矩阵 | 旋转惯量 | 备注 |
|------|---------|---------|------|
| EB2D | 一致 (consistent) | 不含 | 标准教科书形式，无旋转惯量 |
| EB3D | 一致 (consistent) | 扭转含 ρIx | 弯曲无旋转惯量；扭转 DOF (θ_x) 含极惯性矩项 ρIx，否则扭转模态缺失 |
| Timo2D/3D | 一致 + 旋转惯量 ρI | 含 | 对短粗梁高阶模态有影响 |
| GE3D | 待定 (理论文档确认) | 含 | 参考代码用 lumped (显式)；隐式/模态需 consistent，须单独推导验证 |
| GE2D | 待定 (理论文档确认) | 含 | 跟随 GE3D 决定 |

此决定影响 Batch 8 (动力学) 和 Batch 9 (模态)，必须在各单元理论文档中确定。

## 4. 硬约束

- **namespace**: 所有代码在 `beamlib` 命名空间内
- **C++ 标准**: C++17，不使用 C++20 特性
- **依赖**: 仅 Eigen3 (header-only)，无其他第三方库
- **构建**: CMake ≥ 3.16
- **测试**: 纯 `main()` + `ctest`，不用 gtest/catch2
- **代码风格**:
  - 不写注释，除非解释 WHY（非 WHAT）
  - 不加 error handling / fallback，除非在系统边界（如 LinearSolver 奇异检测、SolverStatus 报告）
  - 优先 header-only，只有复杂实现才用 .cpp
  - 变量命名遵循物理约定（E, G, rho, Iy, Iz 等）

## 5. 架构原则

### 5.1 编译期多态 + 运行期 API 边界

内部实现使用模板（编译期 dispatch，零开销）：
```cpp
template <typename ElemType>
class BeamModel { ... };

template <typename ElemType>
class NewtonRaphsonSolver { ... };
```

在 API 边界提供 type-erased wrapper，供 Phase 2 GUI 和 Phase 3 Python bindings 使用：
```cpp
class IAnalysisModel {   // 运行时多态，虚函数接口
public:
    virtual void solve() = 0;
    virtual AnalysisResult getResult() const = 0;
    // ...
};
```

**测试原则**: 数值测试直接用 `BeamModel<ElemType>`，不经过 `IAnalysisModel`。虚函数层仅在 Batch 10 做 API 封装。

### 5.2 结构化结果数据模型

所有分析结果通过 `AnalysisResult` 结构返回，便于程序化访问：
```cpp
struct SolverStatus {
    bool converged = false;
    int iterations = 0;
    double finalResidual = 0.0;
    std::string errorMessage;      // 系统边界的 error reporting
};

struct AnalysisResult {
    std::vector<NodalResult> nodalResults;
    std::vector<ElementForces> elementForces;
    std::vector<ReactionForce> reactions;
    SolverStatus status;
};
```

`SolverStatus` 须包含结构化的错误信息（如 SparseLU 奇异报告），供 GUI/Python 层消费。

### 5.3 可扩展分析类型

分析类型不硬编码，预留接口以后扩展：
```cpp
enum class AnalysisType {
    Static,
    TransientImplicit,
    Modal,
    // future: Buckling, FrequencyResponse, ...
};
```

## 6. 目录约定

```
BeamLib/
├── CMakeLists.txt
├── PROJECT_SPEC.md
├── EXECUTION_PLAN.md
├── REVIEW_NOTES.md
├── docs/
│   ├── theory/                     # 理论推导手册 (LaTeX → PDF)
│   │   ├── 01_euler_bernoulli_2d.tex
│   │   ├── 02_euler_bernoulli_3d.tex
│   │   ├── 03_timoshenko_2d.tex
│   │   ├── 04_timoshenko_3d.tex
│   │   ├── 05_geometrically_exact.tex
│   │   └── 06_time_integration.tex
│   └── tutorials/                  # 应用案例教程 (LaTeX → PDF)
│       ├── tutorial_overview.tex   # 总览 + 软件基本操作
│       └── cases/                  # 各案例独立文件
│           ├── case_01_xxx.tex
│           └── ...
├── include/BeamLib/
│   ├── Core/                       # 基础类型
│   │   ├── Types.h
│   │   ├── Node.h
│   │   └── SectionProperties.h
│   ├── Math/                       # 数学工具
│   │   ├── BeamMath3D.h
│   │   └── Rotation2D.h
│   ├── Element/                    # 单元
│   │   ├── ElementBase.h
│   │   ├── EulerBernoulli2D.h
│   │   ├── EulerBernoulli3D.h
│   │   ├── Timoshenko2D.h
│   │   ├── Timoshenko3D.h
│   │   ├── GeomExact2D.h
│   │   └── GeomExact3D.h
│   ├── Model/                      # 模型和装配
│   │   ├── DofMap.h
│   │   └── BeamModel.h
│   ├── Solver/                     # 求解器
│   │   ├── LinearSolver.h
│   │   ├── NewtonRaphson.h
│   │   ├── ImplicitTimeIntegrator.h
│   │   └── ModalSolver.h
│   ├── Load/                       # 载荷
│   │   ├── LoadTypes.h
│   │   └── LoadManager.h
│   ├── PostProcess/                # 后处理
│   │   ├── InternalForces.h
│   │   └── AnalysisResult.h
│   └── API/                        # Type-erased API 边界
│       └── AnalysisModel.h
├── src/                            # 仅复杂实现
│   └── GeomExact3D.cpp
└── tests/
    ├── test_compile_check.cpp
    ├── test_eb2d_cantilever.cpp
    ├── test_eb2d_postprocess.cpp
    ├── test_eb3d_cantilever.cpp
    ├── test_timo2d_cantilever.cpp
    ├── test_timo3d_cantilever.cpp
    ├── test_ge3d_cantilever.cpp
    ├── test_ge3d_large_deformation.cpp
    ├── test_modal_cantilever.cpp
    └── test_implicit_dynamics.cpp
```

## 7. 元素类型接口约定

每个元素类型是一个 struct，提供 static constexpr 和 static 方法：
```cpp
struct SomeElement {
    static constexpr int nDofsPerNode = N;
    static constexpr bool hasTransformation = true/false;
    static constexpr bool isLinear = true/false;
    
    // 必须: 计算残差和切线刚度
    static ElementResult<N> computeElement(
        const Vec3& xA, const Vec3& xB,
        const VecN<2*N>& dispVec,
        const SectionProperties& props);
    
    // hasTransformation=true 时必须: 坐标变换矩阵
    static MatMN<2*N, 2*N> computeTransformation(
        const Vec3& xA, const Vec3& xB,
        const Vec3& refVector);
    
    // 动力学/模态需要: 一致质量矩阵
    static ElementMassResult<N> computeMass(
        const Vec3& xA, const Vec3& xB,
        const SectionProperties& props);
    
    // 后处理需要: 提取单元内力 (原始数据，两端各一组)
    static ElementInternalForces computeInternalForces(
        const Vec3& xA, const Vec3& xB,
        const VecN<2*N>& dispVec,
        const SectionProperties& props);
};
```

**残差符号约定**:
- 线性单元: `re = ke * dispVec` (内力 = 刚度 × 位移)
- 非线性单元: `re` 从应变能变分直接得出（内力向量）
- NR 求解器: `R_int - F_ext = 0`，校正方程 `K δ = -R`

## 8. 核心数据结构约定

### 8.1 SectionProperties

```cpp
struct SectionProperties {
    double E = 0.0;
    double G = 0.0;
    double rho = 0.0;
    double A = 0.0;
    double Iy = 0.0;
    double Iz = 0.0;
    double Ix = 0.0;
    double kappa_y = 5.0 / 6.0;   // y 方向剪切修正系数
    double kappa_z = 5.0 / 6.0;   // z 方向剪切修正系数

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

    // 2D 便捷: 单一 kappa (取 kappa_z，xz 平面)
    double kappa() const { return kappa_z; }
};
```

### 8.2 ElementConn

```cpp
struct ElementConn {
    int nodeA = 0;
    int nodeB = 0;
    Vec3 refVector = Vec3(0, 1, 0);  // 用户可覆盖截面方向
};
```

3D 变换矩阵用 `refVector` 构建局部 y 轴。当梁平行于 `refVector` 时自动 fallback 到 `(0, 0, 1)`。用户可 per-element 指定。

### 8.3 LinearSolver

```cpp
struct SolveResult {
    VecX x;
    bool success = false;       // SparseLU 奇异时为 false
};

struct LinearSolver {
    static SolveResult solve(const SpMat& K, const VecX& rhs);
};
```

`SparseLU` 必须报告奇异分解（机构或约束不足），不能静默返回垃圾数据。

## 9. DOF 顺序约定

- **2D (nDofsPerNode=3)**: `[u_x, u_z, θ_y]` — xz 平面内
- **3D (nDofsPerNode=6)**: `[u_x, u_y, u_z, θ_x, θ_y, θ_z]`
- **单元 DOF**: `[nodeA DOFs, nodeB DOFs]`

## 10. 坐标约定

- 梁的参考方向为局部 x 轴
- 2D 问题在 xz 平面内（x 水平，z 竖直向上）
- 2D 变换: cos/sin 取自 dx_x 和 dx_z 分量
- 3D 变换: 用 `ElementConn::refVector` 构建局部 y 轴；默认 (0,1,0)，当梁平行于此向量时 fallback 到 (0,0,1)
- 用户可 per-element 覆盖 `refVector` 以定义任意截面方向

## 11. 数值精度要求

分级容差 (不同测试类型使用不同精度要求):

| 测试类型 | 容差 | 说明 |
|----------|------|------|
| 单元素矩阵检查 / patch test | 相对 1e-12 ~ 1e-10 | 直接矩阵元素对比 |
| 多单元解析位移 (线性) | 相对 1e-9 ~ 1e-8 | 含装配和求解误差 |
| 支反力 (大值) | 相对 1e-9 | 与解析解对比 |
| 支反力 / 近零 DOF | 绝对 1e-12 ~ 1e-10 | 不用相对容差 |
| GE 非线性小载荷 | 相对 1e-4 | 与 Timoshenko 解析解对比 |
| GE 大变形 | NR 收敛 + 物理合理 | 必须含一个定量文献基准 |
| 隐式动力学 | 无阻尼能量漂移 < 1% | Generalized-α 数值耗散可控 |
| 模态分析 | 前 3 阶频率相对误差 < 1e-6 | 与解析频率对比 |

**"1 次 NR 迭代" 定义**: 初始残差评估 → 一步校正 → 残差降至 ~machine epsilon。线性问题必须满足。

### 约束处理

- **Homogeneous Dirichlet**: `node.fixed[i] = true`，DOF 不进入方程组。Batch 2A 即支持。
- **Nonzero prescribed displacement**: elimination + RHS correction。Batch 2B 完成设计、实现和测试（Batch 4 硬性前置）。
- **Phase 1 不使用** Lagrange multiplier 或 penalty 方法。

## 12. 文档体系

Phase 1 配套两类文档，均为 LaTeX → PDF：

### 12.1 理论推导手册 (`docs/theory/`)

**每种梁理论在编码前必须完成**。

每篇结构:
1. 运动学假设与坐标系
2. 应变-位移关系推导
3. 虚功原理 → 弱形式
4. 有限元离散 (形函数选择及理由)
5. 单元刚度矩阵完整推导
6. 单元质量矩阵推导 (含旋转惯量决定)
7. 坐标变换推导
8. 已知陷阱 (如 Timoshenko shear locking)
9. 验证方案设计

**每章必须包含三个附表**:
- **符号约定表**: 正方向定义、DOF 编号、力/位移正方向
- **公式 → C++ 索引映射表**: MATLAB 1-based → C++ 0-based 的完整映射
- **验证容差表**: 该理论对应的所有测试及其容差

**关键原则**: MATLAB 代码是参考而非 Golden Standard。所有公式必须从第一性原理推导，不依赖参考代码的正确性。

### 12.2 应用案例教程 (`docs/tutorials/`)

面向用户（学生/工程师），展示 BeamLib 在真实工程分析场景中的使用。

**教程总览** (`tutorial_overview.tex`):
- BeamLib 能力介绍（支持的单元/分析类型/后处理）
- 软件基本操作流程: 定义材料/截面 → 建模(节点/单元) → 施加约束/载荷 → 求解 → 结果提取
- API 快速上手 (C++ 示例 / 便捷层)

**案例文档** (`docs/tutorials/cases/`):
- 每个案例一份独立 `.tex`
- 结构: 工程背景 → 力学建模 → BeamLib 建模代码 → 求解 → 结果分析 → 与手算/理论对照

**首批案例 (Phase 1 优先)**:
1. EB2D 悬臂梁 / 简支梁 — 最基础入门
2. 2D 门式框架侧向荷载 — 多构件组装
3. Timoshenko 深梁 vs EB 细长梁对比 — 理论差异可视化
4. 3D 机械臂末端刚度 — 3D 应用

大变形和模态教程等线性流程稳定后再写。

**编写时机**: 案例教程在对应 Batch 的代码完成且通过验证后编写。

```
理论文档 → 编码 → 验证 → 案例教程
```

## 13. 验证方案 (独立于参考代码)

### 13.1 解析解验证
- 悬臂梁: 集中力 → 端部位移/转角
- 悬臂梁: 分布力 → 端部位移 (等效节点力)
- 简支梁: 均布载荷 → 跨中挠度
- 两端固支梁: 集中力 → 位移和反力
- **纯端弯矩**: EB2D/3D 线性弯矩分布验证
- Timoshenko 含剪切解析解: `u_z = FL³/(3EI) + FL/(κGA)`

### 13.2 Patch Tests
- 常应变 (纯轴向): 所有单元应给出精确解
- 纯弯曲: EB 单元应精确，Timo 单元需验证
- 刚体运动: 零应变/零内力
- **单元素 θ_y 符号测试**: F_z 下 θ_y 的正负号

### 13.3 收敛性研究
- 网格加密 (h-refinement): 位移误差应按理论阶次收敛
- EB: O(h²) for displacements, Timo: 需验证 locking-free

### 13.4 Shear Locking 测试
- **Timo 粗网格 + 极细长梁** (L/h = 1000, 2~5 单元): 不应 locking
- Timo 细长极限 → EB 结果一致

### 13.5 文献基准
- 大变形: **Bathe & Bolourchi 45° cantilever bend** (1979) — **必须定量通过**
- 备选: Simo-Reissner cantilever roll-up
- 模态: 均匀梁解析频率 ω_n = (β_n L)² √(EI/ρAL⁴)

### 13.6 交叉验证
- GE 小载荷 → 与 Timo 结果一致
- Timo (L/h → ∞) → 与 EB 结果一致
- 2D 结果 → 与对应 3D 面内结果一致
- 模态正交性: `Φᵀ M Φ = I`

**2D/3D 交叉验证截面属性映射表**:

| 2D 属性 | 对应 3D 属性 | 验证场景 |
|---------|-------------|---------|
| Iz (2D xz 弯曲) | Iy (3D xz 弯曲) | F_z 悬臂梁位移 |
| kappa (2D) | kappa_z (3D) | Timoshenko 剪切修正 |

## 14. 参考代码（只读，不修改）

- 已有 C++: `C:\_ZW\DEM\01_res_and_dev\FSI_dev\DEM-FEM Coupling\beam_fem\`
  - 3D 几何精确梁完整实现（14 个测试通过）
  - BeamMath.h (Rodrigues rotation), BeamElement3D.cpp, BeamModel.cpp, BeamSolver.cpp
- MATLAB: `C:\_ZW\DEM\01_res_and_dev\FSI_dev\JaphyBeam\`
  - EB2D, EB3D 完整实现，变换矩阵

## 15. 协作模式

- **Claude Code**: 规划、架构设计、理论文档、主要代码实现
- **Codex/GPT**: Review 数值正确性、边界条件、风险点检查
- 每个 Batch 完成后由 Codex review
- 文档优先于聊天记录；有分歧时以仓库内文档为准

## 16. Git 分支策略

```
main                    ← 仅 review 通过后合并
├── feat/batch-1-scaffold
├── feat/batch-2a-eb2d-core
├── feat/batch-2b-eb2d-postprocess
├── feat/batch-3-eb3d
├── ...
```

先创建 baseline commit，再按 batch 开 feature branch。
