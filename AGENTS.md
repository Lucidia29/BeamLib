# BeamLib 协作入口

## 先读这些文件
1. `PROJECT_SPEC.md`: 项目总规范、硬约束、架构原则、数值基准
2. `EXECUTION_PLAN.md`: 分阶段执行总计划、每 batch 的接口/验收标准
3. `REVIEW_NOTES.md`: 审阅意见、变更点、待确认问题

## 协作分工
- **Claude Code (CC)**: 规划、架构设计、理论文档、主要代码实现、自测
- **GPT/Codex**: Review 数值正确性、理论推导、边界条件、风险点检查

## 工作流程
```
每个 Batch:
  1. CC 完成理论文档 (LaTeX PDF) — 如该 batch 涉及新理论
  2. GPT Review 理论文档
  3. CC 实现代码 + 编写测试
  4. GPT Review 代码（重点: 数值、符号、边界条件）
  5. 修订写入 REVIEW_NOTES.md
  6. 通过后合并到 main
```

## 关键原则
- 文档优先于聊天记录；有分歧时以仓库内文档为准
- MATLAB/C++ 参考代码是辅助理解，不是 Golden Standard
- 理论文档在编码前完成，从第一性原理推导
- 验证方案独立于参考代码（解析解、patch test、收敛性、文献基准）

## GPT/Codex Review 聚焦点
- 数值方法假设是否正确
- 边界条件处理是否完备
- DOF 索引映射（MATLAB 1-based → C++ 0-based）
- 刚度/质量矩阵对称性
- 残差符号约定一致性
- Shear locking 风险
- 坐标变换顺序和方向
