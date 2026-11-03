# Plan: LBDD AI 协同核心契约（Capability / Plan / Tools / Review）

用户裁决「开始实现这个想法」= 在独立 worktree 建立 AI 协同第一层，不接模型商、不改 UI：让 R-group 设计计划可校验、经明确批准可编译为现有 analogScan TaskSpec，并以可引用证据支持受控复盘。

## 交付
1. **Capability Registry + R-group Skill（app/agent/）**
   - UMD 纯模块，浏览器与 Node 共用，无 DOM/I/O/新依赖。
   - 首个 capability `analogScan`：revision、输入、证据、预算和科学限制；注册表是模型可见能力的唯一来源。
   - `rgroup-optimization` Skill：所需上下文、允许目标、分级计算顺序、活性/跨分子能量等诚实边界。
2. **DesignPlan + approval-bound PlanCompiler**
   - `rgroupOptimization` v1：kind/capabilityRevision/baseRevision/parentRef/siteRef/fragmentSet/goals/unverifiedGoals/budget/hypothesis。
   - 规范化预算（默认 400/30/10，硬顶不超现有漏斗），未知目标拒绝；实验活性等未支持目标必须显式列入 unverifiedGoals。
   - planHash 绑定完整规范计划；编译时重验 parent/site/baseRevision、候选数与 approval{planHash,baseRevision,approvedBy}，任何不匹配均拒绝。
   - 复用 Tasks.Contract.makeSpec；仅调整 A/B/C/D 确定性预算，保留 nConfs/seedOffset；计划/批准/未验证目标写 provenance，不改冻结 input v1。
3. **受控 Tools + EvidenceBundle / DesignReview**
   - 模型侧仅提供 getDesignContext/listCapabilities/getCapability/getLibrarySummary/validatePlan/getEvidence；不提供 execute。
   - 从真实 spec + fact 构造带 evidenceId 的证据包，保留 accounts、truncated、科学限制和 unverifiedGoals。
   - DesignReview claims 必须引用包内真实 evidenceId；未知引用拒绝，unresolved 保留，nextPlanDraft 再走独立计划校验。
4. **模型无关 Orchestrator**
   - 注入 `plan(request)` / `review(request)` 两个 provider 方法；规划与复盘各一次受控调用，不做自主循环。
   - proposePlan 只返回 awaiting_approval/rejected；reviewResults 校验证据引用和下一轮草稿，不执行计划、不替用户做 Decision。
5. **文档与测试**
   - docs/ai-collaboration.md 冻结首版边界、模块职责、工具清单、批准与复盘契约、后续浏览器/模型接入点。
   - tests/agent/run.js 覆盖注册表、Skill、规范化和哈希、未验证目标、预算、批准、陈旧上下文、TaskSpec 兼容、证据引用和双调用编排。
   - tests/cdp/run.sh 增加 Node agent suite；既有平台 Node 套件零回归。

## 验收
- `node tests/agent/run.js` 全绿，且包含：未声明 experimentalActivity 被拒、超 400/30/10 预算被拒、无/错批准被拒、陈旧 baseRevision 被拒、编译 spec 的预算/provenance/inputRevision 正确、伪造 evidenceRef 被拒、Orchestrator 无执行旁路。
- `node tests/platform/run.js` 既有 153/153 绿；不改 Rust、WASM export、TaskSpec input schema、UI 或依赖。
- `git diff --check` 通过；CODE_STATUS 记录真实完成与未做项，并推送 `feat/lbdd-ai-copilot`。

## 排雷存档
- 初稿把 conformerEnergy 列为 analogScan 证据，但真实结果只有形状/柔性分数与候选描述符；已删除该宣称，Capability 只暴露实际产物。
- 规划上下文若直接透传会把 parentSdf、候选 molblock、buildId 带给远程模型；Tools 改白名单投影，只保留设计身份、目标、audit 和 candidateCount，编译另用本地执行上下文。
- 候选原样 clone 可能把额外字段带入冻结 input v1；编译器改为 idx/smiles/molblock/mw/rotb/sa 显式投影并校验有限数值。
- Review 对非数组 counterexamples 可能在校验时抛类型错误；改为所有畸形输入都返回错误，测试覆盖。
