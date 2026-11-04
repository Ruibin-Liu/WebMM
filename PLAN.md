# Plan: Plan / Review / Decision 项目级持久化与下一轮草稿

用户裁决「继续下一项」= 将 AI 协同产物从 TaskSpec provenance 提升为现有项目命令日志的一等记录，支持重载/导出、人工决策与下一轮计划草稿；不新增数据库或第二 Skill。

## 交付
1. **State 命令与投影（app/platform/state.js，schema v1 向前兼容）**
   - 新增 designPlans/designReviews/designDecisions 及 hidden 覆盖层；旧快照缺省为空，旧应用重放未知命令沿用跳过策略。
   - RecordDesignPlan(draft) → ApproveDesignPlan → CompleteDesignPlan(taskId) → RecordDesignReview → RecordDesignDecision；引用、状态、重复记录和 task/plan 匹配由 reducer 拒绝。
   - Decision action=accept/reject/needsExperiment/next，actor/reason/nextPlanDraft 可选语义明确；hide/unhide 支持 Decision/Review 撤销，Approve 在有完成/复盘/决策依赖时拒绝 undo。
2. **BrowserHost 持久化接口**
   - 注入现有 project.apply，提供 recordPlanDraft/Approval/Complete/Review/Decision 与 artifacts 读取；模型 Tools 不增加写权限。
   - Plan/Review 保存 provider/model/providerRevision；Decision 永远 actor=user:local，AI 不替人做 Decision。
3. **Plan Card 接线**
   - 有效草稿落 RecordDesignPlan，批准落 Approval，任务完成落 Complete，真实或本地复盘落 Review；保存失败明确提示且不伪称成功。
   - 新增 accept/reject/needsExperiment/next 决策与 reason；next 保存含 nextPlanDraft 的 Decision 后生成新的待批准草稿，不自动执行。
   - 项目导出/重载经既有命令日志和快照恢复，不改 format.version。
4. **测试与文档**
   - Node 覆盖状态机、重复/错误引用、hide/unhide、依赖拒绝、序列化与 replay；保留 153 项既有测试。
   - CDP agent_plan 增 Plan/Review/Decision 投影、下一轮新 draft、导出命令和重载恢复；live 仅模型谱系不受影响。
   - docs/CODE_STATUS 更新真实边界；run.sh 不注册 live。

## 验收
- Node agent/platform 全绿：Plan 不能跳状态、Review 必须指已完成 task、Decision 必须 user actor、next 必须带有效 nextPlanDraft、undo 不删事实。
- CDP agent_plan 全绿：draft→approved→complete→review→next decision 全入项目，导出含新命令，重载后记录仍在；legacy 与零 page error 保持。
- m12_tasks 14/14、agent_glm_live 7/7（若模型服务可用）无回归；node --check/bash -n/git diff --check 全过。
- 无新依赖、无新数据库、无自主执行；推送 feat/lbdd-ai-copilot。

## 排雷存档
- 批准前编辑会改变 planHash：先在 UI 保存初稿、再保存批准前规范版本，不能只对旧 hash 调 Approve。CDP 实证应保留初稿/编辑版/下一轮三个 RecordDesignPlan，而非两个。
- Decision 必须引用已完成的旧 plan+task；nextPlanDraft 另存新 hash 和新 draft，不能继承旧批准或自动执行。
- Review/Decision undo 只走 hidden overlay；Complete/依赖后的 Approve 不允许 undo，保持追加日志与依赖语义。
- schema v1 未升：旧快照缺省 design 集合，旧应用重放未知命令仍按既有跳过策略；导出自动包含新命令。
- 测试初稿误把 Node 侧不存在的 agentLastTaskId 当变量、误算计划版本数；分别改为 page.evaluate 返回与三版本语义，最终 agent_plan 15/15、live 8/8。
