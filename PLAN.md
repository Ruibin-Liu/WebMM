# Plan: LBDD 项目历史视图（Plan / Review / Decision / Lineage）

用户确认「好」= 继续做下一项：让已持久化的设计产物在 analog 页面可见、可检查轮次关系，并可从旧计划显式复制为当前上下文的新草稿；不改执行语义。

## 交付
1. **设计历史 UI（platform.html）**
   - 在 AI Plan Card 下新增可折叠设计历史：计划状态、短 planHash、假设、片段类别、模型、taskId、createdAt。
   - 按 planHash 分组显示 Review claims/unresolved 与 Decision action/actor/reason；next 显示旧 plan → 新 plan 的链接。
   - 全部用 DOM/textContent，模型内容不进入 innerHTML；隐藏覆盖层继续由 getArtifacts 过滤。
2. **刷新与复制新草稿**
   - draft/approval/complete/review/decision 后及进入 analog 模式时刷新；重载后无需恢复旧 Plan Card 也能看到项目历史。
   - `agentPlanCopyDraft(planHash)`：仅在当前 parent/site/category 与旧计划一致时，更新 baseRevision、复制假设并生成独立待批准草稿；不复制 approval/taskId，不自动执行。
3. **测试与文档**
   - CDP agent_plan 增加历史分组、next 链接、重载后可见性、复制为第四草稿且任务数不变；既有 15 项零回归。
   - m12_tasks、Node suites、真实 GLM live 按需回归；docs 明确历史语义和复制边界。

## 验收
- CDP agent_plan 全绿：三版本历史可见、Review/Decision 归属正确、重载后仍可见，复制产生第四 draft 且无新 task、无 page error。
- Node agent/platform/proxy 全绿；m12_tasks 14/14；agent_glm_live 8/8（若模型可用）。
- node --check / bash -n / git diff --check 全过；无新依赖、无自主执行、无历史事实改写；推送分支。

## 排雷存档
- 设计历史与 analog 执行历史分开：前者按 planHash 解释设计/证据/人工决策，后者只列 run；next lineage 用规范化的 nextPlanDraft hash 关联。
- 复制旧计划不能盲目 rebase 到另一分子：仅 parentRef/siteRef/fragmentCategory 一致时更新 baseRevision；approval/taskId 不复制，CDP 验证第四 draft 但 task 数不变。
- 重载不恢复旧 Plan Card 内存状态，history 必须直接从 project projection 渲染；getProjectState 在 compat 未就绪时降级空对象，避免启动竞态。
- 所有模型/旧计划文本均 DOM textContent；next/复制后历史即时刷新，进入 analog 模式也刷新。
