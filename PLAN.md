# Plan: LBDD AI 协同浏览器竖切（Headless Host / Approved Run / Plan Card）

用户裁决「看 GLM 能不能用，按计划下一步」= 先实测 GLM 结构化规划/复盘，再把 headless 契约接入 platform.html：复用现有 R-group 枚举与 analogScan 执行，提供可编辑、批准、运行和证据复盘的最小 Plan Card；不接真实模型密钥。

## 交付
1. **GLM 可行性记录**
   - 通过 Sider `glm-5.3-flash` 实测：DesignPlan 首次合法；DesignReview 首次越界字段/推断被 validator 拒绝，纠正提示后合法。
   - docs 记录“结构化输出 + 确定性校验 + 失败重试”要求；Sider 仅测试，产品仍需独立 GLM endpoint/代理。
2. **Headless R-group Host（app/agent/rgroup.js + browser-host.js）**
   - 抽取 platform.html 现有 A 阶段为无 DOM 枚举器，保留价键/重原子/转子过滤、显式 H 清理、SA、审计与 RDKit 对象释放；旧探索按钮改调同一实现。
   - browser host 提供当前身份/执行上下文、库摘要、计划枚举和真实 task+fact 证据；project/round/buildId 沿用现有服务。
   - DesignPlan 增 `fragmentCategory`（默认 all），与 parent/site/fragmentSet 共同绑定 baseRevision/planHash；Tools 白名单同步。
3. **批准计划接现有执行链**
   - runAnalogFromCands 接受可选 {plan,approval}：有则 compileApprovedPlan，无则保留原 makeSpec；token、池、scheduler、task store、Done fact、取消语义完全复用。
   - 返回 taskId 供复盘；legacy 路径不写 plan provenance，agent 路径带 planHash/approvedBy/fragmentCategory。
4. **最小 Plan Card（platform.html，不重设计 UI）**
   - 本地草稿（或可选 window.WebMMAgentProvider）→ hypothesis/预算/活性未验证编辑 → 批准并枚举运行 → 证据复盘。
   - provider 存在时用既有 Orchestrator plan/review 双调用；无 provider 时明确显示本地草稿/本地证据摘要，不伪装 AI。
   - 模型文本一律 textContent；批准前重验当前分子/位点/类别，变化则拒绝。
5. **测试与文档**
   - 更新 Node agent 契约（fragmentCategory、上下文脱敏）；新增 CDP agent_plan：headless tools、Plan Card 编辑批准、spec 预算/provenance、真实 evidence、复盘、legacy 回归、零 page error。
   - docs/ai-collaboration.md 更新浏览器接入与 GLM 实测；run.sh 注册 agent_plan；CODE_STATUS 记录实际完成和未做。

## 验收
- Node：tests/agent、tests/platform 全绿；新增字段进入 hash/provenance，执行 payload 仍不越过模型边界。
- CDP：`tests/cdp/run.sh agent_plan` 全绿，agent 运行候选≤批准预算、C/D caps 正确、fact 可复盘，legacy run 无 planHash；另跑 `m6_analog_ui` 或 `m12_tasks` 证明旧路径零回归。
- `node --check` 新/改 JS、`bash -n tests/cdp/run.sh`、`git diff --check` 全过；不改 Rust/WASM export/依赖。
- 推送 feat/lbdd-ai-copilot；真实 GLM endpoint、密钥 UI、MCP 和自主循环明确留待后续。

## 排雷存档
- 新 worktree 起初缺 gitignored `pkg/webmm.js`，CDP 首跑引擎 404；在本 worktree 执行 `npm run build` 后通过，不是前端逻辑回归。
- `fragmentCategory` 是实际枚举语义，初版 Plan 未绑定会让批准跨类别复用；已加入 planHash/baseRevision/provenance/context 校验，默认 all。BrowserHost.enumerate 必须读当前类别而非用计划类别覆盖，CDP 钉住 STALE_CONTEXT。
- Agent 编译所需 idx 由 runAnalogFromCands 统一补齐；执行上下文候选以该映射为准，编译器再投影冻结 input v1 字段。
- GLM 首次 Review 自由增加字段和推断，证明不能把模型文本直接当证据；只保留 Schema + evidenceRef 校验 + 纠正重试路径。
- 原 :8901 服务属于原工作区，测试未停止它；Agent 用 :8902 临时服务器完成，m12_tasks 新增 WEBMM_CDP_URL 覆盖以便同 worktree 回归。
