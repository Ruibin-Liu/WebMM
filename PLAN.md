# Plan: 接入 pi/magpie GLM Provider（本地代理 / Schema 重试 / 实跑验收）

用户裁决「就用我们 pi 里有的」= 复用本机 pi/magpie `127.0.0.1:3425` 与 `zcode/GLM-5.3-Flash`，不引入模型账号或密钥到 WebMM；增加 localhost 开发代理和真实 Provider，保留本地草稿降级。

## 交付
1. **localhost 开发服务器（scripts/agent-dev-server.mjs，无新依赖）**
   - 仅绑定 127.0.0.1，服务仓库静态文件与 `/api/agent/status`、`/api/agent/chat`；默认 :8902，可用 AGENT_DEV_PORT 改。
   - `/v1/chat/completions` 上游固定 MAGPIE_URL（默认 pi :3425）与 AGENT_MODEL（默认 zcode/GLM-5.3-Flash），请求体白名单/1 MiB/超时，非流式，不记录分子 payload。
   - package.json 增 `agent:dev` 脚本；不改现有 :8901 或 magpie 服务。
2. **Magpie Provider（app/agent/magpie-provider.js）**
   - `status()` + `plan(request)` + `review(request)`，fetch 可注入以便 Node 测试。
   - 严格 JSON 提示、OpenAI content 解析/代码围栏恢复；Plan 用 Plans.validatePlan，Review 用 Evidence.validateReview。
   - Schema 失败仅纠正重试一次（回传具体错误）；二次失败抛 PROVIDER_SCHEMA，不执行；网络失败不盲目重试。
3. **Plan Card 显式连接**
   - 新增“连接 pi GLM”：status 成功才设置 window.WebMMAgentProvider；展示模型名和“上下文经 pi/magpie 路由”。
   - 未连接保持本地草稿/本地证据摘要；连接后沿用现有双调用、编辑、user:local 批准、真实执行和引用复盘。
4. **测试与文档**
   - Node 增 Provider 单测：status、合法 JSON、围栏恢复、Plan/Review 首次错误后纠正、二次失败、无执行方法、payload 白名单。
   - 新增默认不注册的 agent_glm_live CDP：连接真实 pi/magpie → GLM Plan → 编辑小预算 → 真 analogScan → GLM Review；仅手动 WEBMM_AGENT_URL 运行。
   - 原 agent_plan/mock、m12_tasks、Node suites 零回归；docs 记录启动方式、隐私边界和 GLM 实测。

## 验收
- `node tests/agent/run.js`、`node tests/platform/run.js` 全绿；静态/代理单测覆盖路径穿越与上游模型不可由浏览器改写。
- `tests/cdp/run.sh agent_plan` 10/10，`m12_tasks` 14/14 零回归。
- 手动 live：真实 GLM 计划可校验，批准后 task 带 planHash/approvedBy，复盘引用真实 evidenceId；若 GLM 输出非法必须经一次纠正或诚实失败，不得执行非法计划。
- node --check / bash -n / git diff --check 全过；无新 npm 依赖、无密钥入仓、不重启 magpie；推送分支。

## 排雷存档
- pi/magpie :3425 的 OPTIONS 不支持浏览器预检；未改 magpie，新增仅 localhost 的静态+API 同源代理，模型由服务端固定，浏览器不能覆盖。
- GLM 首次真实 Plan 合法并执行成功；Review 空 content 是长推理耗尽 4096 tokens，简单复现显示 reasoning_content 与最终 content 共用预算。Review 上限升 8192 后 live 7/7，空响应错误现带 finish_reason/completion_tokens。
- 真实模型谱系不能只留 planHash：approval 增加 modelProvider/model/providerRevision，随 TaskSpec 与 Evidence 持久；本地草稿记 null。
- 原 :8901 与 magpie :3425 均未重启；测试只启动/停止本 worktree 的 :8902 agent:dev。live 默认不注册 run.sh，避免 CI 隐式消耗模型。
