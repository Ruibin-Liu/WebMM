# Plan: 第二个 LBDD Skill（Bioisostere / aza-deaza Exploration）

用户确认继续 = 增加真实第二个 Skill，但范围限定为现有芳环 aza/deaza 扫描，不接全部 479 条反应库；泛化 Skill→Capability→PlanKind→Enumerator 映射，R-group 全链路零回归。

## 交付
1. **Capability / Skill 泛化**
   - Capabilities 增 `azaScan` revision 1（taskKind=analogScan、planKind=bioisostereExploration、fragmentSet=builtin-aza-v1、category=aza|aza+deaza、真实证据与限制）。
   - Skills 增 `bioisostere-exploration` v1；Orchestrator/Provider 由 skill 决定 capability 与动态提示，不再硬编码 R-group。
   - DesignPlan 按 capability 校验 kind/skill/fragmentSet/category；provenance/Evidence 带 capability；Decision.next 接受两种 planKind。
2. **Headless aza 枚举（app/agent/aza.js + bio-host.js）**
   - 抽取现有 C→N 单/双位点、可选 N→C 单/双位点、规范 SMILES 去重、描述符/SA、__azaCounts 和 RDKit 对象释放；legacy runAzaScan 共用。
   - BioHost 使用 siteRef=aromatic-ring:auto、方向 category 和母体 revision；返回与 R-group 相同的 candidates/audit/runCtx/execution context，批准与任务执行不复制。
3. **Plan Card 与历史接入**
   - 新增 Skill 选择（R-group / Bioisostere），动态切换 host/tools/provider skill；批准、运行、复盘、Decision 和复制均使用当前 plan 的 skill/capability。
   - 历史显示 planKind/capability；复制按当前 plan skill 取对应 host，不允许跨 skill/上下文误套。
4. **测试与文档**
   - Node 增 aza capability/plan 校验、动态 provider、Decision next、Evidence limitations；R-group 既有断言全保留。
   - CDP 新 bio_agent_plan：mock bio Provider → 小预算批准 → 真 aza task → evidence/review → next draft；legacy aza/deaza 与 R-group agent_plan 零回归。
   - docs 明确这只是 aza/deaza 第二 Skill，不代表通用 bioisostere 反应库。

## 验收
- Node agent/platform/proxy 全绿；R-group CDP agent_plan 18/18；新 bio CDP 覆盖 aza 单/双替换、deaza 方向、spec capability/planHash/budget、Evidence/Review/Decision。
- m12_tasks 14/14、真实 GLM live 8/8；无 page error、无新依赖、无自主执行。
- node --check/bash -n/git diff --check 全过并推送分支。

## 排雷存档
- 显式 H 的 3D 母体可产生显式 H SMILES；直接替换原子会错位。Aza 先去 H 并 canonicalize，再**重新解析 canonical SMILES** 生成 surgery base——get_molblock 保留输入序，仅 canonicalize 不重解析仍不够。
- 化学验收不能只查 SMILES 含 N（原分子酰胺 N 也会误过）；bio CDP 要求真实芳环 n 生成物，legacy deaza 要求 N-free 产物。
- BioHost 初始化未写 __azaCounts/__azaRun，已补；完成时 __azaRun 会复位，断言使用 runCtx.mode 与 counts。
- legacy 完成状态在 Done 文本后仍可能 running，测试改为等待无 planHash task=complete，避免提交事务竞态误判。
- Orchestrator/MagpieProvider 不再硬编码 R-group，改由 skill 决定 capability/planKind/提示；两种真实 GLM live 各 8/8。
