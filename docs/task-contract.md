# Task Contract(信封层,冻结)

配套设计:docs/task-executor-design.md(为什么与整体架构;本文是字段级规范)。
冻结范围:**信封 + 状态机 + seed 派生 + 指纹规范化 + 记账 + 落账语义**。
kind 的 input schema 随该 kind 迁移落地时冻结(§9);未冻结的 kind 只能排队不能执行。

## 1. TaskSpec(信封)

```
{
  taskId:   string           // 生成时定,全局唯一;'a'+Date.now().toString(36)+'-'+rand
  kind:     string           // §9 注册表
  version:  1                // 信封版本;input 语义变更=kind inputVersion+1,不改信封
  input:    object           // kind 特定;提交时必须先经 effectiveInput() 展开(§4)
  seed:     uint32 | null    // 任务级随机源(roundSeed);无随机性 kind 为 null
  provenance: {
    projectId: string | null,
    roundId:   string | null,
    molId:     string | null,
    buildId:   string,       // SHA-256(wasm 字节) 前 16 hex;RDKit 侧记 rdkitVersion
    policyRev: string | null // 漏斗策略内容指纹(analogScan 必填,见 §4)
  }
}
```

## 2. 状态机

```
queued → running → complete | failed | cancelled | interrupted
queued → cancelled | abandoned
running → abandoned          // 仅重载对账
complete → stale             // 仅结果应用校验失败(§7)
其余为终态
```

- `interrupted`:时间/内存/取消导致的非完整执行。**interrupted 的结果绝不冒充完整排名**;已计算部分可入事实但必须带 `interrupted: true`。
- `abandoned`:重载后无存活 session 认领(§6)。
- `stale`:任务成功完成,但应用时目标世界的 inputRevision 已变(§7)。

## 3. seed 派生(域分离)

```
deriveSeed({roundSeed, kind, stage, molId, confIdx})
  = fnv1a32(canonical({roundSeed, kind, stage, molId, confIdx}))
```

- 所有键参与;缺失键以显式 `∅` 参与(不省略——防不同调用形状碰撞)。
- `stage` 必须是稳定 stageId(如 'B'、'C')或 op 名,**不得用数组位置**。
- 派生值为 uint32;worker 内如需 int 种子直接取整。
- 兼容锚:analogScan B 阶段 per-candidate seed = `seedOffset + idx`(seedOffset 入有效 stages 快照),迁移时保持 4200/14200 以减小分数漂移;此后新增阶段一律走 deriveSeed。

## 4. 参数指纹(规范化序列化)

对**展开默认值后的有效输入**(effectiveInput)计算,不是对声明输入:

1. 对象键**排序**后序列化;`undefined` 值的键**丢弃**;数组**保序**(顺序是语义)。
2. 数字:IEEE754 位串(16 hex,大端);`NaN/+Inf/-Inf/-0` 各自字面量(`-0 ≠ 0`,位串不同)。
3. 字符串 JSON 串;`null` 字面量;boolean 字面量。函数/symbol/循环 = 契约错误(抛错,不静默)。
4. 指纹 = fnv1a32(canon({kind, inputVersion, effectiveInput, seed}))。
5. `policyRev`(analogScan)= fnv1a32(canon(effectiveStages))——策略身份独立于参数身份。
6. **rerun 等价判断用指纹;用户界面显示可读差异,永不显示指纹本身。**
7. 旧任务 rerun 用**提交时快照的 effectiveInput**,绝不静默套用当前默认。

## 5. 阶段记账(科学诚实,每阶段必填)

```
{ evaluated, filtered, failed, notEvaluated }   // 4-way 计数,入 Done fact
```

- `evaluated`:真算了;`filtered`:确定性规则筛掉(规则入档);`failed`:算了但失败;`notEvaluated`:确定性工作量预算截断未轮到。
- **排名截断只允许确定性预算**(候选数/构象数/阶段 cap);任何时间/内存触顶 → `interrupted`。
- 排名域必须显式:analogScan 的排名域 = evaluated ∪ failed 中有分数者?否——**排名域 = evaluated 且得分数值者**;failed 不参与排名但计数可见。

## 6. 落账语义

- **旁表**(platform DB `tasks` store,upsert by taskId):spec **提交时**写入(state=queued)——唯一 Done 前写入;此后仅状态迁移与耗时/失败原因补写。进度永不入库。
- **Done = 单事务**:command + facts + 旁表 complete 经 `putBatch` 同库跨 store 原子提交;幂等键 `taskId + attemptId + generation`——重复 Done/迟到回调被拒。
- **对账(boot)**:旁表 `queued|running` 的任务:对应 fact 已存在 → `complete`(幽灵完成);无存活 session 认领 → `abandoned`。session 认领 = BroadcastChannel('webmm-tasks') 心跳,响应窗口 400ms;无 BroadcastChannel 环境(旧浏览器/Node)视为单 session,全部 abandoned(除非 fact 在)。
- **写放大红线**:旁表只在 {提交, 状态迁移, Done} 三时点写;facts 每任务一条(阶段摘要内嵌,只含影响结果字段:cap/n/seedOffset/通过数/触顶/4-way 计数)。

## 7. 结果应用校验(inputRevision)

- TaskSpec.input 必须携带 kind 的 `inputRevision`(内容哈希,fnv1a32(canonical)):
  - conformers/optimize:parent molblock 内容哈希 + restraints 指纹
  - shapeAlign:query/target 姿态哈希 + color_scale
  - analogScan:parentSmiles + fragmentSet rev + 有效 stages(即 policyRev)
- 应用前校验:目标当前 revision == spec revision 且请求身份(taskId+generation)匹配;不匹配 → `stale`,UI 提供「对旧结构重算 / 丢弃」。

## 8. 调度与池(语义,非字段)

- 全局并发令牌 = min(4, hardwareConcurrency);background 最多 total−1,interactive 可占满 total。
- 抢占仅在队列层;worker 内不抢占;漏斗是逻辑任务,**其全部 worker 作业逐个持令牌**。
- 硬取消 = terminate + 顶上;协作取消只是优化;迟到 reply 按 `taskId+generation` 丢弃。
- worker 回收:处理满 N 个任务或堆超限 → 空闲时 terminate 重建;每池留 1 预热备用(上限 poolSize+1)。

## 9. kind 注册表

| kind | input schema | 状态 |
|---|---|---|
| analogScan | §9.1 | **冻结**(本轮迁移) |
| conformers | 随迁移冻结 | 排队中 |
| optimize | 随迁移冻结 | 排队中 |
| shapeAlign | 随迁移冻结 | 排队中 |

### 9.1 analogScan input(v1)

```
{
  parentSmiles: string          // 亲本 canonical SMILES
  parentSdf: string             // 亲本 3D(对齐参照;内容哈希入 inputRevision)
  candidates: [{idx, smiles, molblock, mw, rotb, sa|null}]   // A 阶段产物(主线程枚举,豁免路径)
  stages: [ effectiveStages ]   // 提交时已展开默认值的快照,validate 过
  inputRevision: { parent: fnv1a32(parentSdf), fragmentSet: 'builtin-v1', policyRev }
}
```

**effectiveStages v1**(A/B/C/D;字段名即契约):
```
[
  { id:'A', op:'enumerate', caps:{ products:400, rotbMax:10, heavyDeltaMax:+9, heavyDeltaMin:-6 } },   // A 已在主线程完成,cands 即产物;此处为审计快照
  { id:'B', op:'score', nConfs:10, seedOffset:4200, cap:400, inflightPerWorker:2 },
  { id:'C', op:'score', nConfs:30, seedOffset:14200, cap:30, combine:'max' },
  { id:'D', op:'flex', cap:10, skipSelfT:0.9995 }
]
```
- `combine:'max'`:C 对 B 已得分类取两系综 max(strict > 保 B 值,确定性)。
- 最终排序键:`(-tanimoto, smiles, idx)`——**平手打破是契约**,稳定排序。
- 4-way 记账:A 输出 products=evaluated(A 的 filtered=枚举后未过 sanity 者);B/C/D 各记。

## 10. 兼容

- `window.__analogPool` 门面保留,stats() 字段超集不变(m6 契约)。
- CreateAnalogRun 命令与 facts['analog:<runId>'] 结构不变;runId = taskId。
- 引擎(WASM/Rust)零改动;worker 协议零改动。
