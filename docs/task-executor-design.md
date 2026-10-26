# Top-5:统一任务执行器 + 计算漏斗 —— 收敛设计

状态:**设计冻结**(实施前评审稿)。评审:Fable-5.1(设计评审)+ GPT-6-Astra(红队),2026-10-27。
本文与 docs/task-contract.md(实施第一步冻结的契约信封)配合:本文讲**为什么与整体架构**,契约文讲**字段级规范**。

## 0. 问题

计算编排散落 9 处(工作台构象/批量、平台 analog 探索/optimize3D/shape 叠合、主线程嵌入/检索/SA),同一"3D 任务"概念每处重新发明;漏斗 A/B/C/D(枚举过滤→粗 3D→系综重评→flex 精化)硬编码在 ~600 行页面胶水;进度/取消/持久化口径不一;无统一成本/重跑/废弃语义;页面重载=进行中任务静默蒸发。

## 1. 分层架构(绞杀者,非重写)

```
app/tasks/                    ← 共享级(Fable D1:两页薄壳,平台先迁)
  contract.js   类型化 TaskSpec + 状态机 + 参数指纹(规范化序列化)
  scheduler.js  优先级队列 + 全局并发令牌 + 取消 + 对账
  pools.js      每 wasm 模块一池(worker 复用/回收/预热备用/迟到丢弃)
  funnel.js     通用漏斗 runner(barrier 选择/预算/硬顶/流式/检查点)
  store.js      IDB tasks 旁表(upsert;唯一 Done 前写入)
两页薄壳:platform.html / index.html 只做 UI 绑定
```

## 2. 任务契约(信封现在冻,kind input 随迁移冻)

**TaskSpec**:`{ taskId, kind, version, input, seed, provenance: { projectId, roundId?, molId?, buildId, policyRev? } }`

- **v1 kinds**:`conformers`(嵌入+优化)/ `optimize`(含 restraints)/ `shapeAlign`(含 color_scale)/ `analogScan`(整条漏斗管线)。
  Astra Q1:**分 kind 不合并**——optimize 与 conformers 输入前置/输出基数/失败含义/应用行为不同;共享内部执行代码与 worker 协议,不做 mode 判别联合(类型错误推到运行时不诚实)。
- **seed 派生**:`hash(roundSeed, kind, stageId|op, molId, confIdx)` —— Astra:域分离纳入稳定 stageId/op 与随机流版本,防不同操作意外共享随机流。
- **buildId** = SHA-256(wasm 字节)启动时算一次(subtle.crypto,缓存);RDKit 用其版本串。Astra:复现身份须含实际构建身份,不只是显示版本号。
- **参数指纹**:唯一规范化算法——键排序、浮点 IEEE754 位串、显式默认值/缺省/-0/非法值/数组顺序规则;**对展开默认值后的实际执行输入计算,不是对声明输入**(Astra)。
- **状态机**:`queued / running / complete / failed / cancelled / abandoned / stale / interrupted`。
  `interrupted`(Astra #1):时间/内存/取消导致的非完整执行,绝不冒充完整排名。

## 3. 调度器

- **两级优先级**:interactive / background;并发 = min(4, hardwareConcurrency) **全局令牌**(Fable D2,非每池各自),interactive 预留 1 槽(可占满)。
- **抢占只在队列层**;无 worker 内抢占、无阶段间让路(Fable D4:让路=丢弃重算或持久化中间态,两者都不值)。
- **漏斗是逻辑任务,不持令牌**:其全部 worker 作业统一走全局令牌(Astra #5:父任务占槽再申请子作业会死锁/绕上限)。
- **内存**:不做水位降载(Fable)。v1 保护 = 保守输入规模硬顶 + 活跃 worker 与预热总数上限 + 消息批次设界(Astra #6:byteLength 不覆盖堆外/临时副本/积压)。
- **worker 生命周期**:每池处理 N 个任务或内存超阈值→空闲时 terminate 重建;每池留 1 个预热备用。**硬取消 = terminate + 备用顶上**;协作取消(message 边界)仅为优化非保证(无 SAB 时常不可达);取消后迟到 reply 按 `taskId+generation` 丢弃。
- **同步路径豁免**:主线程嵌入/检索/SA v1 不纳入(伪任务=不可取消不可并发的形态)。但 analogScan 若依赖同步 A 阶段,**迁移验收必须含该路径**:输入上限 + 实测最大连续阻塞时间;UI 取消区分"已请求/已停止"(Astra #7)。

## 4. 确定性契约(红线,写进 task-contract.md)

1. 每阶段 **barrier**:全部候选完成后才做 top-n 选择;流式行只喂 UI,不参与选择。
2. **排名截断只允许确定性工作量预算**(候选数/构象数);超时/OOM/取消 → `interrupted`,不得产出"完整排名"(Astra #1)。
3. 稳定排序 + 确定性平手打破(分数→molId→confIdx)。
4. 每候选独立 seed;worker 无跨任务状态(WASM 全局 RNG/缓存/ETKDG 参数对象残留——worker 回收顺带兜底)。
5. **迁移验收 = 并发 1 vs 并发 4 逐位相等自动化测试**;另加"同候选在干净/用过/重建 worker 一致性"测试(Astra #8)。

## 5. 阶段记账(科学诚实)

每阶段显式四分类(Astra #2,呼应"无证据≠零分"红线):
`evaluated`(真算了)/ `filtered`(确定性规则筛掉)/ `failed`(算了但失败)/ `notEvaluated`(预算截断未轮到)。
"通过数+触顶"不够——可用性与质量不得混淆;排名域与失败策略入档。

## 6. 持久化

- **IDB tasks 旁表**(upsert):TaskSpec **提交时**写入(唯一 Done 前写入;在旁表不在 facts,不破坏"Done 前零写 facts"不变量)+ 状态迁移 + 耗时/指纹。进度全在内存。
- **facts(L3)**:Done 时整任务一条(现状 AnalogRun 语义),阶段摘要内嵌且**只含影响结果的字段**(cap/n/seedOffset/通过数/触顶/四分类计数);耗时/worker 编号留旁表(Fable D3:每阶段一条=写放大+下游先聚合才可用)。
- **写入顺序与对账**:Done = 同库跨 store 单事务(fact+旁表 complete;Astra #4:先写后写的窗口会遇重复 Done/取消竞争),taskId+attemptId+generation 幂等。重载时旁表 running 且无存活 session 认领 → abandoned;若对应 fact 已在 → complete(Fable D5 幽灵对账)。
- **多标签**:旁表记 sessionId,BroadcastChannel 心跳(内存),启动只废弃无主任务。
- **断点续跑否决**:中间态是构象坐标(大),持久化破坏不变量且要求引擎逐位一致;abandoned + 手动 rerun(从 spec 重算,成本=算力,换不变量干净)。
- **analogScan input 提交时保存展开默认值后的有效 stages**;旧任务 rerun 绝不静默套用当前默认漏斗(Astra Q2);提交前校验阶段衔接/支持 op/预算合法性——只写本项目所需校验器,**不发展成通用工作流框架**。

## 7. 结果应用校验

TaskSpec.input 含 **inputRevision per kind**(Astra #3:结构哈希不够——重排/坐标/约束/查询变化都该使旧结果失效):conformers/optimize→molblock 内容哈希+restraints 指纹;shapeAlign→query+target 姿态哈希+color_scale;analogScan→parentSmiles+fragmentSet rev。应用前校验目标当前 revision + 请求身份(taskId+generation);不匹配 → `stale`,UI 给"对旧结构重算/丢弃"。

## 8. UI(裁剪后)

- 页面内沿用现有状态行。
- **全局活动入口**(Astra Q3 裁剪,不做桌面式任务管理器):运行/排队计数 → 展开 {任务名, 阶段, 取消, 结果去向};失败原因/参数/指纹入详情;终态历史留在项目结果处(facts),托盘不做第二套结果管理器。
- rerun 显示**可读参数差异**,不显示指纹 diff(用户无解释价值)。

## 9. 数据搬运

坐标数组 Transferable;流式行 worker 侧按时间窗(100ms)合并再发,防 A 阶段把主线程打成每秒数百次更新。

## 10. 迁移顺序

1. `docs/task-contract.md` 冻结信封 + app/tasks/ 骨架 + 并发等价测试基建
2. **analogScan**(最大收益;验收=并发 1v4 逐位相等 + 同步 A 阶段阻塞实测 + 四分类记账)
3. `optimize3D`(restraints 接 inputRevision)
4. `shapeAlign`
5. 工作台批量(m4)与构象——后置,独立立项
每步:Node + CDP + run.sh 全量;契约破坏性变更须在 PLAN.md 显式声明并记 CODE_STATUS。

## 11. 明确不做(两模型一致 + 既往裁决)

DSL 文本编辑器(schema 先行)/ 断点续跑 / worker 内抢占 / 阶段间让路 / 内存水位降载 / 每阶段一条 fact / 跨标签任务同步(心跳仅认领)/ 服务端 / 通用工作流框架 / 主线程同步计算纳入 v1 / 桌面式任务托盘。

## 评审记录

- **Fable-5.1**(2026-10-27):D1 app/tasks/ 共享两页薄壳;D2 模块分池+全局令牌+回收+预热备用;D3 整任务一 fact+旁表;D4 队列层抢占+预留槽;D5 abandoned+对账+sessionId;D6 信封先冻;D7 确定性红线+硬取消=terminate+迟到丢弃+写三类+Transferable+批流。
- **GPT-6-Astra**(2026-10-27):#1 interrupted 状态与确定性预算截断;#2 四分类记账;#3 inputRevision;#4 同库跨 store 事务幂等;#5 漏斗不持令牌;#6 byteLength 不可靠→保守上限;#7 同步路径豁免须入验收;#8 指纹规范化+域分离+构建身份+脏 worker 测试;Q1 分 kind;Q2 有效 stages 快照;Q3 托盘裁剪。
