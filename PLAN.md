# Plan: Top-5 统一任务执行器 + 计算漏斗(已完成)

设计:**docs/task-executor-design.md**(冻结);契约:**docs/task-contract.md**(信封冻结)。

## 已完成(本轮,阶段 1+2)
- [x] docs/task-contract.md:信封/状态机(含 interrupted/stale)/seed 域分离/指纹规范化/canonical 序列化(IEEE754 位串,-0≠0)/4-way 记账/落账语义/analogScan input v1
- [x] app/tasks/ 五模块:contract(纯)/scheduler(全局令牌+预留交互槽+取消注册,acquire 即 release 闭包)/pools(模块分池+令牌+任务数回收+预热备用+迟到丢弃;__analogPool 门面契约保持)/funnel(barrier+确定性平手+预算截断+TASK_INTERRUPTED)/store(旁表 upsert+Done 跨 store 单事务+对账+session 心跳)
- [x] db.js v3(tasks store)+project.applyWithFacts extraOps+compat.getStorage/extraOps 转发
- [x] analogScan 迁移:A 阶段审计(filtered/notEvaluated 入档)→typed spec(提交即入旁表)→调度器+漏斗(B/C/D)→Done 单事务(facts 含 stages+accounts+buildId)→取消=cancelled/interrupted 无部分入档
- [x] 任务托盘最小版(计数+展开:kind/阶段/耗时/取消+终态最近 3)
- [x] boot 无条件对账(abandoned/ghost-complete)
- [x] Node +21(146);m12_tasks 14/14;run.sh 25 套件全绿
- [x] **红线验收:同 spec 池 1 vs 池 4 逐位相等(m12)**
- 修复记录:①scheduler pump 授予不计数→死锁;②pools drain acquire/shift 竞争→令牌泄漏(毒物混战测试抓出);③acquire resolve 数字被当 release;④取消映射以 th.cancelled() 为准非错误码;⑤reconcile 触发点从惰性 initPlatformStore 体内移到无条件 boot

## 阶段 3(已完成,本轮)
- [x] optimize3D → kind optimize(§9.2 冻结):inputRevision 双哈希守卫,stale=complete→stale 两步,rerun 哈希守卫
- [x] shape 搜索 → kind shapeSearch(§9.3 冻结,设计稿 shapeAlign 改名注明):编排壳迁移,farm 令牌按调度器总量封顶
- [x] rerunTask 快照重跑(shapeSearch 回填重跑/optimize 哈希不符即拒)
- [x] Node 151;m13_tasks2 8/8;run.sh 26 套件全绿

## 明确不做(设计 §11)
DSL/断点续跑/worker 内抢占/阶段让路/水位降载/每阶段 fact/桌面托盘/主线程同步任务/工作台批量迁移。
