# Plan: Top-5 统一任务执行器 + 计算漏斗(设计已冻结,待用户批准实施)

设计文档:**docs/task-executor-design.md**(Fable-5.1 设计评审 + GPT-6-Astra 红队,收敛稿已冻结)。

## 实施计划(批准后执行)

### 阶段 1:契约 + 骨架
1. `docs/task-contract.md`:冻结信封(TaskSpec 字段/状态机含 interrupted/seed 域分离规则/参数指纹规范化算法/buildId=SHA-256(wasm)/四分类记账)。
2. `app/tasks/` 五模块骨架:contract.js(类型化构造+指纹)/ scheduler.js(两级优先级+全局令牌+预留槽+取消队列层)/ pools.js(模块分池+N 任务回收+预热备用+taskId+generation 迟到丢弃)/ funnel.js(barrier+确定性预算截断+检查点)/ store.js(IDB tasks 旁表 upsert+同库跨 store Done 事务+对账+sessionId 心跳)。
3. Node 测试:指纹规范化(-0/默认/数组序)、状态机迁移、对账(幽灵 complete/无主 abandoned)、幂等(重复 Done)。

### 阶段 2:analogScan 迁移(验收里程碑)
4. A/B/C/D 胶水 → typed pipeline 数据 + funnel.js runner;有效 stages 快照入 spec;四分类计数;同步 A 阶段输入上限+阻塞实测。
5. **并发等价测试:同 spec 并发 1 vs 4 逐位相等**(CDP);干净/用过/重建 worker 一致性(Node)。
6. 任务托盘最小版:全局活动入口(计数→任务/阶段/取消/结果去向),页面状态行保留。
7. CDP m12_tasks;CreateAnalogRun 走任务层(CreateAnalogRun 契约不变,facts 结构不变)。

### 阶段 3:optimize3D + shapeAlign
8. optimize3D:inputRevision=molblock 哈希+restraints 指纹;stale 路径(结构变更→"对旧结构重算/丢弃")。
9. shapeAlign:query+target 姿态哈希+color_scale 入 revision。
10. 各自 kind input schema 冻结进 task-contract.md。

### 不做(设计 §11):DSL/断点续跑/worker 内抢占/阶段让路/水位降载/每阶段 fact/桌面托盘/主线程同步任务/工作台批量迁移。

## 验收
- Node + 新增任务层套件全绿;CDP m12 + 既有 24 套件零回归;analogScan 并发等价逐位相等;
- CreateAnalogRun/facts 既有契约零破坏(m6_analog_history 全绿即证);
- clippy/fmt/cargo test 不受影响(纯 JS 层,引擎零改动)。
