# Plan: 多构象 Shape 检索性能优化——worker 并行 + 按查询 screen 缓存

## 背景

多构象 shape 检索(55 库@10 构象)冷跑 14.6s:生成 ~4.1s + screen
550 对 9.6s + phase 2 ~0.9s(缓存热后 8.0s)。三段全部按条目独立、
embarrassingly parallel;同查询重跑(调阈值/pharm 过滤——最常见迭代
流)仍重复 screen 7.1s。

## 实施(app 端,零引擎改动)

1. **app/shape.worker.js**(新 worker,仿 conf.worker.js:init 消息
   携 wasm URL、自含 sdfLines/buildSdfFromCoords):
   - 'prep' {idx, mb, n, seed, qsdf, qSitesJson, inject}:批量生成
     N 构象(seed 42+idx、MMFF94s iter 100——与现行参数逐字相同)+
     逐构象 buildSdf + inject 时构象 0=qsdf + 逐构象
     shape_align_wasm(screen) → 返回 {confSdfs, proxies}
   - 'full' {idx, ci, esdf, tSitesJson, useColor}:缓存 qsdf/
     qSitesJson,shape_align(_color)_wasm(random_starts 8) → 返回
     {tanimoto, color_tanimoto, transform}
   - 与主线程顺序版参数/选项串逐字一致 → 结果确定性不变
2. **runShapeSearchEnsemble 重写为农场编排**:主线程只做 RDKit 侧
   (重 molblock、colorSites[phase2 胜者候选]、pharmMatch、
   applyTransformToSdf)与排序/渲染;W=min(hardwareConcurrency,8)
   worker 动态领活(逐条目派发);运行令牌防陈旧结果渲染;结束
   terminate。状态行格式不变(m5 断言依赖)
3. **按查询 screen 缓存**:条目上缓存 {screenQ: 规范查询 SMILES,
   screenProxies};命中则跳过生成+screen(复用 e.confs);库对象
   重建自然失效;inject 状态由 screenQ 键涵盖
4. 回退不做(worker 全 app 已依赖,file:// 下 wasm 本就不可用——
   文档化惯例);单条目失败跳过(现行语义)

## 验收

- 结果不变性:同库同查询 worker 版与顺序版 top-10(名称+ShapeT+
  Conf)全等(实施中用顺序版快照对照一次后移除);m5 32/32 全绿
- 性能:55 库@10 构象冷跑 ≤5s(原 14.6);同查询重跑 screen 0s
  (仅 phase 2);@5 与 N=1 不回归
- m1/m2/m4 抽查;390px;零 page error;node --check(含新 worker)

## 验收结果(实施后)

- 性能(55 库@10 构象,8 worker):冷跑 14.6s→**5.0s**(3×);同查询
  重跑 11.9s→**1.4s**(screen 0s,仅 phase 2);换查询 **2.4s**(构象
  复用只重筛,不再重新生成);N=1 旧路径不变
- 确定性:同查询多次运行 top-3 全等(100.0/71.9/65.4);worker 参数
  串与顺序版逐字一致
- 过程修掉三个真 bug:①worker 'prep' 未携带 wasmUrl(瞬间全败成
  空结果);②phase-2 resolve 挂在 worker 消息而非 waiter 注册表
  (Promise.all 永挂);③**async onmessage 不排队**——'full' 在
  'initq' 的 wasm import 完成前并发执行致 wasm=null 全败(消息队列
  链串行化修复;全缓存路径另以 'initq' 会话引导 + 失败哨兵兜底)
- 缓存语义:screen 按(条目,N,规范查询)缓存;换查询时恢复构象 0
  的纯生成副本(genSdf 保存)并对新查询成员重新注入;worker 数按
  库规模(全缓存重跑仍并行 phase 2);运行令牌丢弃陈旧结果
- m5 **32/32**、m1 10/10、m2 11/11、m4 32/32;390px shape 模式无
  溢出;零 page error;node --check 全部(含新 worker)
