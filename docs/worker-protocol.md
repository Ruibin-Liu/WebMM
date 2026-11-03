# Worker 消息协议(冻结于 v1.12.0)

本文件钉死三个 ES-module worker 的消息面。任何增删改都是破坏性变更,须按
docs/STABILITY.md 的纪律执行(golden 更新 + 版本升级)。守卫测试:
`tests/cdp/m22_api_surface.test.js`(精确相等,增删改皆红)。

## conf.worker.js(构象/优化/metad 家族)

| 消息方向 | type | 载荷 | 回帖 |
|---|---|---|---|
| 主线程→worker | `init` | `{ wasmUrl }` | `{ type:'ready', seedBase? }` |
| 主线程→worker | `optimize` | `{ sdf, engine, maxIter, restraints? }` | `{ type:'optimized', E, converged, iterations, message, engine, nAtoms, coords, charges, termsJson }` / `{ type:'error', message }` |
| 主线程→worker | `run` | `{ sdfHeavy, seedBase, count, engine, maxIter }` | 流式 `{ type:'meta', seedBase, nConfs, template }` + 每构象 `{ type:'conf', seed, E, coords, sdf }` + `{ type:'done' }` |
| 主线程→worker | `imtdRun` | `{ id, runIndex, sdf, specJson, optSets:[[setter,value],…] }` | `{ type:'imtdMinima', id, runIndex, minimaJson }` / `{ type:'imtdError', id, runIndex, error }` |
| 主线程→worker | `imtdNeb` | `{ id, index, sdf, jobJson, optSets }` | `{ type:'imtdNebDone', id, index, resultJson }` / `{ type:'imtdNebError', id, index, error }` |

注:`imtdRun`/`imtdNeb` 的 `specJson`/`jobJson`/`resultJson` 为引擎侧
bit-exact JSON(f64 走 17 位指数字符串);`optSets` 为 ImtdOptions setter
对,worker 保持 schema 无关。

## shape.worker.js(形状叠合 farm)

| type | 说明 |
|---|---|
| `initq` | 初始化(wasmUrl + 查询态) |
| `prep` | 系综准备(每 run 专属) |
| `screen` | 打分(流式回帖) |
| `full` | 完整叠合(联合优化,colorOpts) |

## analog.worker.js(SAR 类电扫描 farm)

| type | 说明 |
|---|---|
| `init` | 初始化(wasmUrl) |
| `score` | 单候选打分 |
| `heap` | top-n 堆维护 |

变更流程:先改本文件与 m22 golden,同轮升 minor 版本,再动实现。
