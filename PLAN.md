# Plan: WASM 体积缩减可行性调查(诚实关闭,零代码变更)

## 背景

用户问:是否有必要尝试减小 wasm 体积(不损失功能与性能)。动机线索:
引擎加载慢于 RDKit 的观察。先量化、再实验、后结论。

## 现状

- pkg/webmm_bg.wasm = 1138 KB raw / **433 KB gzip9**(Pages 实际传输口径);
  webmm.js 60 KB / 8 KB gz。
- app/vendor/RDKit_minimal.wasm = 7161 KB / 2324 KB gz(第三方,Workbench
  加载两者,引擎只占传输的 16%)。
- Cargo [profile.release] 已是 opt-level=3 + lto + codegen-units=1
  (速度优先,体积零调优);机器此前无 binaryen,wasm-pack 一直静默
  跳过后端优化,当前 wasm 为纯 rustc LTO 产物。
- 内嵌 JSON(gfnff 109KB + mmff 10KB)非大头;数据段约 329 KB,
  其余为代码(698 个 LTO 后函数)。

## 实验(node 基准:30 次优化累计计时 × 3 轮取最优;奇偶 caffeine/
ethanol/butane 对冻结参考;实验后 Cargo.toml 与 pkg 均已逐字节还原)

| 变体 | raw | gz | aspirin | ibuprofen | 结论 |
|---|---|---|---|---|---|
| base(现役) | 1138 | 433 | 4.3ms | 12.5ms | — |
| wasm-opt -O4 | 1137 | 434 | 4.0 | 12.1 | 体积无收益 |
| wasm-opt -Oz+strip | 1131 | 433 | 3.8 | 12.6 | 体积无收益 |
| panic=abort+opt3 | 1135 | 433 | — | — | 无收益(−3KB) |
| panic=abort+opt-s/z | 1032 | 404 | 4.8 | **17.1** | −9% 体积换 −12~37% 速度,**否决** |

## 结论(诚实关闭)

1. **不建议做**:现役 433 KB gz 已接近该代码库的自然体积——所有
   无损杠杆(binaryen 后端优化、panic=abort)收益 <1%;唯一显著的
   opt-level=s/z 用 12–37% 优化性能换 9% 体积(29 KB gz),与仓库
   逐版本积累的性能基线(v1.2.4→v1.3.1 的 24×/5× 等)直接冲突。
2. 加载体验的真实瓶颈不在此:Workbench 侧 RDKit 2.3 MB gz 是引擎的
   5.4 倍(第三方 vendored);Demo/Playground 侧 433 KB gz 一次缓存,
   且加载态 UX 已处理(版本占位/引擎失败兜底)。
3. 若未来确有强需求,候选方向(均有代价,需单独立项):serde_json
   换更轻解析(数据段/派生代码占比需先做 twiggy 剖析);功能特性
   门控裁剪 demo 构建(违反"WASM 导出是公共契约"的稳定性承诺);
   brotli(不在 Pages 控制范围)。
4. 附带发现:harness 中复用 OptimizationOptions 对象会在第二次调用
   触发 "null pointer passed to rust"(wasm-bindgen 按值传参=move,
   非引擎 bug,生产页面每次新建 options 不受影响)。

## 验收

- 零代码变更;Cargo.toml 与 pkg/ 实验后逐字节还原(git diff 空、
  cmp 一致);cargo test 281/281;demo 页加载还原 pkg 正常(v1.3.1,
  零 page error)。
