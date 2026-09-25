# Plan: 单点 E+G 求值成本探索 → 目标 v1.2.5

## 背景

与 RDKit 基准:opt1 差距 2–11×,最大嫌疑是每次迭代的 MMFF94s
能量+梯度求值成本。需要先分离:(a) 算法/实现成本(原生 Rust 单点
E+G 多少 μs)、(b) wasm 开销、(c) 求值次数(线搜索)。

## 任务

1. **测量基线**:`examples/diag_eg.rs` — 解析一次、建 FF 一次、
   循环 E+G 1 万次(ethanol/aspirin/ibuprofen × MMFF94s/GFN-FF),
   报 μs/次。RDKit 参照:python `Minimize(maxIts=1)`×1000 摊销。
2. **剖析热点**:`sample` 采样原生 E+G 循环,定位前十行。
3. **候选优化**(按剖析结果,逐项做,逐项交错 A/B):
   - 每次调用分配(Vec::new 在热路径)
   - powf/exp 调用(MMFF vdW 7 次幂)
   - 非键双重循环(无截断?邻接排除逻辑)
4. **门禁**:benchmark_mmff 230/230(打印 6 位小数不得变化);
   cargo test 275;clippy 0;fmt
5. **发布**:1.2.4→1.2.5(若有可测收益);CODE_STATUS/PLAN 更新

## 验收(实施后实测记录)

- `cargo test` 277/277 全绿(+2 新 FD 一致性测试);benchmark 230/230
  与 v1.1 基线逐字节一致(2 个既有离群不变——能量函数未动);
  clippy 0;fmt;wasm(node 冒烟 1.2.5);API 零变化
- **单点 E+G(原生交错 A/B,3 轮,load ~20)**:
  - ibuprofen:266–341 → 53–67 μs(**~5×**)
  - aspirin:204 → 41 μs(**~5×**,负载校准后)
  - ethanol:8.1 → 3.1 μs(2.6×)
- **端到端优化(原生交错 A/B)**:ibuprofen 154–199 → 31–35 ms
  (**~5×**);与 RDKit 原生 C++ 差距 11× → ~2×
- **wasm opt1**:ibuprofen 183 → 64 ms(差距 11× → 3.7×);
  aspirin 56 → 16.5 ms(9× → 3.0×);ethanol 1.4 → 0.83 ms
  (2× → 1.15×,接近持平)
- 三项修改:vdW 解析梯度(原为前向 FD:每对 3 次全坐标克隆 + 6 次
  能量求值);torsion 解析梯度(原为 12 次 FD 求值 + 克隆/项);
  oop 分配修复(12 次克隆→1 次/项,保留中心差分——项数少不值得
  解析化)。GFN-FF 不受影响(独立梯度代码)。
