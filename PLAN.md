# Plan: MMFF 键合项 E+G 融合(几何量单算)→ 目标 v1.2.9

## 背景

v1.2.8 后 E+G @ibuprofen 14.9μs。驱动对每个键合项分别调用
X_energy + X_gradient——r/θ/φ 等几何量算两遍(acos ~30-60ns ×
~100 次/求值是主要浪费)。E-only 路径 4.9μs 说明融合上限可观。

## 任务

1. **逐项融合**(bond/angle/stretch_bend/torsion/oop):新增
   X_energy_and_gradient 计算几何中间量一次,能量与梯度表达式
   逐位复刻既有两函数;驱动改调融合版。E-only 路径不动。
2. **逐项验证**:benchmark 230/230 逐字节(能量算式不变)、
   既有 FD 一致性测试、全量测试
3. **交错 A/B**(E+G μs + opt ms,vs v1.2.8)
4. **发布**:1.2.8→1.2.9;CODE_STATUS/PLAN;commit+tag;冒烟
   必须命中引擎输出行

## 验收(实施后实测记录)

- `cargo test` 280/280 全绿;benchmark 230/230 与基线逐字节一致
  (五项融合的能量/梯度表达式逐位复刻);clippy 0;fmt;wasm(node
  冒烟 1.2.9,引擎输出行确认);API 零变化
- **实施**:bond / angle / stretch_bend / torsion / oop 全部融合
  为 X_energy_and_gradient(几何中间量 r/θ/φ/cos 单算;E-only
  路径不动)。angle 融合版修掉了初版仍算两次 acos 的残余。
- **严格交错 A/B(3 轮,load ~12–16)**:
  - E+G:aspirin 10.6–11.6→8.7–9.0 μs(**~1.25×**);ibuprofen
    17.9–20.6→14.3–15.6 μs(**~1.25×**)
  - opt:aspirin 3.4→2.9 ms(1.19×);ibuprofen 9.3–9.7→8.1–8.4
    ms(1.15×)
- 累计(v1.2.4→v1.2.9):E+G @ibuprofen ~24×;opt 原生 15.3→8.2 ms
