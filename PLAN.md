# Plan: DIC 变换增量化 — 免分配热路径 + 无效功剔除 → 决策门复审 → v1.1.2/v1.2.0

## 背景

v1.1.1 后 DIC(opt-in)迭代数 0.44–0.66× 笛卡尔,但墙钟仅 aspirin GFN-FF
打平,其余 1.5–1.8× 慢:每迭代变换开销 ~0.5–1 ms,而 MMFF E+G 本体仅
~76–200 μs。开销构成(代码审查):

1. `f_and_g` 每迭代重算 `q(x)` 做漂移检查,但 DRIFT=1e9 下**永不触发**
   ——纯无效功(O(n_prim) 三角函数);
2. `project_out_tr` 每次调用 12 次 Vec 分配(6 基向量 + Gram-Schmidt);
3. `back_transform` 每次调用 2 个 Vec 分配 × 6 轮 Baker 迭代,线搜索每
   迭代调 2–3 次;
4. map/gx/gq 各 1–2 次分配。

## 任务

1. **微基准分解**:临时例程实测 back_transform / q / project_out_tr /
   grad_q / map 在 aspirin/ibuprofen 尺度的单次耗时,确认开销排序。
2. **无效功剔除**:DRIFT_REBUILD 为"实际上永不"时跳过 q_true 漂移
   检查(编译期常量分支;重建路径保留,DIC_DEBUG 不受影响)。
3. **回变换降本**:Baker 轮数 6→3(残差收缩测量为准,roundtrip 精度
   目标 ≤1e-4 Å/rad——现有 roundtrip 测试阈值同步);分配削减
   (scratch 复用或就地写入)。
4. **TR 投影降本**:改写为顺序修正 Gram-Schmidt——两个 scratch 缓冲
   交替,消 12 次分配;数学不变(正交投影顺序实现等价)。
5. **决策门复审**:全矩阵(3 fixtures × 2 引擎)重测迭代数/墙钟;
   门禁:GFN-FF 任一 fixture 墙钟 ≥1.2× 改善且其余不劣于 1.1×,
   MMFF 不劣于 1.5× → **默认切 internal**(v1.2.0,行为变化注
   CODE_STATUS/README);未达标 → 维持笛卡尔默认(v1.1.2),数据记录。
6. **发布**:版本/README/CODE_STATUS;wasm 构建 + node 冒烟;
   commit + tag(推送另请示)。

## 验收(实施后实测记录)

- `cargo test` 273/273、clippy 0、fmt 干净、benchmark_mmff 输出与 v1.1.1 一致(phosphirane/cyclobutene 既有离群点)
- 微基准(ibuprofen n=33):back_transform 138→~50 μs(B_q 扁平化行主流向 + 3→2 轮 + 预除 Λ)、q 漂移检查编译期剔除(原 24 μs/迭代纯无效功)、试探点→接受点映射缓存(每迭代省一次回变换);TR 投影实测仅 3–6 μs,**按测量数据放弃改写**(计划偏差:非瓶颈,不值得风险)
- 决策矩阵(同机同热状态两次运行稳定):
  - DIC 迭代数 0.43–0.87× 笛卡尔;墙钟:aspirin GFN-FF ≈持平(1.0×),ethanol/aspirin MMFF 1.4–1.6× 慢,ibuprofen 两引擎 1.2–1.25× 慢
  - 决策门(GFN-FF ≥1.2× 胜出且其余 ≤1.1× 劣化、MMFF ≤1.5×)**未过** → 默认维持笛卡尔,internal 维持 opt-in,版本 v1.1.2(纯性能优化,无行为变化)
  - 结论:DIC 现为“迭代数减半、墙钟仅贵 0–60%”的可用 opt-in(柔性大分子/FF 成本占比高时接近或优于笛卡尔);进一步反超需 FF 侧无关的结构性削减(增量基维护或缓存 b_matrix),另立项
- wasm 重建(node 冒烟 1.1.2 双路径);API 零变化
