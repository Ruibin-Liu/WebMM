# Plan: 构象管线第三轮 — ETKDG 嵌入内部提速 + 集成首轮迭代协议 → v1.2.1

## 背景(实测)

ibuprofen 30 构象管线:嵌入 3.2s(106ms/构象,55%)+ MMFF 优化 2.5s
(43%)。剖析嵌入:~30% 在 minimize_etkdg 的逐步 Vec 分配(闭包内
x_new/g_new/s/y × 400 迭代 × ≤25 线搜索试探);~17% 在扭转 snap 循环
(每个 snap 做一次完整再最小化);~3% 逐构象重建 bounds 等(拓扑不变)。
优化侧:30/30 撞满 250 迭代上限、0 收敛退出——迭代数由上限决定,
首轮松弛降到 ~100 迭代是协议层免费收益(app 已有 ΔE 窗口完全收敛
重优化兜底)。

## 任务

1. **minimize_etkdg 缓冲复用**(`src/etkdg/mod.rs`):闭包逐步分配
   改为预分配 scratch(x_new/g_new/s/y/方向),数学与接受逻辑不变
   ——能量序列逐位一致由种子回归测试锁(既有 ETKDG 测试集 +
   ensemble_stats_vs_rdkit 6/6 门禁不回归)。
2. **扭转 snap 重排**:逐 snap 完整再最小化 → snap 序列在轻量局部
   松弛(或直接快照评估)后一次最终 minimize;实现方式:snap 循环
   里 minimize_etkdg 的迭代预算调小或跳过中间最小化,最终一次完整
   收尾。以 ensemble 统计门禁(同盆地/收敛数)裁决具体形态。
3. **bounds 等拓扑预计算共享**:generate_optimized_conformers_wasm
   的逐构象循环里,bounds/平滑/平面约束/手性/扭转偏好构建提出循环
   外(embed_impl 加预计算结构参数或拆分 API);单次调用路径行为
   不变。
4. **首轮迭代协议**:conf.worker.js 的构象路径 maxIter 250→100
   (ΔE 窗口重优化已有,首轮只需排序质量);单结构 Optimize 路径
   不动(仍自适应 40×原子数)。
5. **实测门禁**:嵌入微基准(ms/构象,目标 ≤50,争取 ~30)、
   端到端构象管线(ibuprofen 30 构象 node A/B,目标整体 ≥2×);
   ensemble_stats_vs_rdkit 6/6、ETKDG 种子回归、cargo test 全绿。
6. **发布**:1.2.0→1.2.1(纯性能+协议默认值,无 API 变化);
   CODE_STATUS/README;wasm 构建;commit + tag(推送另请示)。

## 验收(实施后实测记录)

- `cargo test` 275/275 全绿(ETKDG 能量序列逐位不变——pair 预计算
  保持相同顺序与累加;minimize_etkdg scratch 缓冲不改变数学);
  clippy 0、fmt 干净
- ensemble_stats_vs_rdkit 6/6 不回归;benchmark_mmff 230/230 一致
- **诚实性能结论(交错 A/B,同机同热)**:
  - **ETKDG 内部(pair 预计算 + scratch 缓冲):在药物分子尺度(33
    原子/528 对)无可测收益**(< 3%,在噪声内)——skip 测试太便宜,
    对列表太短;收益需 100+ 原子分子才可能显现。保留实现(正确、
    架构合理、逐位一致),但预期收益如实记录为零。
  - **首轮迭代协议 250→100:1.6× 优化阶段提速**(端到端 30 构象
    ibuprofen:3250→2010 ms)。30/30 全部由 max_iterations 决定
    退出,力阈值放宽无效果(7500→6762 迭代,10% 降幅)——协议
    层改 cap 是唯一杠杆。
  - **嵌入仍是最大单项**(~40% 管线时间,41-89 ms/构象,随分子
    大小超线性);内部 L-BFGS 的 400×4D + 300×3D 迭代是本质
    成本——pair 预计算没改变这一点。
- wasm 重建(node 冒烟 1.2.1);API 零变化
