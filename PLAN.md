# Plan: 构象管线原生化 — 批量 embed→attachH→optimize 单次 WASM 调用闭环 → v1.2.0

## 背景

现行多构象管线(ETKDG 嵌入 + 逐构象优化)算法全在 WASM,但编排在前端
JS:每个构象跨 WASM 边界两次 SDF 文本往返(JS 拼 SDF → Rust 解析 →
优化 → JS 再拼)+ 逐构象重复解析/建 FF + N 次边界调用。worker 农场
并行模型保留,消字符串搬运。

## 任务

1. **批量原生 API**(`src/lib.rs`):新 WASM 导出
   `generate_optimized_conformers_wasm(sdf_heavy, n, seed_base, engine,
   max_iterations) -> OptimizedConformers`:
   - 重原子 mol 解析一次;逐构象:ETKDG(seed_base+i,复用现行种子语义)
     → `add_hydrogens`(几何感知 H 放置)→ 优化(与 optimize_dispatch
     同语义:MMFF94s/MMFF94/GFN-FF 分支,GFN-FF 的 conv 放宽;
     **MMFF FF 构建一次全批复用**(拓扑与坐标无关),GFN-FF 逐构象
     重建(拓扑来自几何));
   - 返回:flat coordinates、energies、converged、iterations(u32)、
     seeds、n_atoms/n_heavy/n_confs/success/error,以及全氢模板
     molblock(conformer 图相同,仅坐标异——供 JS 惰性构 SDF);
   - 沿用 n∈1..=500 与分块语义(与 generate_conformers_wasm 一致);
     附加导出 additive,d.ts 随构建更新。
2. **等价性测试**(Rust):批 API 与逐步组合(现行公开路径:批量嵌入
   → molblock_with_h → 优化 dispatch)对同一 seed 集**能量逐位一致**、
   n_confs/形状一致;MMFF FF 复用路径与逐构象建 FF 结果逐位一致。
3. **worker 改造**(`app/conf.worker.js`):`run` 分支改调新 API
   (分块 ≤500,块内逐构象 postMessage 保持流式);消息协议:
   `meta` 增带全氢模板,`conf` 改传 coords(主线程用自己的
   buildSdfFromCoords 构 SDF);主线程 `runConformers`/`conf` 处理适配
   (批量 3D 模式走同一 worker 自动受益);单结构 optimize 消息与
   ΔE 重优化路径不动(小 N,另议)。
4. **性能 A/B**(node,pkg 同模块):ibuprofen(33 原子)50 构象
   MMFF94s 与 GFN-FF 各一轮,旧路径(generate_conformers_wasm +
   attach_hydrogens_3d_wasm + optimize_from_sdf 逐构象)vs 新 API,
   记录墙钟;如无净胜(边界节省 < FF/优化本体)如实记录并保留 API
   (结构收益:消双份 buildSdfFromCoords 维护面)。
5. **发布**:版本 1.1.2→1.2.0(additive WASM 导出);README 构象管线
   描述更新;CODE_STATUS;wasm 构建 + node 冒烟;commit + tag(推送
   另请示)。CDP:M2 浏览器套件本会话无 harness,标注需下次站点
   会话复跑。

## 验收

- `cargo test` 全绿(+等价性/形状/复用一致性);clippy 0;fmt 干净
- `benchmark_mmff.py --no-speed` 输出与 v1.1.2 一致
- node 冒烟:新 API 50 构象成功、能量与逐步路径一致、模板 SDF 可被
  buildSdfFromCoords 往返
- A/B 数据入 CODE_STATUS;API:WASM 导出仅新增(additive)
