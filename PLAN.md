# Plan: ETKDG 嵌入求值次数削减 → 目标 v1.2.6

## 背景

SIMD 已诚实关闭;嵌入(55–60 ms/构象 @33 原子)的杠杆是减少
求值次数/阶段成本。此前剖析行级归因粗糙(~40% 在内联黑洞),
需要先拿到**每阶段的迭代数 × 能量/梯度求值次数 × 墙钟**的
精确地图,再定刀。

## 任务

1. **插桩**:ETKDG_ITERS 环境变量门控,在 embed_impl 各阶段
   (4d_first / 4d_collapse / flatten / 3d 主最小化 / 每个 snap
   重最小化 / H-only / trilaterate / 验收检查)打印
   [迭代数, 能量求值次数, 梯度求值次数, μs]。
   minimize_etkdg/lbfgs_minimize 返回或累计求值计数。
2. **定刀**(按实测地图,候选):
   - snap 重最小化预算/触发条件(占 17% 墙钟)
   - 线搜索试验次数分布(过多次 Armijo 拒绝?)
   - h_bond / torsion_pref / dihedral FD 等小项的每求值成本
   - 设置复用(lr_pairs/scratch 跨 minimize_etkdg 调用)
3. **实施 1–2 项**,逐项交错 A/B + 门禁(ensemble 6/6 裁决质量)
4. **发布**:1.2.5→1.2.6;CODE_STATUS/PLAN;commit+tag

## 验收(实施后实测记录)

- `cargo test` 277/277 全绿;ensemble 6/6 通过;clippy 0;fmt;
  wasm(node 冒烟 1.2.6);API 零变化
- **求值地图(新插桩工具,ETKDG_ITERS 门控,ibuprofen seed 42)**:
  主 3D 最小化 300/300 跑满(37ms,47%)+ 3 个 snap 各 ~95 迭代
  (34ms,43%)——90% 在 4 次 minimize_etkdg;4D 阶段仅 1.6ms;
  H-only 1 迭代;trilaterate 2μs。
- **实施**:snap 重最小化预算 300→25(局部弛豫职责,全场景收敛
  属主最小化与验收门禁);**主 3D 300→150 试验失败回退**——工作
  转嫁给 H-only(1→50 迭代,总时间反升),300 承重。
- **收益(插桩单机对比,79.6→64.1 ms/embed,~-20%;迭代总数
  900+→339)**。⚠️ 机器后被外部任务饱和(load 30→190),严格
  交错 A/B 未能完成——数字为临时值,待安静机器复核。
- 新增结构化插桩(阶段 × 迭代/能量/梯度计数 × μs)保留为
  ETKDG_ITERS 门控的常驻审计工具(后续轮次的地图生成器)。
