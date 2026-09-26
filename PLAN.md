# Plan: 嵌入 E+G 内核逐项计时 → 定刀并实施 → v1.2.7

## 背景

工作守恒确立后,嵌入的唯一杠杆是**单次求值成本**。主 3D 300 迭代
= 47% 嵌入时间,其中 lr 对循环理论只需 ~3-6ms(16k flops × 300 ×
2)但实测 16.8ms——~60% 在其它项(K12/chiral/planarity/h_bond/
torsion_pref 的 FD dihedral)。先逐项计时定刀,再实施。

## 任务

1. **逐项计时**(ETKDG_ITERS 门控,etkdg_gradient_with_pairs /
   etkdg_energy_with_pairs 内):lr_pairs / K12 bonds / K13 angles /
   chiral / planarity / h_bond / torsion_pref 各项累计 μs/次与占比
2. **按数据实施 1–2 项**(候选:torsion_pref 的 dihedral FD 12×
   atan2/项 → 解析梯度(v1.2.5 同型已验证);h_bond 邻接优化;
   其它)
3. **门禁**:cargo test、ensemble 6/6、benchmark 230/230、
   clippy 0、fmt;原生 + wasm 交错 A/B
4. **发布**:1.2.6→1.2.7;CODE_STATUS/PLAN;commit+tag;冒烟
   必须命中引擎输出行

## 验收(实施后实测记录)

- `cargo test` 279/279 全绿(+2 解析-FD 一致性测试);ensemble 6/6;
  benchmark 230/230 与基线逐字节一致;clippy 0;fmt;wasm(node 冒烟
  1.2.7,引擎输出行确认);API 零变化
- **逐项计时地图(修复前,ibuprofen 单次嵌入)**:
  - H-bond 梯度 FD:51.1%(分子级中心差分:99 次全坐标克隆 +
    198 次 h_bond_energy 拓扑扫描/梯度调用)
  - planarity 梯度 FD:27.2%(ring/exocyclic 二面角 FD 每次 12 个
    atan2 + impropers/sp1 中心差分)
  - torsion_pref 梯度 FD:11.7%(同型二面角 FD)
  - lr 对循环(此前以为的瓶颈)合计仅 ~5%
- **实施**(能量函数逐位未动,仅梯度路径):
  1. H-bond:三元组预计算(几何无关拓扑)+ 解析梯度 → 51%→0.1%
  2. dihedral_gradient_contrib 解析化(cos φ 雅可比 + −dedphi/sinφ,
     sinφ<1e-12 跳过——Fourier 能量在该处 dedphi 同阶为零)→
     torsion_pref 11.7%→~2%
  3. impropers 解析化(χ=asin|q| 链式,q=(v1·N)/|N| 未除 |v1|,
     饱和区梯度零)与 sp1 线性解析化 → planarity 27.2%→~4%
  4. **修出两处真 bug**:角链式导数分母误用 |u||v|(正确为 |u|²)
     ——H-bond 与 linear 初版均有;FD 一致性测试抓出
- **严格交错 A/B(5 轮,load ~25)**:
  - 原生:aspirin 中位 24.4→8.8 ms(**2.8×**);ibuprofen
    62.3→17.8 ms(**3.5×**)
  - wasm(5×2 轮中位):aspirin 34.9→10.9(**3.2×**);ibuprofen
    49.7→17.6(**2.8×**)
  - 叠加 v1.2.6 snap 修复,对 v1.2.5:嵌入累计 ~4×
- 逐项计时改为 Option 门控(clock() 关闭时热循环只付一次分支)
