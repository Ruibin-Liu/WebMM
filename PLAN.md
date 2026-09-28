# Plan: 联合 shape+color 刚体对齐 —— 颜色进入优化目标(v1.5.0)

## 背景

颜色 v1 的姿态仍由 shape 单独决定,颜色只在对齐后打分——当氢键特征
方向与体积最优姿态冲突时会留分。ROCS 实践是联合优化。本增量把颜色
重叠加入引擎优化目标,并在引擎侧返回 Color T / Combo(JS 事后打分路径
退役)。

## 设计

1. **引擎**(`src/shape`):
   - `ColorSite { c, alpha, type_id }`(type_id 稳定整数枚举 donor/
     acceptor/pos/neg/hydrophobe/ring)
   - `color_overlap_grad(qs, ts, w, t)`:同型位点对高斯重叠 + 解析
     6-DOF 梯度(与 pairwise 形状代理同数学形式,K/β 公式复用)
   - `ShapeObj` 增加可选双份位点表:f_and_g = −(O_shape + w·O_color),
     w=1.0(v1 文档化;颜色位点稀疏,量级天然小于形状)
   - `align_colored(query_atoms, target_atoms, q_sites, t_sites, opts)`:
     与 align 同骨架(帧起点/随机/polish),重打分 = shape T(全量)
     + color T(位点对,同型)@最优姿态;返回 {tanimoto, color_tanimoto,
     combo, transform, ...}
   - 单测:颜色梯度 FD 奇偶(种子锁定);自对齐含色 = combo 2.0;
     颜色恒等式(自色 1、对称)
2. **WASM(additive)**:`shape_align_color_wasm(query_sdf, target_sdf,
   query_sites_json, target_sites_json, opts_json) -> JSON`;版本
   1.4.0 → **1.5.0**;d.ts 同步;`shape_align_wasm` 签名不动
   (sites JSON:[{i, t}] 原子索引+类型名;ring 位点 {atoms:[..], t})
3. **页面**:Shape 检索改调新导出(sites 已在感知层就绪);移除 JS
   事后 colorTanimoto 路径;列/排序不变(Combo)
4. **parity**:引擎 color T @固定姿态 vs JS 公式 @同姿态(同常数同公式,
   求和序不同 → 容差 1e-9);端到端:aspirin 自匹配 Combo 200;salicylic
   居次席不回归;**联合 ≥ 事后**统计(top-10 combo 平均不降——姿态为
   颜色让步的净效应应为非负,如个别下降如实记录)

## 验收

1. Rust 单测新增全绿(289 → ~293);clippy 1.98 = CI;fmt
2. node 冒烟:1.5.0;benzene-H 自对齐 shape 1.0;aspirin 自对齐
   combo 2.0
3. Playwright:Shape 检索全流程(引擎侧颜色)、top-k 化学直觉排序
   保持、390px、零 page error、sim/sub 回归
4. 既有回归:289 测试、LBDD parity 51/51 不受影响(纯增量)

## 边界与不做

- 不做:类型间权重表/可调 w、颜色项的全量包含-排斥、多构象联合
