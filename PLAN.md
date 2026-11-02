# Plan: Metad 构象搜索设计基线(iMTD-GC 迭代路线)+ 能量景观探索方向记录

用户指示记录两项设计决策(本轮纯记录,不改引擎代码)。

## 交付
1. **Metad 构象搜索设计基线 = CREST iMTD-GC 迭代路线(排期实施时的设计约束)**
   - 组织方式:多轮**短程** metadynamics——轮换偏置强度与 CV 组合(DihedralCV/DistanceCV 的不同选法),每轮从轨迹收集构象;对低能构象做**遗传交叉**(genetic crossing)产生新种子;再进入下一轮,直至构象系综**收敛**。
   - 明确不做单次长程 metad 充当构象搜索基线。
   - 现有 well-tempered metad 引擎**原样保留**(src/metad/mod.rs:MetaDynamics/DihedralCV/DistanceCV,lib.rs MetaDOptions/run_metadynamics_from_sdf/MetaDLive);迭代-交叉调度加在其**外层**(轮次循环 + 构象收集/去重 + 交叉产种 + 收敛判定),不改内层 metad 语义。
2. **远期潜在应用方向(仅记录,不排期):浏览器端能量景观探索视图**
   - 基于现有 metad/优化(L-BFGS 局部极小化)/轨迹能力:
     - 极小值网络(minima network);
     - 断连图(disconnectivity graph);
     - 盆地统计:轨迹帧周期性淬灭(quench 到所属局部极小),给出盆地占据分布与盆地间跳转分布。

## 验收
两项指示原文语义完整落入本文件对应章节;引擎与前端零改动、零行为变化(纯文档轮)。

## 排雷存档(本轮)
- 无(纯记录轮,无代码变更)。
