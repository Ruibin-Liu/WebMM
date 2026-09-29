# Plan: 药效团层 —— 逐特征匹配(Pharm k/n)+ 3D 特征可视化(纯 app 端)

## 背景

颜色层给的是整体 Color T;药效团实践关心的是**逐特征命中**:"这个供体
有没有对上受体、苯环有没有贴上疏水区"。ROCS 的 feature map 即此。本
增量在 Shape (3D) 检索上叠加逐特征匹配报告与过滤,并在 3D 视图可视化
特征球。零引擎改动(对齐变换与位点感知均已就绪)。

## 设计

### T1 逐特征药效团匹配(Shape 结果列 + 过滤)

- 对齐后(引擎 transform 施于目标位点,sitePos 现成):对每个查询位点
  i,match_i = max_j(同型目标位点)O_ij / O_ii(高斯对重叠归一,封顶 1);
  **命中** = match_i ≥ 0.5
- 结果列「Pharm」= "k/n"(命中数/查询特征数);自匹配 = n/n
- 搜索栏过滤器:「Pharmacophore」下拉(any / ≥60% / ≥80% / all,默认
  any)作用于命中比例 k/n,先于阈值渲染;状态行报告过滤前后计数
- 排序仍按 Combo(Pharm 为过滤维度,不抢排序——文档化)

### T2 3D 特征可视化(单分子视图)

- 3D 面板增「Features」开关:当前分子的 colorSites 以小球渲染进
  3Dmol 视图(donor 红 / acceptor 蓝 / pos 紫 / neg 橙 / hydrophobe
  灰 / ring 黄绿;ring 画在质心,半径稍大;球半径 ~0.5 Å)
- 开关随视图重置/新分子加载保持状态; legend 用 title 提示(色-型
  对照写在开关 title)

### T3 验证

1. 自匹配:aspirin 查询自身 Pharm = n/n(全部命中)
2. 化学直觉:aspirin 查询下 salicylic 高命中、glucose 低命中;过滤
   「≥60%」后 glucose 出局而 salicylic 保留
3. 可视化:截图目检特征球(颜色-型对应、ring 质心位置)
4. 回归:sim/sub/shape 列不变(Pharm 列仅 shape 模式);390px;零
   page error;门禁 sanity(293 测试,纯 JS)

## 边界

- 单构象语义(与 shape 层一致,文档化);不做距离容差手工查询编辑器
  (后续);不做目标侧逐特征高亮进 3D(后续)
