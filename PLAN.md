# Plan: 药效团查询编辑器 —— 特征点选取 + 距离容差筛选(纯 app 端)

## 背景

药效团实践的核心形式:**特征类型 + 3D 位置 + 点对距离约束(带容差)**。
已有逐特征匹配是对齐驱动的(先对齐再比对);经典 pharmacophore query
是**纯几何约束求解**:目标分子(构象)中存在同型特征指派使所有点对距离
落在容差窗内。本增量在 Search 标签加查询编辑器面板。

## 设计(单构象语义,文档化)

1. **查询构建**:「Build from current molecule」取单分子视图 3D
   (sdf3d;无 3D 则提示先 Embed)的 colorSites;特征复选列表
   (type @atom idx + 坐标),最多选 6 个;容差输入 ± Å(默认 1.5);
   实时距离矩阵预览(所选特征间当前距离,只读)
2. **筛选算法**(页面 JS):逐库条目(ensureEntry3D 全档构象)——
   - 候选:同型目标特征(位点感知已就绪)
   - **回溯指派**:特征按候选数升序排列,逐个指派并即时检查与已指派
     特征的距离约束(|d_query − d_target| ≤ tol),剪枝;全指派成功 =
     Match
   - k ≤ 6、候选项少,组合爆炸可控
3. **结果表**:#/Name/SMILES/Match 徽章;行点击回载条目(无对齐 SDF);
   状态行 matched/total;**不做 CSV(后续)**
4. 语义边界(面板脚注 + CODE_STATUS):单构象(与 shape 层一致);
   纯几何(不做对齐);类型口径同颜色力场六类

## 验收(Playwright)

1. aspirin 查询选 4 特征(donor/双 acceptor/ring)tol 1.5:salicylic
   Match(COOH+环几何同族)、glucose 不 Match;自匹配 Match 且 0 违约
2. 极端容差 0.1:自匹配仍 Match(距离精确同源);乱序距离约束
   (tol 极小 + 类型不匹配组合)正确拒绝
3. 回归:三模式/RGD/骨架面板不受影响;390px;零 page error;
   门禁 sanity(295 测试,纯 JS)
