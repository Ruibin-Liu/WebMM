# Plan: 目标侧逐特征高亮进 3D(纯 app 端,零引擎改动)

## 背景

药效团层已有逐特征匹配(Pharm k/n 列),但点行回载对齐构象后 3D 特征球
无差别渲染——看不到"哪些特征对上了"。本增量把匹配状态渲染进 3D 视图。

## 设计

1. **pharmMatch 增目标侧分数**:对称计算 tBest[j] = 每个目标位点的最佳
   同型查询匹配(按目标位点自重叠归一);返回值增 `tBest`
2. **行点击携带匹配数据 + 自动 Embed**:loadSearchRow 增 pharm 参;
   process() 入口清 stale 状态,loadSearchRow 在 process 后回填
   `pendingFeatureMatch`;对齐构象(自带 3D 坐标)自动 embed3D() 并自动
   勾上 Features 开关(可手动关)
3. **特征球按匹配态渲染**:addFeatureSpheres 存在 pendingFeatureMatch 时
   — tBest[j] ≥ 0.5:类型色 alpha 0.55 正常半径(命中)
   — 否则:类型色 alpha 0.12 半径 0.35(暗)
   3D 状态行报告 "target features x/y matched";索引对应关系依赖
   applyTransformToSdf 保全原子序(已成立)与 colorSites 确定性序
4. **清理**:新分子 process 自动失效;离开检索上下文不残留

## 验收(Playwright + describe_image)

1. aspirin 查询 → shape 检索 → 点 salicylic 行:自动 Embed、特征球多数
   亮(6/8);点 glucose 行:明显亮暗混排(4/8)——目检
2. 程序断言:亮/暗球计数与 pharmMatch.tBest 一致(测试钩子读回)
3. 手动 Embed 普通分子(无 pending 数据)→ 特征球回归统一亮
4. 390px、零 page error、三模式回归;门禁 sanity(295 测试,纯 JS)
