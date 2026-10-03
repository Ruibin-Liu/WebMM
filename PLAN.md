# Plan: Features 可视化三项改进(纯显示层,A+B+D;零算法/引擎改动)

## 背景

用户观察:Features 半透明球无方向感且不可读。方向性(C,投影位点)涉
Color T 语义已单独搁置;本轮做 A 图例、B 类型化视觉、D 类型清单。

## 实施(app/index.html 单文件)

1. **A+D:图例 = 类型清单**(带计数):选项行 Features 复选框后加
   `<span id="featLegend" class="hint">`,chips = 色点 + "type ×n"
   (仅当前分子存在的类型;六类降序按计数);addFeatureSpheres 渲染,
   入口/未勾选/无位点时清空;onFeatureSpheres 未勾选时清空
2. **B:类型化视觉**:FEATURE_COLORS.hydrophobe 灰 #94a3b8(贴碳隐形)
   → 琥珀 #a16207(与橙 neg、黄绿 ring 拉开);默认 alpha 按类型
   (donor/acceptor/pos/neg 0.45→0.55,ring 大球保持 0.45);两档亮度
   匹配逻辑(0.6/0.12)不动;复选框 title 文案同步("amber
   hydrophobe");pharm 查询编辑器色点同源自动跟随
3. 零引擎/算法/wasm 改动;colorSites、检索、打分全部不动

## 验收

- Playwright:图例 chips 与 colorSites 计数一致(aspirin:acceptor×4/
  hydrophobe×6/donor×1/neg×1/ring×1 序按计数);未勾选不显示;
  390px 不溢出;零 page error;hydrophobe 球颜色 #a16207
- m0-m5 全回归(Features 相关断言不受影响);node --check

## 验收结果(实施后)

- A+D:featLegend 图例 chips(色点+type ×n,按计数降序,仅现存类型)
  挂 Features 复选框后;addFeatureSpheres 渲染,无分子/解析失败/
  未勾选清空;aspirin 实测 "acceptor ×4 donor ×1 neg ×1 hydrophobe
  ×1 ring ×1",未勾选清空 ✓
- B:hydrophobe 灰 #94a3b8 → 琥珀 #a16207(FEATURE_COLORS 单源,pharm
  编辑器色点同源跟随);类型化默认 alpha(donor/acceptor/pos/neg
  0.45→0.55,ring 保持 0.45);两档亮度匹配逻辑不动;复选框 title
  重写(含单点模型与匹配亮度说明)
- 截图目检:图例五色点清晰、与选项行排版协调;琥珀球与碳原子区分
  开(此前灰色隐形);390px sw==cw==390;零 page error
- m0-m5 全绿 37/10/11/10/32/41;7 内联脚本 node --check;零引擎/
  wasm/算法改动
