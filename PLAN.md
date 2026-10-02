# Plan: GFN-FF 的 CH···π 行为 UI 提示(app 端标记级)

## 背景

上轮诊断:GFN-FF 把芳环 C-H···π 建模为吸引性"氢键"(HB 受体即低
qa 的 π 碳),优化后 H 可压在环面上方 ~2.4 Å——用户看到"氢往环内挤"
实为力场设计,MMFF94s 则无此项(2.98 Å 常规距离)。用户要求加 UI
提示避免误读为缺陷。

## 实施(app/index.html,纯标记+两行 JS)

1. GFNFF `<option>` 加 title(悬停可见完整说明)
2. 3D 面板工具条下方的选项行内加 `<span id="gfnffHint" class="hint"
   style="display:none;">`:引擎选 GFN-FF 时显示——"GFN-FF: aromatic
   C–H···π is attractive by design — H's may press onto ring faces
   (~2.4 Å); MMFF94s relaxes these contacts"
3. engineSel change 监听 + 初始同步,切回 MMFF94s/MMFF94 隐藏

## 验收(Playwright)

- 默认 MMFF94s 提示隐藏;切 GFNFF 显示;切回隐藏;文案含 π 说明
- 390px 提示显示时无溢出;m1/m2 抽查回归;node --check 内联脚本;
  零 page error

## 验收结果(实施后)

- Playwright 7/7:默认 MMFF94s 隐藏 → 切 GFNFF 显示(文案含 π 与
  MMFF94s 对照)→ 切 MMFF94 再隐藏;GFNFF option 悬停 title 在位;
  390px 提示显示时无溢出;零 page error
- 回归:m1 10/10、m2 11/11;7 内联脚本 node --check
