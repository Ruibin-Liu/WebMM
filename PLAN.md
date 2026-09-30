# Plan: Workbench 提示文字字体归位——静态说明散文从 mono `.version` 拆出为 sans `.hint`(纯 app 端)

## 背景

用户质疑 "click a row to load it into the single-molecule view" 与
"pick 2–6 features, …" 等文字的字体。审计 `.version`(ui-monospace
栈)全部 16 处用法,两类混居:

- **mono 合理(10 处,不动)**:版本行(832)、动态状态读数
  searchStatus/status3d/batchStatus/searchResultStatus/rgdStatus/
  scaffoldStatus/pharmQStatus、confChartAxis 数字轴、pharmQDist
  距离矩阵(white-space:pre 必须等宽)
- **mono 误用(6 处,静态说明性散文 → sans)**:batchInfo 功能摘要
  (853)、"click a row…" × 3(1147/1161/1196)、RGD 说明(1180)、
  药效团说明(1218)

字体原则:静态引导/说明文字 = caption,应用正文同款比例无衬线
(body 的 system-ui 栈)+ 小字号 + 弱化色;等宽只留给代码/SMILES/
版本号/状态读数/对齐数字。仓库既有先例:`.spiral-note`(0.75rem 弱化
sans,继承 body)就是正确形态。

附带收益:batchInfo 是 flex 行内 span,`.version` 的 margin-top:0.5rem
把它压低于同行按钮;换类后自然对齐。

## 修复(app/index.html,零引擎/零 API)

1. CSS 新增 `.hint { color: var(--muted); font-size: 0.75rem; line-height: 1.45; }`
   (sans 靠 body 继承;不带 margin——各用法自带内联 margin 或在 flex 行内)
2. 六处 class="version" → class="hint":batchInfo、三处 click-a-row、
   RGD 说明、药效团说明(id/内联样式/aria 不动)
3. 不动:其余 10 处 `.version`、`.spiral-note`、状态行 aria-live 语义

## 验收(Playwright,localhost:8901)

1. 六处计算样式 font-family 不含 monospace、含 system-ui;字号 12px
2. 十处保留处仍 ui-monospace
3. batchInfo 与同行按钮基线对齐(top 差 ≤2px)
4. 截图目检 batch/search(含检索结果)/pharm 面板;390px 无溢出;
   零 page error;node --check 内联脚本;m4/m3 回归;cargo test sanity

## 验收结果(实施后)

- 计算样式实测:6 处 .hint 全部 system-ui 栈(无 monospace)、12px;
  10 处保留 .version 仍 ui-monospace;batchInfo 与 Run batch 同行
  top 差 = 0(换类顺带消掉 .version margin-top 造成的下沉)
- 截图目检:batch/search/pharm 面板说明文字为比例 sans、版面干净
- 390px 三模式 scrollW==clientW==390;零 page error
- 回归:m4 7/7;m3 9/10(axe nested-interactive 存量,文献在案);
  7 内联脚本 node --check;cargo test 295 全绿(引擎零改动 sanity)
