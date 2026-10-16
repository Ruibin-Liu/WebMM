# Plan: CSV 导出加固(公式注入中和 + RFC4180 引号)

## 背景(最初平台页审查发现 #5,一直未修)
- 两页共 9 处 CSV 导出器各自实现 esc():platform.html 3 处(exportResultsCSV
  的 q()/exportRGDCSV/exportScaffoldCSV——scaffold 版连 \n 都不查),
  index.html 6 处(exportResultsCSV/exportRGDCSV/exportScaffoldCSV/
  exportHistoryToCSV/exportBatchCSV + 一处内联)。
- 风险:①字段含换行破行;②公式注入——导入库的分子名以 = + - @ 开头时
  Excel/Sheets 打开即当公式执行(名字是用户导入数据,真实攻击面)。

## 改动
- 两页各加共享 `csvField(v)`:危险前缀中和(`'` 前缀,Office/Sheets 惯例)
  + 按需引号(含 , " \n \r 才引,引号转义 "" )。
- 9 处调用点全部替换为 csvField(调用点形如 `const q = csvField;` 最小 diff)。
- 新增 tests/cdp/m6_csv_export.test.js:两页 csvField 直接单测
  (公式前缀/逗号/引号/换行/正常名)+ platform 真实导出流
  (敌意名分子入库→sim 检索→exportResultsCSV 下载→文件内容断言)。

## 不做
- 不改列结构/文件名/状态文案;不动 m4/m5/m6 既有断言(朴素 split 兼容)。

## 验收
- 新测试全绿;m4/m5/m6/platform Node/loop 复跑全绿(两页回归)。
