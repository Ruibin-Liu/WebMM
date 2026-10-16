# Plan: Platform 代码审查修正(已完成)

## 范围
- app/platform.html:删除逐字重复的 CSS 区块(1188–1358 = 1017–1187 副本);
  移除未使用的 JSME 加载脚本与 meta 文案中的 JSME(Ketcher 已替换)。

## 实施结果(完成)
- CSS 副本删除前 diff 复核零差异;删除后 design-v2 段仅存一份,级联顺序不变。
- jsme.nocache.js 加载与 meta 中 JSME 文案移除;openJSME/jsmeModal 内部名保留,
  vendor/jsme 保留(app/index.html 仍在用)。
- 验收:platform Node 39/39、m6 32/32、import/reload loop 10/10、
  jsME 网络请求 0、零 page error。
