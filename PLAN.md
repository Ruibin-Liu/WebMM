# Plan: 前端 Review 应修项修复(6 项,单轮打包)

## 范围(承接 frontend-design + web-artifacts-builder 双基准审查)

1. **CSS 同名选择器合并**:两页各删三块死码/被覆盖块(.ctl-group .mini-field
   死变体、组选择器中已被 v2 覆盖的 .mini-field:disabled、与 v2 逐字相同的
   .action-btn.primary:disabled 早期副本);platform 第二 :root 加注释(并存
   token 族,非死码,防未来误判)。
2. **焦点环**:两页 textarea:focus 补 box-shadow(outline:none 单换边框色
   不达 WCAG 2.4.11;探针实证 #input 焦点不可见 1→0)。
3. **390px 溢出**:index.html .topnav .links 加 flex-wrap(探针 47px→0)。
4. **prefers-reduced-motion**:四页 style 尾部注入减弱动效块。
5. **全局错误面**:两页挂 unhandledrejection → showError(#error 本有
   role=alert/aria-live)。
6. **本地部署暂存**:scripts/stage-pkg.sh(清坟场 + wasm-pack + index/
   playground/vendor/app 全量暂存,镜像 pages.yml);package.json "stage"
   脚本;AGENTS.md 构建命令行同步。

## 不做(审查记录在案)
- landmark/region(axe moderate,低价值);platform prompt/confirm(项目
  已记 v1.1);内联 onclick(内部 id 无现实风险);conf/shape 农场持久化
  (文档化范围决策)。

## 验收
- 前端探针复跑:focus 不可见 0、390px 溢出 0、零 page error。
- run.sh 全量 19/19(m6_projects 首轮导航竞态 flake,单跑 12/12 + 全量
  复跑 19/19 确认)。
