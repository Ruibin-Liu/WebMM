# Plan: Search/batch 主按钮收编 `.action-btn` 体系(纯 app 端)

## 背景

用户质疑 Search 标签 "Load library" 等按钮 → 实测(Playwright 量几何)
三类症状,均为存量、非 536652b 引入:

1. **图标失控**:btnBatchRun/btnLibLoad 标裸 `class="primary"`,游离于
   `.action-btn` 体系外 → svg 尺寸约束(`.action-btn svg{14px}` 等)
   不匹配,`<svg viewBox>` 无 width/height 渲染 74.8×74.8/65.8×65.8,
   按钮 107.6/98.6px(正常 ~33px),flex stretch 拉齐同行
2. **hover 白字白底**:`button:hover{#f8fafc}` 晚于 `button.primary`
   同特异性 → 6 个裸 primary 按钮 hover 白字近白底
3. **行高被长标签折行绑架**:`.action-btn{width:110px}` 固定宽下
   "Demo library" 折两行(48px)、"Build from current molecule" 折三行
   (63px),stretch 拉齐整行

用户第二轮指出正解:**Single 标签按钮的 CSS 才是参照**——
`#singleActions`(Draw/Edit=`action-btn primary` 等,图标/hover 全由
体系提供)与 `.ctl-group .action-btn{width:auto /* no fixed 110px —
labels differ */}`、`btnBatchCancel` 内联 `width:auto` 三个既有先例。

## 修复(根因:收编体系,替换第一轮的症状补丁)

1. 6 个裸 `class="primary"`(btnBatchRun/btnLibLoad/btnSearchRun/
   btnRGD/btnScaffold/btnPharmQ)→ `class="action-btn primary"`
   ——svg 14px、hover #1d4ed8、内边距/字重/阴影全部由既有体系提供
2. 撤掉第一轮症状补丁:svg 的 width/height 属性(类规则接管)、
   新增的 `button.primary:hover` 规则(死代码;`.action-btn.primary:hover`
   覆盖)
3. 折行标签按仓库既有先例(ctl-group 注释 + btnBatchCancel)加内联
   `style="width:auto"`:btnBatchRun、btnLibLoad、btnLibDemo(存量
   同病同修)、btnRGD、btnScaffold、btnPharmQ、药效团
   "Build from current molecule"(存量同病同修)——各所在行回单行
   33px;短标签(Search/Example/Auto core/Export CSV/Clear/
   Upload file)维持 110px 网格不动

**不做**(记录理由):`.version` 等宽脚注字体——全局样式决策另行
拍板;不改 `.action-btn{width:110px}` 全局规则(全 app 影响面)。

## 验收(Playwright,localhost:8901;desktop+390px)

1. 六钮 + 所在行(batch/library/RGD/scaffold/pharm/search)全部
   单行 33px、无拉伸;btnBatchRun/btnLibLoad svg 实测 14px
2. 六钮 hover 实测 `rgb(29,78,216)`(#1d4ed8)+ 白字(transition
   0.2s 结束后读数)
3. 390/1440 无横向溢出;零 page error;node --check 全部内联脚本
4. cargo test 295(引擎零改动 sanity)

## 验收结果(实施后)

- ✓ 全部行单行 33px:batch 118×33/svg14、library 四钮 33、
  RGD/scaffold/pharm/search 全 33(药效团行存量 63px 拉伸同愈)
- ✓ 六钮 hover 全 `rgb(29,78,216)`+白字(btnBatchRun 单独在 batch
  模式下实测——混序读数曾因按钮被 switchMode 隐藏而得静息色,测试
  伪影)
- ✓ 390/1440 无溢出、零 page error、6 内联脚本 node --check、
  cargo test 295 全绿
- 注:m3_site_nav 的 axe `nested-interactive`(#viewer3d)为存量
  (git stash 前后 A/B 同现,第一轮已验证),非本次回归
