# Plan: 主区真重构 + 新设计语言 v2(不复用 workbench UI 组件)——已完成

## 承认问题

前几轮 = 侧栏/顶栏/stepper 套在**原面板堆叠**外;主区仍是 textarea→按钮行→表格
自上而下。本轮重构 `app/platform.html` 主区本身(index.html 保持经典页不动;
引擎/worker/存储零接触)。

## A. 单分子模式 → 三栏 dashboard(最大视觉转变)

```
┌──────────┬──────────────────────┬────────────┐
│ 参考配体  │ 2D 结构图卡           │ 检查器卡    │
│ (输入卡)  │ ──────────────────── │ 性质横条    │
│ textarea │ 3D 查看器卡(含能量)   │ Properties │
│ 绘制/拖放 │ ──────────────────── │ Rules      │
│ 历史      │ 构象系综卡            │ Alerts     │
└──────────┴──────────────────────┴────────────┘
```

- 实现:`#dashSingle` 三列网格 + `applyDashLayout(single)` 模式切换时
  **物理搬移**共享节点(textarea/#singleActions/#drop ← inputPanel;
  2D 子块/3D panel/#confPanel ← output;props 子块 → 检查器列),
  离开 single 模式按启动时记录的锚点还原(batch/search 共享这些节点,
  不能复制)。纯位置搬移,零行为改动;全部 id 保留。
- **性质横条(检查器头部)**:诚实指标条——QED↑(0-1)/ SA↓(1-10 反向)/
  MW(÷500 Lipinski 上限)/ cLogP(|x|÷5)/ TPSA(÷140)/ RotB(÷10),
  超限变琥珀;数据取自现有 descriptors 渲染路径,#propBars 由同一处填充。

## B. 检索模式 → 查询卡 + 库卡双列 + 筛选条

- `#searchBar` 内部重排:`.dash-search-row` 双列网格——左"查询"卡
  (query/mode/指纹/阈值/pharm/Search),右"配体库"卡(装载/Demo/上传/
  清除 + 状态行);useCurrent 并入查询卡。
- MW/cLogP/PAINS 过滤控件 + filterStatus 从查询行抽出,贴到结果卡头部
  成筛选条(仍只触发 renderSearchResultsCurrent)。
- 探索器/RGD/scaffold/药效团/工作集面板:卡片头统一(图标+标题),
  不改行为。

## C. 通用

- 卡片头组件(.card-head:SVG 图标 + 标题 + 说明);面板圆角/阴影已有。
- 窄屏(<1100px)三栏降级:输入/结构/检查器纵向堆叠;<920px 侧栏已图标化。
- 零 page error;m6 30/30(全部 id/函数不变);两页 7 脚本 check;
  390px 无横向溢出。

## 不做

- index.html(经典页刻意保留);引擎/worker/IndexedDB;表格列结构;
  假 QSAR/ADMET。

## 验收

m6 30/30 + 平台 Node 34/34 + m5(index 不受累)44→47 项全绿;
单分子全流程(process→embed→optimize→conformers)在三栏内工作;
检索全流程(库→查询→过滤→行点击回载单分子三栏)无断裂;
目检截图(宽/窄)+ 零 page error。


## 实施与验收(完成)

1. **单分子三栏**(290px 参考配体卡 | 弹性结构列[2D/3D/构象] | 340px
   检查器卡):`applyDashLayout` 模式切换时物理搬移 7 个共享节点
   (textarea/singleActions/drop ← inputPanel;2D 子块/3D 面板/confPanel
   ← output;props 子块 → 检查器),离开 single 按启动锚点逆序还原;
   **元素引用首查缓存**(排雷:二次调用重查空容器 → undefined →
   appendChild 抛错 → render 静默中断——症状:结构列不显/横条 0);
   空态(未处理分子)只显示输入卡,结构/检查器隐藏。
2. **检查器性质横条**:QED(0-1)/MW(÷500)/cLogP(|x|÷5)/TPSA(÷140)/
   RotB(÷10),teal 渐变填充,超限琥珀+⚠;render() QED 处填充
   (paracetamol:0.60/151/0.00/49/1)。
3. **检索双卡 + 筛选条**:initSearchCards 一次性重排——查询卡(1.45fr:
   query/mode/指纹/阈值/pharm/Search)+ 配体库卡(1fr:Load/Demo/
   Upload/Clear);MW/cLogP/标PAINS + filterStatus 抽到结果卡头部
   虚线筛选条;结果卡 = 图标卡头 + 漏斗 chip(searchResultStatus 迁入)。
4. 响应式:<1150px 三栏/双卡降单列;390px 溢出 0(实测)。

## 验收数字

m0-m6 37/10/11/10/32/47/**30** 全绿;平台 Node 34/34;7 脚本 check;
目检三张(boot/单分子三栏含横条/检索双卡+筛选条/窄屏堆叠)通过;
零 page error;引擎/worker/存储/id 契约零改动。


## 设计语言 v2(用户授权"不需要使用原 workbench 的任何 UI 元素")

1. **新组件系统**(CSS 层,.btn 系列**映射取代** .action-btn 旧皮肤):
   btn-primary(teal 实心)/ btn-ghost(描边)/ btn-danger(红描边)/
   btn-accent(teal 描边);统一圆角/按压反馈/焦点环。
2. **.field 字段模式**:小号大写标签在控件上方 + focus teal 环——检索卡
   全部控件字段化(查询结构/检索方式/指纹/阈值/药效团/形状选项)。
3. **分段控件**:检索方式 select → 相似/子结构/形状3D 分段按钮
   (segModeSet/syncSegMode;隐藏 select 保留契约——m6 程序化设值+
   onSearchModeChange 仍生效,seg 经钩子同步)。
4. **粘贴字段**:单分子输入 = 虚线等宽 textarea(focus 转实线 teal 环)+
   瘦身拖放条 + ::after 提示行;批量态经 body[data-mode] 区分皮肤。
5. **表格 v2**:小号大写表头/紧凑行距/tabular-nums;空态组件。
6. 检索行重构为 qgrid(auto-fit 字段网格)+ 检索按钮加大(检 索)。

验收:m6 30/30、m5 47/47、平台 Node 34/34;390px 溢出 0;
目检(boot/单分子 paste-field/检索分段+字段化)通过;零 page error;
id/函数契约零改动(index.html 经典页不动)。
