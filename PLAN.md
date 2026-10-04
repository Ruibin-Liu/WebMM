# Plan: LigandLab 风格平台外壳(纯 UI 重构,能力零改动)——已完成

## 目标

按竞品截图(LigandLab,中文 LBDD 平台设计稿)的视觉语言,把
`app/platform.html` 的外壳重做:**深蓝左侧导航 + 顶栏(面包屑/版本/任务
指示)+ 四步流程指示器(stepper)+ 卡片化主区 + 徽章/漏斗条**。所有既有
功能、元素 id、JS 函数、worker、IndexedDB 模块**零改动**——只动壳、CSS
与少量展示层 JS。诚实适配:不做 QSAR/pIC50/适用域/云端库的假 UI。

## 范围(仅 app/platform.html + tests/cdp/m6_platform.test.js)

### A. 设计系统 CSS(追加层,不重写既有样式)
- 变量:侧栏 navy `#16202e`、主色 teal `#14b8a6`、底 `#f4f6f9`、白卡
  圆角 10px + 轻阴影;`.panel`→卡片、`.action-btn.primary`→teal。
- 徽章/风险标签 pills(绿/橙/灰)、stepper(圆点+连线)、任务条。

### B. 外壳 DOM:body → `.app-shell`(grid:侧栏 216px + 主列)
- `<aside class="sidenav">`:品牌、导航 6 项(内联 SVG 图标 + 文案)
  - 项目概览 `tabProject`(新,mode 'project')· 分子工作台 `tabSingle`
    (id 原样,按钮从 mode-tabs 迁入)· 批量处理 `tabBatch` ·
    相似性检索 `tabSearch` · 骨架探索(锚点:search+analogPanel)·
    检索历史(锚点:search+queryHistoryPanel)
  - 页脚:Workbench 链接、本地计算声明、存储徽章
- `.main-col` = 顶栏 + stepper + 既有 container
  - 顶栏:面包屑(随模式更新)+ 任务指示 `taskTicker`
    (MutationObserver 监听 searchStatus/analogStatus/batchStatus,
    运行态 teal 脉冲)+ RDKit/engine 版本(id 迁入,原处删除)
  - stepper 四步(可点击,诚实映射):结构导入→single · 候选检索→
    search · 骨架探索→analogPanel · 性质与精选→batch;switchMode 同步
    done/current 态

### C. 项目概览面板(新 `overviewPanel`,mode 'project')
- 统计卡:库规模/工作集(读既有 state)/引擎与 RDKit 版本;
  存储状态(navigator.storage);
- 项目导入/导出按钮**迁入**(`projectImportFile` 等 id 原样保留,
  setInputFiles 兼容);
- 诚实横幅(对应竞品"模型预测≠实验结果"):"本地确定性引擎——不含
  QSAR/ADMET 预测模型;分数用于优先级排序"。

### D. switchMode 扩展(增量)
- `projectMode` 布尔;single = !batch&&!search&&!project;
  overviewPanel 显隐;tabProject active;面包屑+stepper 同步;
  modeInputs/modeShown/MODE_SECTION 增 project 键。

### E. 窄屏
- <920px 侧栏缩为 56px 图标列;stepper 文字缩短。390px 无横向溢出。

## 不做
- 引擎/wasm/worker/IndexedDB 模块改动;新能力(ECFP4 导出、PAINS 等);
  右侧浮动详情卡(单分子工作台即检查器,映射声明);假 QSAR/ADMET 列。

## 验收
1. m6 既有 27 项全绿 + 新增 3 项(外壳渲染 6 导航项;overview 可达且
   统计卡+项目按钮在场;stepper 第一步点击切到 single)。
2. m0–m5、平台 Node 34/34 回归全绿;内联脚本 node --check;零 page error。
3. 目检:侧栏/顶栏/stepper/卡片/徽章符合截图视觉;窄屏 390px 无溢出。


## 实施与验收(完成)

1. **外壳**:body → `.app-shell`(216px 侧栏 + 主列);深蓝渐变侧栏 6 导航
   项(内联 SVG 图标;tabSingle/tabBatch/tabSearch id 迁入 + 新 tabProject;
   骨架探索/检索历史 = search 锚点 sideGoto);页脚 Workbench/GitHub/本地
   声明/存储徽章。顶栏 = 面包屑(随模式)+ 任务指示 taskTicker
   (MutationObserver 四状态源,运行态 teal 脉冲)+ RDKit/engine 版本(id 迁入)。
2. **Stepper**:四步可点击,诚实映射(结构导入→single · 候选检索→search ·
   骨架探索→analog 锚点 · 性质与精选→batch);switchMode 同步 cur/done 态;
   project 模式全中性(初版"全 done"语义错误已修)。
3. **项目概览(mode 'project')**:5 统计卡(库/工作集 pin/轮次
   [triageStore().projectApi().state.rounds 同源]/引擎/存储 estimate);
   项目导入导出按钮**迁入**(id 保留,setInputFiles 兼容);诚实横幅
   "本地确定性引擎——不含 QSAR/ADMET 预测模型";输入面板该模式隐藏。
4. **switchMode 扩展**:projectMode 布尔;modeInputs/modeShown/MODE_SECTION
   增 project 键;syncShellChrome() 统一驱动面包屑+stepper。
5. **控件协调 CSS 层**:teal 主按钮、中性描边 edit/history/spatial、
   红描边 clear、teal 描边 export(旧金黄/淡紫废除外观);select/range
   accent-color;表格表头着色+行 hover;dropzone 强化。
6. **窄屏**:<920px 侧栏缩 54px 图标列;390px 溢出 0px。

## 验收数字

- **m6 30/30**(+3:外壳 6 导航项+4 tab+stepper+ticker;overview 统计卡+
  项目按钮迁入+stepper 中性;step1 点击回 single);m0–m5 37/10/11/10/32/44
  全绿;平台 Node 34/34;7 内联脚本 node --check;零 page error;
  390px 溢出 0px;目检三轮(overview/search/窄屏)通过。
- 引擎/wasm/worker/IndexedDB 模块零改动;全部既有 id/函数保留。


## 视觉轮 2(继续调;已完成)

1. **sticky 缝隙修复**:topbar 实测 34.6px vs stepper top:46px 写死 → topbar
   定高 46px,滚动时 gap=0 实测验证。
2. **stepper 连线**:line1-3 id + done 态(已完成段 teal 填充,flex 可伸缩
   18-70px);syncShellChrome 同步。
3. **平台身份**:h1 "WebMM Workbench" → "WebMM Platform",副题改中文
   本地配体设计流描述(侧栏/顶栏/标题三级一致)。
4. **漏斗 chips**:searchResultStatus → teal 胶囊(内容不变)。
5. **探索器空态**:analogStatus 初始"待命 — 选定位点后 Explore,或
   aza-scan 一键骨架跃迁"(不含 m6 终态正则词)。
6. **行距**:analogPanel/pharmQueryPanel 按钮行 wrap+row-gap;separator
   字距微调。

验收:m6 30/30、m5 44/44、平台 Node 34/34;目检(search 视图:连线 teal/
chips 可读/标题正确/无重叠);零 page error。
