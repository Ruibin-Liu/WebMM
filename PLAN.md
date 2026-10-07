# Plan: 分子工作室重设计 —— 彻底去除 workbench UI 元素(ui/studio-redesign)——已完成

## 范围(用户指令:完全去除旧 UI 元素,重新设计;逻辑/契约不动)

1. **分子工作室卡**:#output 首面板重写——页签(三维/二维)+ 收纳式工具栏
   (嵌入3D/engine/优化/构象 + "高级"/"导出"下拉:构象参数/着色/特征 +
   全部导出按钮);3D 视图居中放大;能量表转折叠;2D 导出芯片行/控件摊开
   全部消灭。dashEls 重构(studio-main→结构列,studio-inspect→检查器列,
   身份面板经 #smiles.closest 全局定位入检查器)。
2. **结果行 2D 缩略图**(workbench 从未有):molThumb(SMILES→64×44 SVG,
   Map 缓存)入检索/探索器/批量三表;**排雷:SVG 杂原子 <text> 污染
   children[1].textContent 契约**——m6 四处读取改 .sm/.nm span 优先。
3. **头部块/页脚删除**:h1/副题/徽章/©/Playground 引导全部移除;
   Playground/Demo/查看源码入侧栏足部(fixSiteLinks null 守卫既有)。

## 审计结论

- 静态文本 224 条,~80 条英文残留(批量区 separator/按钮、库卡按钮、
  单分子面板标题与控件、探索器/RGD/scaffold/药效团/工作集/历史面板、
  顶栏徽章、Playground 引导)。
- **动态状态串是 m6/m5 测试契约**('55'/'at T'/'restored'/'aza-hop'/
  'pinned'/'contain the core' 等)——本轮**不动**,定位为"系统消息"
  (双语可接受);静态标签全量中文化(技术词 SVG/PNG/InChI/MW/flex 保留)。
- 三个 legacy separator("-- or --" 等)、旧面板标题(strong 裸文本)、
  部分旧按钮类需 v2 化。

## 改动(app/platform.html only;逐项先 grep 测试契约)

1. **静态标签中文化**(约 60 处):批量(Run batch→批量处理/Cancel→取消/
   3D 复选文案/separator 重写);库卡(装载库/Demo 库/清除/上传文件);
   单分子(绘制/历史/二维结构/性质/更多性质/规则/更多规则/结构警报/
   三维结构/嵌入 3D/优化/构象/按部分电荷着色/特征/导出系综 SDF/
   能量标题/查看器占位/分子式/规范 SMILES/复制/表示);徽章(本地运行/
   查看源码);Playground 引导句。
2. **旧面板 v2 化**:Neighborhood explorer→骨架探索卡、RGD、scaffold、
   Pharmacophore Query→药效团查询、Working set→工作集、Queries→检索历史
   ——.card-head(图标+中文)+按钮 btn 化(Screen library→筛选库/
   Build from current→由当前分子构建/Auto core→自动核心/Analyze
   scaffolds→分析骨架 等)。
3. **separator → 现代分隔**(CSS:细线居中小字)或删除(卡片结构已替代)。
4. **模态框**(源码/历史):v2 卡片样式 + 按钮中文化。
5. 探索器表头/操作列(★/pin)保持;表格 v2 全局样式已就位。

## 安全检查

- 每条待改文本先 grep m6/m5/m4-m0(index 不动,但同源文本)确认无断言;
  命中的(如 aria-label、占位串)只改可见文本不动 JS 引用值。
- 涉及 JS 读取的标签(label for)保留 for/id 关系。

## 验收

m0-m6 + 平台 Node 全绿;7 脚本 check;目检六状态(boot/单分子全流程/
检索全家族面板/批量/概览/模态);390px 溢出 0;零 page error。


## 实施与验收(完成)

1. **静态标签中文化 ~70 处**(逐条先 grep 测试契约):批量区(批量处理/
   取消/3D 复选/提示/批量结果);库卡(装载库/Demo 库/上传文件/清除);
   单分子(绘制/历史/二维结构/性质/规则/结构警报/三维结构/嵌入 3D/优化/
   构象/特征/导出系综 SDF/能量标题/查看器占位[HTML+JS 两处]/分子式/
   规范 SMILES/复制/表示);徽章(本地运行/查看源码);Playground 引导;
   探索器族(骨架探索/R-Group 分解/骨架频次/药效团查询 四卡头 v2 化 +
   探索/用当前分子/筛选库/由当前分子构建/自动核心/分析骨架/取消 +
   片段类别与作用域下拉全库/当前结果);表格表头(名称/分子数/占比/
   结果,RGD 动态表头同步);行提示三处;Export×5 → 导出。
2. **保留边界(文档化)**:动态状态串 = 测试契约 + "系统消息"语义
   ('55'/'at T'/'restored'/'aza-hop'/'pinned'/'scaffold'/'/3 molecules'
   等)保持英文;技术词(SMILES/MW/Tanimoto/ECFP4/flex/InChI)保留。
3. **修正过程排雷**:按钮文本带换行缩进——`>X<` 锚点失配,改按
   行内容锚定;三次补丁分批落盘。

## 验收数字

m0-m6 37/10/11/10/32/47/30 + 平台 Node 34/34 全绿(服务器未起的两次
假阴性已甄别);7 脚本 check;390px 溢出 0;目检四张(单分子/检索
家族×2/批量)通过;零 page error。


## 工作室轮验收(完成)

m6 30/30(4 处读取点适配缩略图)、m0-m5 37/10/11/10/32/47、平台 Node
34/34;7 脚本 check;390px 溢出 0;目检:工作室卡(页签/工具栏/3D 居中/
检查器)、缩略图行(55/55)、头部页脚无残留;零 page error。
排雷存档:命名 IIFE 外部不可见(检测误报);身份面板在 #output 同级
非子级(closest 全局定位取代子级过滤);dashEls 的 viewer3d.closest
('.panel') 在工作室结构里命中工作室面板自身(显式六元素表取代)。


## 排查优化轮(继续;已完成)

审计法:六状态截图(高级下拉/历史模态/RGD+探索器/药效团/批量/JSME)
逐张目检 + 文本盘点。修复清单:

1. **页脚真身**(crosslinks + copyright 两个 div,非 <footer> 标签——
   上轮正则没打着)+ 遗留 year 独立脚本(页脚删后 null 崩,pageerror
   抓获)→ 全删。
2. **JS 按模式设置的 placeholder**(批量/检索/单分子三态 + 查询框两态)
   → 中文——静态盘点抓不到的动态英文。
3. 探索器/RGD/骨架/药效团作用域与容差 select **剥离内联旧样式**
   (CSS .panel select 统一接管);过滤条 MW/cLogP 数字输入样式化;
   复选框全局 accent-color=teal。
4. 标签:身份表 Name→名称;工作集 ↶undo→撤销/hide excluded→隐藏已排除;
   批量过滤器 全部/Lipinski 通过;药效团说明段中文;JSME 模态 取消/
   应用到工作台;源码模态 Copy→复制。
5. 间距与容器:检索历史/工作集面板底色+边距;批量态拖放条瘦身;
   模态全局升级(blur 背景/14px 圆角/渐变阴影/页脚右对齐)。

验收:m0-m6 37/10/11/10/32/47/30 + 平台 Node 34/34;6 脚本 check;
390px 溢出 0;零 page error;批量 placeholder 实测中文。
视觉审计器误判甄别:过滤器/取消按钮只在运行态出现,截图态"缺失"
非缺陷;placeholder 为 JS 动态设置,静态盘点盲区。


## 续磨轮(JSME 包壳 + 批量列选择器;已完成)

1. **JSME 模态包壳**:标题"绘制结构"+副题;画布居中约束 780px;页脚
   .modal-footer 右对齐(取消红描边/应用 teal 实心)。画板本体第三方
   组件不动(文档化边界;空输入 ERROR 为其自有行为)。
2. **批量列选择器**(新 UI 元素):card-head 化批量结果卡(batchStatus
   迁入 card-sub);"列"下拉六项可选列(分子式/TPSA/HBD/HBA/RotB/
   E MMFF94s);**排雷:初版 has-hide 单类会一刀切隐藏全部六列**——改
   per-column hide-N 类(nth-child 稳定因行模板固定);实测全隐后
   9 列舒适密度 + 缩略图 + sticky 表头正常。
3. 历史模态 Clear All → 清空全部。

验收:m0-m6 37/10/11/10/32/47/30 + 平台 Node 34/34;6 脚本 check;
390px 溢出 0;零 page error。


## Ketcher 替换 JSME(已完成;platform 绘制模态)

1. **分发探明**:npm 无免构建单文件——ketcher-standalone(Indigo WASM
   服务,21MB ESM)+ ketcher-react(UI,需 React)分体;**一次性 esbuild
   打包为 IIFE 提交 vendor**(app/vendor/ketcher/:31MB bundle + 179KB
   CSS + Apache-2.0 全文),运行时零 CDN,与 JSME 同为提交库文件。
   Indigo 单线程(零 SharedArrayBuffer/crossOriginIsolated)——
   GitHub Pages/裸 http.server 无 COOP/COEP 也可用。
2. **入口 API**:WebMMKetcher.mount/getMolfile/setMolecule/unmount;
   **排雷①:banner `var process=…` 覆盖页面全局 process() 函数**
   → 守卫式 shim(typeof 检查,不覆盖既有);**排雷②:setMolecule 在
   onInit 前调用被静默丢弃**(画布空)→ readyPromise 队列,所有 API
   await onInit。
3. **接线**:绘制模态保留 id/函数名(openJSME/applyJSME/closeJSME),
   内部换 Ketcher 懒加载(首次 31MB,状态行提示);apply = getMolfile
   → input → process(错误不关模态防丢);失败态显示在模态状态行。
4. LICENSE 第三方清单 +5(Ketcher/Indigo Apache-2.0)。

验收:模态内编辑器完整(工具栏/元素面板/模板)+ 当前分子载入画布
(视觉确认 paracetamol)+ apply 闭环(molblock→process→2D 重渲染→
模态关闭);m0-m6 + 平台 Node 全绿;零 page error;零外联请求。


## 细节轮(已完成)

审计面:概览(带数据)/历史模态内容/源码模态/色权重弹层/390px 高级下拉。
修复:

1. **真 bug:概览工作集卡 "undefined pinned · undefined excluded"**
   ——getPins/getExcludes 返回对象(Object),原代码取 .length →
   Object.keys().length;实测 pin 后 "1 pinned · 0 excluded"。
2. 历史模态:空态中文(empty-state 组件);条目卡新 token(悬停 teal
   边+阴影)、SMILES 等宽字体灰底 pill;头部按钮层级化(导出 CSV
   teal 描边/清空全部 red 描边);删除 X 为 hover 显示(既有设计,
   视觉审计误报甄别)。
3. 色权重弹层 5+1 排布与 z-index 目检通过(无裁剪/无 sticky 冲突)。

验收:m0-m6 + 平台 Node 全绿;390px 溢出 0;零 page error。


## 侧栏审美重做(已完成)

1. **品牌区**:teal 渐变圆角磁贴(分子六边形 glyph,内嵌高光+外投影)
   + WebMM Platform 字标(MM teal,PLATFORM 小号大写副行)。
2. **导航分组**:工作区(项目概览/分子工作台/批量处理)/ 发现
   (相似性检索/骨架探索/检索历史)——小号大写组标签。
3. **激活态**:teal 13% 药丸底 + teal 文字/图标 + 左缘发光指示条
   (3px 圆角+glow);hover 中性提亮;过渡 0.15s。
4. **背景**:深海军蓝双色——右上 teal 4% 径向辉光 + 纵向渐变 +
   右缘 10% 分隔线(取代纯平面渐变)。
5. **页脚**:状态卡(呼吸灯 teal 圆点 + sideStorage 文本)+ 链接行
   + 本地计算注记。
6. **收窄 rail(56px)**:只留磁贴/图标/指示条/圆点;组标签/文字/
   链接全隐。侧栏 CSS 块整体重写;`.side-item`/tab id 契约保留
   (m6 30/30 零改动)。

验收:m0-m6 + 平台 Node 全绿;宽/窄/激活态三截图目检通过
("production-ready polish");500px 溢出 0;零 page error。


## 侧栏图标对齐与间距修正(用户指出;已完成)

放大 2× 裁剪逐项审计证实:分子工作台图标过重偏高 / 批量处理过轻
偏下 / 骨架探索偏空,行距偏紧。修复:

1. **图标集统一**:六枚全部 stroke 1.8 + round caps;分子工作台换
   六边形+节点(与品牌磁贴同构);批量处理=带前导圆点的列表行
   (补重量);骨架探索=角括号+实心六边形核;概览格/检索/时钟微调
   viewBox 填充;svg display:block 消基线缝。目检复验:同轴同光学
   尺寸,检索镜面再 -4% 微调。
2. **间距**:nav gap 2→4px;条目 padding 0.56→0.68rem;组标签上距
   0.85→1.05rem;图标 17→18px。

验收:m6 30/30;复检评语"统一成功,无进一步动作必要"。


## 图标移除 + 回归 workbench 配色(用户指令;已完成)

1. **图标全去**:六枚侧栏导航 SVG 移除(文本导航);品牌磁贴删除,
   回归 workbench 文字品牌(WebMM 蓝色 MM + PLATFORM 小标)。
2. **配色回归 workbench**:--ll-accent=#2563eb / --ll-accent-2=#1d4ed8
   (原 teal 变量重定义,CSS 变量引用全站自动换装);字面量 sed 全文
   (含内联属性与 JS heat/pinned 染色):rgba(20,184,166→37,99,235)、
   #2dd4bf→#3b82f6、#f0fdfa→#eff6ff、#b5d8d3→#bfdbfe、#0f766e→
   #1d4ed8 等;热力单元格/性质横条/分段控件/stepper/药丸态全部转蓝。
3. **侧栏浅色化**:白 92% + blur(与 workbench topnav 同工艺),右缘
   1px 边线;hover #f1f5f9;激活 #eff6ff 药丸 + 蓝字 + 蓝指示条;
   状态卡浅底蓝点;组标签 slate。
4. **窄屏**:图标既除,rail 方案不可行 → 侧栏转**横向滚动芯片条**
   (static,随页滚动;组标签/品牌/页脚隐);**排雷:sticky top:46
   在未滚动时把条带推进 topbar(z40)底下完全遮挡**——改 static,
   未滚动时 条带0→topbar42→stepper112 层序正确。

验收:m0-m6 + 平台 Node 全绿;390px 溢出 0;宽/窄截图目检通过
(浅侧栏"production-ready theming",零 teal 残留)。


## CSS 误删修复 + 全面普查(用户发现"结构导入大图标";已完成)

**根因**:图标移除轮的窄屏媒体块替换用了**区间切片**(`首个 @media 920 →
design-language v2 标记`),把两者之间的全部 CSS 连带删掉——round 2/3
(控件协调/viewer3d/能量卡/批量 sticky 表头/滚动条)+ **内容重构全段**
(.dash/.dash-card/.card-head[含 svg 16px 规则]/.pbar/.qgrid/.filter-strip
等)——card-head 图标失去尺寸约束爆到 1052px(即用户所见"大图标")。

**修复**:从 1c02ab5 完整重建——老 CSS 全量 → 色板 sed → 换浅侧栏块 →
**定点**替换 920 媒体块(不再区间切片)→ 删死规则 → 拼回 v2 之后幸存段
(v2/studio/thumbs/audit/hist)。定义唯一性核对(.dash/.pbar/viewer3d)。

**普查**(用户要求的全面检查):SVG 尺寸扫描 boot/处理后 0 个超限(内容
区除外);六状态截图(boot/单分子/检索/批量/概览/390px)目检:三栏/
工具栏/横条/蓝主题全部就位;390px 溢出 0;零 page error;m0-m6 + 平台
Node 全绿。

**教训存档**:①区间切片替换必须先核对切片内含的无关规则;②CSS 大改后
加"图标尺寸普查"回归探针(本次即由尺寸普查定位)。


## 侧栏层级修正(用户指出;已完成)

问题:组标签与条目同左边线、字号差小(0.62 vs 0.86rem)、颜色相近
(#94a3b8 vs #64748b)——七个元素读起来平级。修复(三重层级线索):

1. **缩进**:条目 margin-left 0.8rem(width calc 补偿),组标签贴左缘;
2. **字号差拉大**:标签 0.58rem vs 条目 0.9rem;
3. **颜色拉开**:标签 #b3bfcc(更浅)vs 条目 #3f4c60(更深)。

验收:放大目检"层级立即可读,无歧义";m6 30/30。


## 概览/工作集中文收尾(用户指出;已完成)

- 概览卡:工作集 pin→工作集钉选;detail 'N pinned · M excluded'→
  'N 已钉 · M 已排除';'DAG 轮次 · hit-as-query ⇄'→'轮次 DAG · ⇄
  命中即查询'。
- 工作集计数行(renderWorkSet)同步中文化——**m6 五处断言同步**
  (/1 pinned/×4、/1 excluded/×1 → 已钉/已排除)。
- 行悬停 title:pinned to working set→已钉入工作集;excluded from
  triage→已排除出甄别。
- 边界维持:瞬时状态行(检索/批量/探索器进度)与 CSV 列名
  (consensus_mean_rank/pinned/note/excluded 为导出数据格式)保持
  英文,已文档化。

验收:m6 30/30(断言同步后);m0/m5 抽检绿;实测概览卡五标题全中文。
排雷:python heredoc 语法错是编译期——整段不执行,勿以为前半已落盘。


## 侧栏层级 v2(中文正确的层级语法;已完成)

**认知修正**:小型大写+宽字距是拉丁文的层级语法,对中文无效(0.58rem
中文不可读)。v2 用对中文有效的三要素:

1. **组标签 0.7rem 可读灰 + 右延细线**(hairline rule,给组一个
   "视觉地板",标签不再悬浮);
2. **条目 0.95rem 深色,缩进 1.15rem**(≈18px,明确子列);
3. 组间呼吸 1.15→1.4rem(复检 nit 落实)。

验收:放大目检"层级一眼可辨、细线是帮助不是杂讯、缩进恰好在
甜点";m6 30/30。


## Stepper 移除(用户质疑"一定要常驻吗";已完成)

判定:不必常驻——点击目标与侧栏完全重复,占一条 sticky 横条,纯装饰
(LigandLab 模仿残留)。移除:

1. markup(#stepperBar 整块)、syncShellChrome 的 step 同步与 stepGo、
   CSS(.stepper/.step/.step-line + 媒体引用)。
2. m6 三断言改写:shell 检查改"stepper removed"守卫(元素不存在);
   "step1 点击回 single" 改为"侧栏 tabSingle 点击回 single"。
3. 验证:滚动后 topbar 0..46 无缝;390px 溢出 0;m6 30/30;零 page error。


## 侧栏终极简化:平铺(用户三回不满意;已完成)

判定:三轮回退说明**分组本身(工作区/发现)不匹配用户心智模型**——
任何编码它的缩进/字号/细线都读作噪音。执行预留方案:

- 两个组标签删除;六项**平铺**(同字号 0.92rem、同缩进 0、同距
  gap 3px + padding 0.64rem);激活药丸 + 左缘指示条不变。
- m6 零改动(side-item 契约未动);390px 溢出 0。
- 放大目检:"零缩进、均一尺寸间距、无分组残留"四项全过。

教训存档:当用户反复否定同一元素的多种视觉编码时,答案往往是
删掉该元素承载的**概念**(分组),而不是换第三种画法。


## 逻辑审计(用户问"其它地方有不合逻辑的吗";已完成)

**发现(概念↔行为不匹配):**

1. **侧栏六项不等价(实锤)**:骨架探索/检索历史不是模式,是相似性
   检索模式内的滚动锚点——实测点"骨架探索"高亮的却是"相似性检索";
   从批量模式点它还会先切模式。与分组问题同构:视觉平级、行为不平级。
   → **修复**:删除两项(sideGoto 一并清除),导航 = 4 项 = 4 个模式,
   所见即所得;面板仍在检索模式内且卡片醒目。
2. **面包屑常量前缀**:"WebMM 项目 · X"——没有项目切换,"WebMM 项目"
   是永不变的噪音。→ **修复**:只显示节名。
3. **记录不改**(收益/扰动比低):药效团查询/骨架探索面板依赖单分子
   视图的分子却住在检索模式(历史布局,搬迁牵测试);批量处理是工具
   而非流程节点但作为模式成立;瞬时状态行英文(文档化边界)。

m6:shell 断言 items>=6 → ===4。验收:m6 30/30;导航四项/面包屑
纯节名实测;390px 溢出 0。


## 逐模式内部审查轮(用户指令;已完成)

**程序探针**(:mode 可见性矩阵 / Enter-拖放行为核实 / 标签 grep):
- 悬空承诺排查:输入卡提示"Enter 处理"**属实**(keydown Enter→
  maybeProcess);拖放**有模式感知**(batch 填框不自动跑 / search 填框
  +提示装载)——非缺陷。
- 四模式 × [drop/singleActions/batchBar/searchBar] 可见性矩阵全部
  正确(无泄漏控件)。

**目检修复(静态标签 + 对齐):**
1. RGD 提示与核心 placeholder 中文化;药效团 ± 容差 (Å)/构象数;
2. 探索器空态措辞("在「位点」下拉选择…"——与实际控件对应);
3. 骨架表 分子数/占比 统一右对齐(tabular-nums);状态行/底部提示
   呼吸间距。

**边界维持**:动态状态行(含骨架汇总 "19 scaffolds across…")保持
英文——与检索 chip 同一文档化边界,不在本轮扩张。

验收:m6 30/30、m5 47/47;骨架卡复检三项生效;零 page error。


## B 方案落地:骨架探索独立成模式(用户选定;已完成)

概念摆正:相似性检索/药效团查询 = 筛库(留检索模式);骨架探索 =
从当前分子生成新候选(独立第五模式)。

1. **模式管线**:analogMode 布尔;modeInputs/modeShown/MODE_SECTION/
   SECTION_LABEL 增 analog;switchMode 五路(data-mode/面包屑/tab
   active/面板显隐);该模式隐藏共享输入面板。
2. **侧栏**:tabAnalog 排在相似性检索**之后**(先筛选后生成的顺序);
   导航 = 5 项 = 5 模式。
3. **药效团查询卡**迁至检索模式内结果卡之后、检索历史之前
   (DOM 顺序 compareDocumentPosition 验证)。
4. **探索器空态**措辞:输入来源 = 「分子工作台」的当前分子。
5. m6 同步:items===5、tabs+tabAnalog、M3 段 switchMode('analog');
   排雷:crumb 模式解析初版漏 analog 分支(显示"分子工作台")。

验收:m0-m6 37/10/11/10/32/47/30 + 平台 Node 34/34(m1 一次假阴性
=/tmp/caff24.sdf 夹具被清,按 README 配方重生 24 原子后过);
390px 溢出 0;零 page error。


## 激活指示条移除 + 左边距收紧(用户确认所见后定位;已完成)

**定位过程**:侧栏五项 computed 实测零缩进(3× 目检复核)——"缩进感"
的真实来源是**激活态左缘指示条**(3px 蓝线钉在侧栏最左缘,药丸内缩
12px,激活项被读作"缩进")。修复:

1. `.side-item.active::before` 指示条删除——激活 = 药丸 + 蓝字加粗,
   无任何附加装饰;
2. 条目水平 padding 0.7→0.5rem、品牌 0.6→0.4rem:文字左缘 ~20px,
   与品牌字对齐(实测 textLeft=20);
3. 复检:"零指示条、文字与品牌同线、无任何缩进观感、ready to ship"。

m6 30/30。


## 侧栏左缘像素级对齐(用户"还有差距";已完成)

**差距坐实**(亚像素测量):品牌字 18.4px / 条目 20px / 链接 17.6px /
注记 17.6px——四条不齐的左缘,最大差 2.4px。修复:品牌 padding-left
0.4→0.5rem、.side-links/.side-note 0.35→0.5rem——**实测四类文本全部
= 20.0px**;3× 目检"像素级干净,零可检偏移"。m6 30/30。


## 侧栏真凶:mode-tab 的 all:unset(用户点破"字号字体不一样";已完成)

**实锤**(五项全量测量,纠正此前只量第一项的盲区):项目概览 14.72px/
左缘 20px;其余四项 **16px/左缘 12px**。根因:四项带 mode-tab 类,
遗留规则 `.mode-tab { all: unset; }`(×2 处)把 .side-item 的字体与
padding **全部重置**(回落默认 16px、padding 0);tabProject 恰无此类。
上轮"四类文本 20px"的测量只量了 items[0]=tabProject——测量盲区存档。

**修复**:删除两处 `.mode-tab { all: unset; }`。五项全量复测(两种
激活态):左缘 20 / 14.72px 全一致;目检"除颜色/药丸外零差异"。

m6 30/30。


## 侧栏条目间距放宽(用户指出挤;已完成)

item padding-block 0.64→0.72rem;nav gap 3→6px;顶部 margin 0.35→
0.45rem。实测四条间隙均 6px;目检"呼吸充分、节奏均匀、对齐字号不
受影响"。m6 30/30。


## 排版转科学计算惯例(用户指令;已完成)

1. **正文字体**:system-ui → `"Helvetica Neue", Helvetica, Arial,
   "PingFang SC", "Hiragino Sans GB", "Microsoft YaHei", sans-serif`
   (MATLAB/Jupyter 一系 + 中文回退);
2. **等宽栈**统一:`"SF Mono", SFMono-Regular, Consolas,
   "Liberation Mono", Menlo, monospace`(Consolas 优先覆盖 Win);
3. **去 SaaS 风微排版**:.f-label 与 .panel table th 的小号大写+宽字距
   → 普通大小写 0.75rem/600(×2 处副本同改);品牌 PLATFORM 副标同。
4. 实测:body=Helvetica Neue;标签/表头 transform:none/12px;目检
   "RDKit/Jupyter 类工具视觉语言,中文字形清晰,混排协调"。

验收:m0-m6 + 平台 Node 全绿;390px 溢出 0;零 page error。


## 概览"偶发打开工作台"修复 + 检查器长值截断(用户报告;已完成)

1. **复现锁定**:概览模式下**任何重载**(手动刷新/项目导入后的显式
   reload)都落回单分子模式——模式从不持久,即用户所见"项目概览有时
   打开分子工作台"。修复:switchMode 写 sessionStorage('webmm-mode');
   启动恢复(URL 深链 ?molecule/#hash 优先单分子;非法值拒绝;
   try 守卫私隐模式)。实测:概览重载→项目概览,骨架探索重载→骨架探索。
2. **检查器"表示"长值截断**(选截断方案,弃整页下移方案):规范
   SMILES/InChI/InChIKey/Murcko 单元格 max-width 168px + ellipsis;
   复制按钮读 textContent = **完整值不受影响**(实测 fullLen 65、
   截断生效、copyStillFull)。

验收:m0-m6 + 平台 Node 全绿;零 page error。
