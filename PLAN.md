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
