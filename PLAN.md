# Plan: 骨架探索 UI 三件(缩略质量 / 流程提示 / 结果导出)

## 1. 缩略图质量(用户:RDKit 先设定图大小,web 严格同尺寸)
- 根因:molThumb 单一 64×44 渲染被 CSS 拉到 140×96;小画布下 RDKit
  最小字号与键距比例失衡。
- 修:molThumb(smiles, w, h) 按尺寸缓存 —— 行内 128×88 渲染显示 64×44
  (2× retina),悬停弹层 400×272 渲染显示 200×136(严格同比 400/272=
  200/136);.thumb-pop svg CSS 200×136,弹层夹紧常量同步 210×146。

## 2. 流程标签悬停说明(详细、专业)
- 单位点交换卡副题 '位点 → 片段 → 形状漏斗 B·C·D':title 描述 A 枚举
  (SMARTS 位点检测 6 类末端基 → 片段拼接 → 价键/转子/质量过滤)→
  B 粗筛(10 构象 ETKDG+MMFF94s,刚性形状叠合 Tanimoto,流式)→
  C 精化(top-30 扩 30 构象,两轮系综取 max)→ D 柔性(top-10 MCS-
  Kabsch 核心 + 位置约束 MMFF 松弛重打分)。
- aza 卡副题:生物电子等排 C→N 说明(苯→吡啶/二嗪,单/双取代,
  规范 SMILES 去重,邻位 N-N 过滤,同漏斗,母体构象为参照)。
- 检索分段控件三钮:各自相似度度量说明(Tanimoto/ECFP4、子结构、
  构象系综形状 Tanimoto)。

## 3. 结果导出
- analogResultsWrap 表上方加 导出 CSV / 导出 SDF 两按钮。
- CSV(csvField 复用):SMILES,Name(枚举标签),ShapeT,FlexT,MW,
  RotB,SA。
- SDF:逐行 2D molblock 记录 + 数据字段(ShapeT/FlexT/SA/MW),
  $$$$ 分隔,RDKit 可读。

## 验收(新增 tests/cdp/m6_analog_ui.test.js)
- 悬停弹层:svg naturalWidth=400/height=272,显示盒 200×136;行内
  svg naturalWidth=128。
- 两个卡副题与三个分段钮 title 非空(≥60 字符,含关键术语)。
- 小跑一次探索后:导出 CSV 文件头+行数+数值列;导出 SDF 含 V2000 与
  $$$$ 计数 == 行数;零 page error。
- 回归:m6_platform/loop/pool/cancel + run.sh 全量绿。
