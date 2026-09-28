# Plan: 形状检索颜色力场 v1 —— Color Tanimoto(ROCS 隐式颜色规则)+ Combo 评分

## 背景

3D 形状 v1 的已知短板:纯 shape 对化学特征不敏感(glucose 会排在含氢键
特征的分子之上)。ROCS 的实践是 shape + color 组合打分。本增量给 Shape
(3D) 检索加颜色层,纯 app 端(不改引擎/版本):

- **特征感知(页面 JS,规则表)**:六类特征位点,多标签不互斥
  - donor:N/O 带氢([NX3;H1,H2]、[OX2;H1])
  - acceptor:中性 N/O(排除阳离子;含酰胺 O/吡啶 N)
  - pos:formal charge ≥ +1;neg:formal charge ≤ −1
  - hydrophobe:与杂原子不相邻的 C/F/Cl/Br/I
  - ring:芳香环质心(单点,等价 ROCS ring 特征)
  - 权重 v1 统一(donor/acceptor/hydrophobe/ring/charge 全 1.0,文档化
    简化;ROCS 官方权重属闭源实现细节)
- **颜色重叠(页面 JS)**:同型位点对的原子高斯重叠(复用引擎同款
  GCI=2√2/逐元素 α 公式,纯解析 10 行);Color Tanimoto =
  O_c/(O_cA + O_cB − O_c);Combo = Shape T + Color T
- **流程**:shape_align_wasm 拿到最优姿态变换 → JS 对目标特征位点施加
  同一变换 → 颜色打分(对齐仍只优化 shape,combo 联合优化为后续,如实
  记录)
- **UI**:Shape 结果表加 Color T、Combo 两列(按 Combo 降序开关:列头
  点击切换排序键,默认 Combo)

## 验证

1. **规则表对照**(scripts/color_rules_check.py):donor/acceptor/pos/neg
   位点集合 vs Python RDKit fdef 特征(Donor/Acceptor/PosIonizable/
   NegIonizable 家族)在 LBDD refs 51 分子 + demo 库上的一致率;差异
   逐类列出(hydrophobe/ring 定义本就不同,不参与对照,文档说明)
2. **单测式验收(Playwright)**:自匹配 Color T = 1;对称性;色零分子
   (己烷)Color T 分母保护;combo = shape + color
3. **端到端**:aspirin 查询 → salicylic/glucose 排序变化符合化学直觉
  (带 COOH/酯受体的分子 combo 提升);390px;零 page error;sim/sub
   回归
4. 门禁:289 测试、clippy(已对齐 CI 1.98)、fmt;纯 JS 无需 wasm 重建

## 边界与不做

- 不做:联合 shape+color 优化、fdef 递归宏整体移植、LumpedHydrophobe、
  ZnBinder、可调权重 UI(后续)
