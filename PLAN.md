# Plan: Workbench 全功能审查(重点 Search 标签,真实浏览器端到端)

## 背景

两轮 UI 修复(按钮收编/模式隔离/字体归位)上线后,用户要求全面验证
各项功能"是否真如预想中一样",重点 Search 标签。本轮为审查驱动:
先实测、发现缺陷再修。冻结参考:tests/fixtures/lbdd/search_refs.json
(2 相似查询 × 5 指纹 × 55 库全精度 Tanimoto + 4 子结构命中集)与
rgd_refs.json(对位苯核心 6 命中含片段 SMILES)。

## 审查矩阵(Playwright,localhost:8901,repo root 伺服使引擎可用)

### A. Single
process ibuprofen:canonical SMILES/MW/cLogP/TPSA/HBD/HBA/RotB/
QED mean+max/PAINS/BRENK/Murcko/InChIKey 正确性;History 入史/恢复;
JSME 弹窗;Embed 3D(引擎就绪)→ Optimize MMFF94s 收敛 → GFN-FF 切换
→ Conformers N=10 图表;charge color 与 Features 开关。

### B. Batch
3 分子批跑(描述符+QED/PAINS 列);MW 排序双向;Lipinski 过滤;
CSV 导出 download 事件;行点击回载 single;2 分子 3D 批量 E 列填充。

### C. Search(重点)
1. Demo 库 55 装载;刷新后切 search 自动恢复(localStorage);
   Clear 后不再恢复
2. **相似性奇偶**:aspirin/benzene × morgan/rdkit/maccs/atompair/
   topologicaltorsion,全 55 条 Tanimoto 与冻结参考逐一精确相等;
   UI 路径跑一次(阈值过滤 + 降序 + 自命中居首)
3. **子结构奇偶**:pyridine/phenol/carboxylic/amine_smarts 命中名集合
   与参考精确一致
4. **Shape (3D)**:aspirin 查询自命中 100% 居首;Color/Pharm/Combo
   列出现;Pharm ≥60% 过滤 glucose 出局 salicylic 保留;行点击回载
   对齐构象自动 embed + Features 自动勾
5. **RGD**:Example 核心 Decompose → 6/55,命中名+片段 SMILES 与
   rgd_refs 精确一致;Auto core → [*:1]c1ccc([*:2])cc1
6. **Scaffold**:Analyze 计数和=34(环状分子数);行点击回载骨架
7. **Pharmacophore**:aspirin embed 后 Build → 默认 4 特征勾选 +
   距离矩阵;±1.5/N=1 Screen → 命中恰为 {aspirin, salicylic acid}
   (glucose 负控制蕴含其中);Conf 列=1;行点击回载

### D. 横切
全程零 page error;模式隔离抽查;390px 三模式无溢出。

发现任何缺陷:定位根因 → 修复 → 复测,再收口。

## 验收结果(实施后)

**53/53 通过,零 page error,未发现缺陷,无需修复。**

- A. Single 17/17:ibuprofen canonical SMILES/InChIKey/MW/cLogP/TPSA/
  HBD/HBA/RotB/QED mean+max/Murcko 全对、PAINS/BRENK clean;History
  入史;JSME 弹窗;引擎链全通(embed 33 原子 1.05s → MMFF94s
  25.29 kcal/mol 收敛 211 iters → GFN-FF → Conformers N=10 图表;
  charge color/Features 开关)
- B. Batch 7/7:描述符/QED/PAINS 列(绿色 0=无警徽,初次误报为审计
  断言偏差)、MW 双向排序、Lipinski 过滤、CSV download 事件、行点击
  回载 single(初次误报:乙醇 canonical=CCO 仅 3 字符触发断言下限)、
  3D 批量 E 列填充
- C. Search 24/24:demo 库 55、刷新自动恢复、Clear 清存储;
  **相似性奇偶 550 值全精确**(aspirin/benzene × 5 指纹 × 55,逐一
  === 冻结参考);UI 路径自命中居首/降序/阈值生效;**子结构 4 查询
  命中集全精确**(3/7/10/10);Shape 自命中 100% 居首、Color/Pharm/
  Combo 列、54/55@T≥0.30 用时 1.5s、Pharm ≥60% 过滤后恰为
  {aspirin, salicylic acid}(glucose 出局);行点击回载对齐构象
  自动 embed+Features;RGD 6/55 命中名与片段 SMILES 全精确、
  Auto core=[*:1]c1ccc([*:2])cc1;Scaffold 计数和=34 环状;
  药效团 aspirin 查询 ±1.5/N=1 命中恰 {aspirin, salicylic acid}、
  Conf=1、行点击回载
- D. 390px 三模式无溢出
