# M5 设计:结构-数据底座(身份分层 + 实测活性)—— 冻结稿

> 评审:Fable-5.1 一轮(7 必修/9 建议全部采纳或记录)。
> 本文档为实施契约;「冻结」条款变更需重新评审。

## 0. 定位
补上"看一组分子"的地基:化学身份分层(重复家族/QC 旗标,绝不改写用户结构)
+ 实测活性数据(assay 语义,qualifier 与数值分离,来源可溯)。纯客户端,
库规模几十~几千。不做(本期):MPO、立体枚举、互变感知引擎、盐剥离改写、
kNN、单位换算。

## 1. 身份分层(L2 派生,不进命令日志;渲染时计算 + Map 缓存)
键 = 完整 InChIKey(27 位);注入点校验 `^[A-Z]{14}-[A-Z]{10}-[A-Z]$`,
失败降级键 `SMI:<canonical smiles>` 并打「非 InChI 键」旗标(Node 测试桩
必须产出格式合法假键)。

- **块语义(诚实命名)**:块1=骨架(标准 InChI 仅规范化杂原子 1,3-氢迁移
  类互变;酮-烯醇/环-链不覆盖——面板固定说明行);块2=非主层哈希(立体/
  同位素/电荷,后 2 位为版本标志不参与比较);块3=质子化层(N 中性,M/L/K
  去 1/2/3,O/P/Q 加 1/2/3)。
- **重复分层**:L1 同物(全键同)/L2 块2 异(**立体/同位素/电荷**,交叉
  「立体已指定」旗标细分:双方已指定→立体异构对;单方未指定→指定/未指定
  关系;块2==UHFFFAOYSA 恒=无任何非主层)/L3 质子化异(块1+2 同块3 异)。
- **分组键可切换**:完整键 / 母体键(最大片段)。游离碱 vs 盐酸盐的重复
  家族只在母体键视图出现。母体判定:重原子数最大 + 常见反离子白名单
  (Cl/Br/F/I/Na/K/乙酸盐/富马酸盐/盐酸盐等)+ 平局或反离子不在白名单
  →「母体判定不确定」旗标。**仅展示,绝不改写。**
- **QC 旗标(每行)**:盐/多组分(SMILES 含 '.')、立体已指定(规范
  SMILES 含 @ 或 /\ —— provenance 旗标,非完整性声明)、非 InChI 键、
  母体判定不确定。缓存键 = raw smiles + rdkit 版本串;requestIdleCallback
  分批回填,表格先渲染。
- **被拒行**:ImportLibrary 命令保留 `rejected:[{row, raw}]`(L1,小);
  审计面板可持久列出。

## 2. 实测活性数据(L1 用户输入,进命令日志)
### 数据模型
`state.activities: assayId -> {meta:{name,endpoint,unit,source?,direction?,
importedAt}, records: Map(molId -> Array<{value:number|null, qualifier,
valueRaw, note?, replicateId?, matchTier?}>)}`

- qualifier ∈ `= < > <= >= ~ ND`;ND 时 value=null;valueRaw 保留原始单元格
- **records 数组化**:同 molId 多条 = 重复测量/复孔/批次,常态;渲染
  "8.2 (n=3)";不定义聚合(留 MPO 期);排序键=数组内最小 value,qualifier
  次键;reducer 按 molId 排序 records 保 logHash 稳定
- assayId = UUID(创建时生成,进命令);重导同名 CSV → UI 按 meta.name 匹配
  弹「覆盖现有/另建」;NaN/Infinity 拒绝(canon 层已有,复检)

### 命令
- `SetActivityData {assayId, meta, records}` —— 按 assayId 全量替换(幂等);
  逆命令由 reducer 依 apply 前状态携带先前完整快照(镜像 checkpoint 模式)
- `RemoveActivityData {assayId}` —— 逆=SetActivityData 先前快照

### molId 生命周期(评审必修)
- `IdentityMerged(a→b)`:a 的记录数组**并入** b(拼接,冲突共存);
  逆命令还原(快照含 activities)
- `IdentitySplit`:记录留在原 molId,新 molId 无记录;审计面板标
  「来自拆分,活性未继承」
- `RemoveLibrary`:**留孤儿**(不级联删;删库重导场景活性数据存活),
  审计面板列孤儿记录;「molId 必须存在」仅写入时校验(reducer 拒绝整单并
  返回错误,重放确定性一致),读路径容忍 dangling

### CSV 导入(两步预览)
1. **dry-run 报告**:匹配 N / 歧义 M(一键匹配多个 molId → 拒绝或让选,
   不默认第一个)/ 未解析 K / 同名多行 J(归重复测量)——用户确认后才提交
2. **匹配阶梯**(等级写进 record.matchTier,入命令 L1 —— 派生结果一旦影响
   用户决策即固化):导出 CSV 带 molId 列(最高优先)→ 精确键 → 母体键 →
   块1(后三级需用户确认);name 兜底(NFKC/大小写/空白/连字符规范化,
   同 name 多 molId → 歧义拒绝)
3. 解析器纯函数 + 表驱动测试:>`≥`≤` µ/μ 双码点、`1.2e-6`、欧逗号、
   ND/inactive/N/A、BOM;name 列绝无数值化(Excel 1E3 陷阱)

### 导出
- 面板 CSV(molId/name/smiles/assay/qualifier/value/valueRaw/matchTier/
  note;走既有 csvField 中和)
- 项目导出 v1:命令自动携带;**前向兼容 = 旧应用重放遇未知命令跳过
  (既有行为),活性命令被忽略其余正常** —— 文档化,不升版本号

## 3. UI(检索模式内,不新增模式)
- **活性区**(检索表下方卡):导入(CSV 粘贴/文件,两步)、assay 列管理
  (默认显示最近导入,列头带单位)、孤儿记录列表、assay 删除入口、面板导出
- **检索表**:选定 assay 追加活性列("≥ 前缀 + value",tooltip 显示全部
  重复测量与 valueRaw)
- **身份审计卡**(检索模式):重复家族(分组键切换 完整/母体键;家族内
  标 L1/L2 细分/L3;点击行跳转)、QC 旗标汇总、被拒行、拆分未继承标记、
  InChI 互变覆盖范围固定说明行

## 4. 测试(冻结)
- Node:reducer(新命令/逆/重放/Merged 并入/RemoveLibrary 孤儿/撤销还原/
  未知 molId 整单拒/NaN 拒/records 排序稳定 logHash);identity 纯函数
  (块解析/27 位校验/降级键/UHFFFAOYSA/母体判定含反离子白名单);CSV
  解析器表驱动(≥/µ/e-notation/欧逗号/ND/BOM/name 数值化陷阱)
- CDP 新套件 m7_activity(套件序号 m7 = 本 M5 里程碑,头注释说明):
  种子库含已知重复(盐/立体/质子化家族)→ 审计面板家族断言;CSV 导入
  dry-run 报告 → 提交 → 列渲染 → 撤销 → 重载持久;歧义匹配路径;
  assay 删除;孤儿(删库后活性面板仍列出)
