# Plan: M5 结构-数据底座(Top-1 展开;Fable-5.1 评审一轮收敛)

## 范围(docs/assay-foundation.md 冻结稿,7 必修全采纳)
- **身份分层**(L2 派生):identity.js 增 isInchiKey/inchiKeyBlocks(块2 可比
  部分=前 8 位)/dupeTier/isLayerless(UHFFFAOY 常量);m5.js 纯函数
  identityAudit(完整键/母体键双分组 + 立体交叉细分 + 降级键不崩);
  QC 旗标(盐/立体已指定/非 InChI 键/母体判定不确定,页面侧 L2 缓存)。
- **实测活性**(L1 命令):SetActivityData(assayId 全量替换,records
  [molId, 数组] 对,qualifier ∈ =,<,>,<=,>=,~,ND,valueRaw 保留,
  matchTier 入档)/RemoveActivityData;撤销走 pre 快照(首导入逆=
  RemoveActivityData);IdentityMerged 并置记录(冲突共存);
  RemoveLibrary 留孤儿(重写允许孤儿 id,幽灵仍拒);ImportLibrary
  保留 rejected 行;快照序列化 activities Maps;旧快照前向兼容默认空。
- **CSV 两步导入**:m5.parseCsv(BOM/引号/;嗅探)+parseActivityCell
  (µ/μ、≥≤、e-记数、欧逗号、ND 族、单位尾缀)+normalizeName
  (NFKC/大小写/空白)+matchLadder(molId→精确键→母体键→块1→名称,
  歧义绝不静默取首);预览报告(匹配/歧义/未解析/松匹配警告)→确认提交。
- **UI**(检索模式):活性面板(导入/列选择/导出带 molId 列/删除/孤儿
  提示)+ 结果表活性列(qualifier 显示,tooltip 全部重复测量+raw+tier)
  + 身份审计卡(家族/旗标/互变覆盖说明);模态注册入 MODALS(空守卫)。
- **测试**:Node +40(reducer 全语义/往返/撤销/快照/重放/纯函数表驱动);
  CDP 新套件 m7_activity 14 断言;run.sh 20 套件。

## 不做
MPO/聚合、立体枚举、互变感知、盐改写、kNN、单位换算(评审记录在案)。

## 验收
Node 113/113;m7_activity 14/14;run.sh 全量 20/20;零 page error。
