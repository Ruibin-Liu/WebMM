# Plan: 双轮 —— Top-4 SAR 决策工具 + Ketcher bundle 扩展(已完成)

## 交付
**Part B:Top-4 SAR 工具**
- sar.js 纯模块:V2000 图解析/小核心回溯匹配(芳香键 4↔1/2 双向容差,
  度剪枝)/附着点检测/R 基团提取(核外连通片→环感知 kekulé SMILES 写出器,
  [*] 附着,先分类树边/回边后发射——祖先开环后代闭合)/透明 MPO(梯形
  desirability 窗口/双向;空值 null 不计入;Lead-like 预设)。
- UI(检索模式「SAR 分析」卡):核心 SMILES/MCS 建议(MolList API,MCS
  只是建议非真相)/R 表(列=附着位,RDKit 再规范化片段+tooltip 保留原始
  [*] SMILES,未匹配行明示不静默丢)/MPO 表(逐准则 原始值(得分)+总分
  +n/N;活性准则可选=所选 assay 最小值,IC50 类越低越好)/CSV 导出。
- 测试:Node +11(124/124);m10_sar 15 断言。

**Part A:Ketcher bundle 扩展**
- ketcher-entry.js 增 api()/getSmarts()/getSelection()(防御式,缺失时
  明确报错);build-ketcher.sh --install 重建 30MB bundle(sanity 过)。
- 运行时探针实证:getSmarts 返回完整 SMARTS(aspirin);getSelection 形状
  {atoms,bonds};editor 有 highlightAtomById(后续可用)。
- 桥接 UI(模态两钮):「SMARTS → 检索」画布含查询特征导出→填查询→切
  sub 模式(画图即查询,m11 断言真实检索命中);「选区 → 3D」molfile 序
  经原子账本核对→3D 高亮+拾取(空选区教学提示不报错)。
- 测试:m11_ketcher 7 断言。

## 验收
Node 124/124;CDP 24 套件(run.sh 含新 m10/m11);零 page error;
Ketcher 既有闭环(mount/getMolfile/setMolecule/apply)零回归。
