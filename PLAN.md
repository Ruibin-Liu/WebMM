# Plan: 文档轮——引擎侧能力总览(docs/ENGINE-CAPABILITIES.md)+ README 能力表补全(已完成)

现状盘点结论:README 能力表停在景观栈之前(缺 iMTD-GC/能量景观/NEB 攀岩/振动热力学/简并度);细节埋在 CODE_STATUS 轮次日志;task-contract 与 worker-protocol 是协议文档不是能力参考。本轮纯文档,零代码改动,冻结期合法。

## 交付
1. **docs/ENGINE-CAPABILITIES.md**:按模块的能力参考——每个模块列方法、参数语义、输出含义、精度/验证状态、性能数字、限制与诚实降级语义。模块:解析加氢 / ETKDG / MMFF+GFN-FF / 优化(含 restraints)/ MD(NVE/NVT/GBSA)/ well-tempered metad / iMTD-GC 搜索(调度+CV 池+遗传交叉+两层去重+verify 梯子)/ 能量景观(梯子/disconnectivity/网络图/盆地-转迁表)/ NEB 与攀岩(FIRE/冻结/钉住 CI/回退守卫)/ 振动热力学(FD Hessian/QH 地板/玻尔兹曼/变体几何简并)/ 3D 形状;总则(确定性/诚实原则/冻结面)与文档地图。
2. **README 能力表补行**:iMTD-GC 构象搜索、能量景观、NEB/鞍点精修、振动与热力学四行 + 指向新文档;Performance 节补 NEB/iMTD 数字。
3. 数字全部取自已实测的档(perf_probe / climb_molecules / README 原表),不造数。

## 验收(实测)
docs/ENGINE-CAPABILITIES.md 入库:14 节 + 文档地图(总则/解析/ETKDG/力场/优化/MD/metad/iMTD-GC/景观/NEB 攀岩/vib-thermo/形状/性能速查/限制与降级/验证状态),数字全部取自 perf_probe、climb_molecules、README 原表与各轮实测档(未造数;GFN-FF 拒 NEB 的报错原文核对 lib.rs:1954)。README 能力表 +5 行(iMTD-GC/景观/鞍点精修/振动热力学)+ Performance 节补 iMTD/NEB 数字与新文档指针。零代码 diff;m22 不受影响。

## 排雷存档(本轮)
- 写能力文档前先 ls docs/ + grep README 找到"已有什么/缺什么"——本仓库 README 能力表与实际引擎能力脱节约五个月(景观栈四块全缺),脱节点即文档债。
- 文档数字必须逐个回溯源(perf_probe/climb_molecules/lib.rs 报错原文),不凭记忆写——本轮核对中发现 GFN-FF 拒 NEB 的报错原文、QH 地板常量名(QH_CUTOFF_CM)等细节,凭印象写必错。
