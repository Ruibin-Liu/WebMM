# Plan: 颜色力场逐类型权重可调(引擎 v1.6.4→v1.6.5 + app UI)+ ANCopt 评估结案

## 背景

形状检索的颜色层(六类特征)自 v1.5.0 起权重统一 1.0(ROCS 官方权重
闭源,文档化简化)。用户要求可调。ANCopt 为候选清单评估项。

## 引擎(v1.6.5,additive)

1. `ColorSite` 增 `w: f64`(默认 1.0);lib.rs 的 sites JSON 解析接受
   可选 `"w"`(0–3,越界 clamp;缺省 1.0)
2. `color_overlap_grad` 同型对项乘 `w_i·w_j`(梯度按链式同乘——
   线性因子,FH 命名空间不变);`color_tanimoto_at` 的 O_ab/O_aa/O_bb
   走同一函数自动一致
3. **零漂移**:w=1 时 `1.0·e` 与现值逐位一致(乘一不改浮点)——
   现有全部 shape 测试不迁改通过即为证
4. 单测:①默认权重与现输出逐位;②单类型双位点几何,w=√2 → 该类型
   对偶恰为 2×(线性语义锚点);③带权重 FD 梯度;④全零权重 →
   colorT=0(den 0 守卫已在)

## app

5. shape 模式下搜索栏加六输入(donor/acceptor/pos/neg/hydrophobe/
   ring,0–3 步长 0.1 默认 1.0,折叠 details 防 390px 溢出);
   `sitesToEngineJson(sites, weights)` 写 `w = sqrt(u)` 使 **用户权重
   u 对对偶项线性**(w_i·w_j = u);目标侧位点同权重(对称,类型级)
6. m5 增补:①六输入在位、默认 1.0;②全零权重 → 所有行 Color T = 0%
   且 Combo = Shape T;③donor=3 至少改变一行 Color T(与默认对比);
   ④默认下既有 shape 断言全不回归(零漂移端到端)

## ANCopt 评估(决策记录,不实施)

7. 结论写入 CODE_STATUS:动机(鲁棒性)已被 v1.6.4 三层修复+线搜索
   地板关闭;性能动机不成立——v1.1.1 DIC 矩阵实测迭代数 0.44–0.66×
   但墙钟持平或更差(E+G 主导、坐标变换 O(N³) 吃掉收益),ANC 的
   Lindh 模型 Hessian+特征分解属同一成本类;xtb 对齐不可达(不同
   最小是混沌而非算法差);实施成本 = 完整里程碑(模型 Hessian/信赖
   域/收敛判据/重建策略/语料验证)。**重开触发条件**:大规模柔性
   系统的紧优化 API、或迭代数成为瓶颈的证据。备选廉价路径已在库:
   internal_opt(DIC)为 API opt-in,如需可暴露为 app 选项

## 验收

- cargo test 全绿(301+新增)、clippy -D warnings、fmt;wasm 重建 +
  node 冒烟 v1.6.5
- m5 新断言过 + 六套件回归;390px;零 page error;node --check

## 验收结果(实施后)

- 引擎:ColorSite.w(additive),sites JSON 可选 "w"(0-3 clamp),
  对偶项乘 w_i·w_j(梯度线性因子);4 新单测(默认逐位/√u 线性语义
  锚点/带权 FD 梯度/全零→colorT=0);305 测试全绿(既有 shape 14 项
  不迁改通过 = 零漂移)、clippy -D warnings 0、fmt;v1.6.5 wasm
  重建,node 冒烟 1.6.5
- app:shape 模式折叠 "color weights" 六输入(0-3 步长 0.1 默认 1.0,
  reset)+ sitesToEngineJson(sites, weights) 写 w=√u(用户权重对对偶
  项线性),两条路径(单构象/系综)全接
- m5 **36/36**(+4:六输入在位默认 1.0;全零→Color T 全 0% 且
  Combo=Shape T;donor=3 改变 Color T;回默认结果不变);六套件
  全绿(37/10/11/10/32/36);ANCopt 评估结案记录于 CODE_STATUS
