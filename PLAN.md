# Plan: 冻结期两件套——接口稳定化(golden 钉面)+ 测试加固(覆盖缺口与破坏性自验)(已完成)

用户指令「稳定接口」+「加强测试」。接口 golden 守卫本身即最强测试;同时补齐本轮功能演进暴露的覆盖缺口,并要求新守卫**能抓人为破坏**(自验红)。

## 交付
1. **接口稳定化**(原方案):
   - 版本 1.11.0 → 1.12.0(Cargo.toml + package.json;冻结线标记);
   - **m22_api_surface(纯 Node)**:pkg/webmm.d.ts 导出函数全集+各类方法全集 vs golden(精确相等);三 worker `msg.type` 集合 vs golden;contract.js KINDS/ENVELOPE_VERSION vs golden;window.IMTD 键集 vs golden;版本一致性;
   - docs/STABILITY.md(稳定域/semver 纪律/破坏性流程)+ docs/worker-protocol.md(三 worker 消息逐条)。
2. **测试加固——引擎(原生)**:
   - **系综容量逐出路径**(高风险未测):max_conformers=3 强制逐出——五统计向量+变体+ID+pair 重排一致性、landscape 对齐、refined 按 ID 存活;
   - GFN-FF 热力学路径 + 精确坐标缓存正确性(同坐标两次 thermo 逐位相等;缓存命中不改结果);
   - boltzmann 短 degeneracy 向量(chain+repeat 补 1 路径);
   - JSON 编解码恶意输入(坏浮点串/错形 xyz/缺字段)诚实报错不 panic;
   - NEB 退化输入守卫(端点原子数不等→debug_assert;同端点已测);
   - 性质测试(proptest 增补):boltzmann(和=1、能量单调、g 乘子单调)、kabsch(平移/旋转不变式)。
3. **测试加固——CDP(m23_hardening)**:
   - conformers kind 的 ↻ 重跑(m21 只测了 imtd);
   - **重载遗弃**:verify 进行中刷新页面 → boot 对账置 abandoned(阶段记录进薄壳语义);
   - verify 中途取消(扣消息法确定性触发)→ 记录 failed + 原因;
   - 恶意 minima JSON 经 UI 路径被拒(gather→false 不崩);
   - 较大分子端到端(ibuprofen SMILES)imtd+验证+精修冒烟——性能包络不退化。
4. **破坏性自验**:开发中临时破坏(m22 删一条 golden、引擎改一处行为)确认红灯,再还原——守卫的守卫。

## 验收(实测)
m22 24/24(含破坏性自验:删 golden 条目红、篡改 worker 消息红,还原后绿);m23 10/10;cargo test 369(lib,+7:逐出一致性/GFNFF 热力学缓存逐位/短筒并向量/恶意 JSON 拒绝三则 + proptest 三性质);run.sh 36/37 CDP(m1 xtb 环境项)+ Node 153/153;fmt/clippy 双 0;版本 1.12.0 三处对齐;docs/STABILITY.md + docs/worker-protocol.md 落档。**修复一处真实取消缺陷(见排雷 4)。**

## 排雷存档(本轮)
- **golden 生成必须从现实提取**:手写 golden 与真实 d.ts 面差 20 处(漏 ecfp4/tanimoto/sa_* 等旧导出、多写幻觉方法、漏 `msg.type !== 'run'` 负向守卫)——守卫先抓的是我的猜测而非代码;正道 = 从构建产物提取 → 固化 → 精确相等。wasm-bindgen 内部符号(init/initSync/free/toJSON…)须在提取侧过滤。
- **proptest 语义边界即浮点边界**:「g 乘子严格抬升布居」在 f64 下有饱和反例(1−1e-35 == 1.0)与下溢反例(4·e^{−494}==0)——性质断言改为非严格 + 饱和分支,数学上正确的形式。
- **JSON 字面量 null 是合法载荷**:refine_result_from_json("null")=Ok(None) 是"无鞍点"的诚实语义,不是畸形输入;恶意测试要打 "{"、坏浮点串,别打 null。
- **真实缺陷(测试加固抓出):换分子不取消在途 verify/refine**——reset3DState 只清 handle 引用,cancel 闩未扳,Promise.race 永等(令牌泄漏 + 记录悬空);修复:clear 前先 handle.cancel('molecule changed'),驱动 verify 增 cancelled 分支(记录 cancelled+原因,状态行诚实)。这正是 m23 扣消息法要抓的类目。
- m22 为纯 Node 套件(无浏览器):解析 d.ts/worker 源/contract.js/imtd.js,秒级运行;run.sh 无需服务器但共存无害。
