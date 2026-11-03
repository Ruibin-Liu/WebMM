# Plan: 能量景观视图接入 workbench(引擎轨迹记录 + 断连图 + 盆地统计)(已完成)

用户裁决(设计页三决策全按建议):①景观数据**始终记录**(契约零变更,UI 按需展示);②v1 主视觉=**断连图**;③纯 MD 盆地模式**不排,记 v2**。远期项就此进入实施。

## 交付
1. **引擎轨迹记录(src/metad/imtd.rs + mod.rs 一个观测方法)**
   - `MetaDynamics::underlying_energy(&self, coords)`(无偏势能观测;零语义改动,不触山体状态)。
   - `Minimum` 增字段:`pe_max_since_prev`(同 run 相邻淬灭帧窗口内**无偏** PE 峰=观测转移势垒上界)、`snapshot`(run 内淬灭序号);execute_run 每步采无偏 PE 维护窗口峰(重置于快照点);minima JSON 增 pe/snap 字段(解析容缺省)。
2. **折叠统计与系综对齐**
   - Ensemble 提取两级同构象判据为 `find_conformer(&coords,&energy)->Option<usize>`;`try_add_quenched` 返回改 **Option<usize>**(成员位置;创建即计 1 次访问);容量逐出同步收缩统计向量。
   - ImtdSearch.fold:逐帧分类(重复→访问+1;新→try_add;窗/容量拒→pruned);同 run 相邻**已指派**帧:异盆→转移计数+势垒取历史最小上界,同盆→驻留+1;跨 run 不计转移;pruned 帧断链(记录语义)。
   - finalize 重构出 `finalized_members()->(Vec<Conformer>, 旧→新映射)`:result() 与 landscape() 共用,统计按最终序重排(窗滤/容量截断丢弃的成员统计随之丢弃)。
3. **断连图(确定性合并树)**
   - `DgNode{threshold, basin:Option, children}`;构造=Kruskal 升序(势垒,from,to)+union-find,每次合并生成内节点;观测不连通⇒森林(诚实)。叶 threshold=盆地能量。
   - 类型:`BasinStats{visits,residence,transitions_in/out}`、`TransitionStats{from,to,count,barrier}`(有向计数,势垒无向取最小)、`Landscape{basins,transitions,frames,pruned}`。
4. **JSON/WASM(加法)**
   - `landscape_to_json`(能量/势垒/阈值走 f64 字符串通道,同 spec/minima 纪律);`ImtdSearch::landscape()`(原生,测试面);wasm `ImtdSearch.landscape_json()`(d.ts +1 方法);驱动 result 增 `landscapeJson`(gather 后于协调器上取)。
5. **UI(index.html,confPanel 景观区)**
   - `landscapeWrap`:诚实横幅(metad 偏置采样,占据定性;势垒=观测上界非鞍点)+ `dgraphSvg`(手绘 SVG:叶=盆地按能量排布 y=E、内节点 y=阈值 x=子均值;叶可点/键盘 ←→)+ 盆地表(E/占据%/驻留/转入转出,行点击联动)。
   - 叶子点击→`showEnsembleConf(i)`(复用 3D/能量面板/状态行);iMTD 完成时 `lastLandscape=JSON.parse(res.landscapeJson)` 并绘制;ETKDG conformers 路径(finishConformers)与 reset3DState 隐藏景观区。
6. **测试**:原生(轨迹字段有限/序号连续;frames 守恒 Σvisits+pruned≤frames;全局 Σin==Σout;势垒≥max(E_a,E_b)-ε;断连树叶数==系综数;顺序==散射-收集 landscape 逐位;JSON 往返);CDP 新套件 **m17_landscape**(乙醇:景观区可见、叶数==系综数、占据列和≈100%、横幅文案在、点叶→selConf 联动、表行数、零 page error);run.sh 挂 m17。
7. **不做(记录)**:纯 MD 盆地模式与力导向网络图(v2)、鞍点精搜、占据的玻尔兹曼重构。

## 验收(实测)
cargo fmt/clippy 双目标 0(sort_by→sort_by_key 一处修复);cargo test 347(lib,+4 景观)+ 7(sascore_golden)全绿;wasm 重编 d.ts 导出集合零变化(ImtdSearch 类 +1 方法 landscape_json);run.sh 29/30 CDP 全绿(m1 xtb 环境项)+ m17 11/11 + Node 153/153。

## 排雷存档(本轮)
- **势垒必须采无偏势能**:runner.potential_energy() 含 metad 偏置(山体填井后峰被压平)——给 MetaDynamics 加零语义观测方法 underlying_energy(&self)(走内层 ff.energy),每步一额外能量求值;execute_run 持 Rc 克隆保句柄(原实现 move 进 MDRunner 后无法再采样)。
- **逐出路径的位置重排**:容量逐出收缩五条平行统计向量外,pair_stats 里的盆地号也要 >i_max 减一——首版占位写法被编译器抓住,重写为 shift 闭包;返回位置=pos-1(新成员永不被逐,cap 检查保证 strictly-below-max)。
- **转移链语义**:pruned 帧(窗/容量拒)断链——相邻已指派帧之间才计转移,跨 run 不计;驻留=同盆连续快照对数。frames 守恒为不等式(Σvisits+pruned≤frames,因 finalize 窗滤/截断丢成员丢访问)。
- **锚点差异再犯**:finishConformers 的实际代码序(showSdfIn3D 在 display=block 前)与记忆不符——脚本断言失败即退出不写盘,重读真实块后 5 锚全中;同轮 m17 的 ETKDG 隐藏路径用真实断言兜住。
- **JS 闭包递归**:SVG 叶计数首版用闭包 let count_leaves = |n| {... map(count_leaves)}(Rust 闭包不能自引用)——改 fn 声明;浏览器端 drawLandscape 的树遍历无此问题(JS 闭包可递归)。
- dgraph 森林语义:观测不连通的盆地各为一根(诚实"observed disconnected"),非断言真不连通;Kruskal 序 (barrier, from, to) 确定性。
