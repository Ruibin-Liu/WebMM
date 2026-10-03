# Plan: 投影位点法(v1.8.0→v1.9.0,Color T 打分语义变更,立项记录)

## 背景

调研结论(知识库+aspirin 实证):ROCS 隐氢 donor 投影至缺氢位置 ~1 Å、
LigandScout donor@H/acceptor@孤对;Catalyst 的 2.4-3 Å 搭档位置是口袋
查询口径不适用于叠合打分。我们采纳短投影:donor=氢实位(显氢构象零
估计)、acceptor=虚拟孤对 1.0 Å(羰基反轴/角平分线反向);hydrophobe/
ring 不变。解决两个问题:色重叠无方向敏感性 + Features 球无指向。

## 引擎(v1.9.0,additive)

1. lib.rs parse_sites:原子位点条目接受可选 `"off":[dx,dy,dz]`,
   c = atom.c + off(α 仍取锚原子——保守,自重叠量级与现状可比);
   ring 条目不变;无 off 逐位不动。单测:off 精确生效 + 缺省零漂移

## app(colorSites/sitePos/sitesToEngineJson/消费方)

2. colorSites:邻接解析增键级;donor/acceptor 位点带 `proj` 描述
   (逐构象可重算,不烘焙绝对坐标):
   - donor:`{h:[H 索引]}` → p = 氢质心(多 H);空则回退原子
   - acceptor 羰基(C=O×N):`{ext:cIdx}` → p = i + 1.0·unit(i−c)
   - acceptor 角平分线(≥2 邻居):`{nb:[j,k]}` → p = i − 1.0·unit(ûj+ûk)
   - acceptor 单邻居(腈 N):`{ext:nbr}`
   hydrophobe/pos/neg 无 proj(原子);ring 质心不变
3. sitePos(s, coords, transform):proj 分支先算 p 再过 transform;
   pharmMatch/药效团编辑器/Features 可视化自动获得投影位置
4. sitesToEngineJson(sites, weights, coords):带 proj 的原子位点
   发 `{i,t,w,off}`(off = p − 锚原子坐标,1e-6);查询/目标/flex
   三处调用点接 coords(qCoordsCache/sdfCoords(esdf)/sdfCoords(flexed))
5. Features 复选框 title 文案更新(投影语义);零 UI 结构变化

## 语义变更与参考(诚实清单)

6. Color T 全体数值漂移(位点移 ~1 Å):m5 中**钉死绝对色值的断言**
   按"新值即化学正确值"更新并在提交说明逐条列出理由;检索排序可能
   微动,人工过 demo 库 top-10;药效团 ±1.5 Å 命中集可能变化(约束
   距离现在量在相互作用位置)——如实记录新命中集
7. 自匹配 colorT 恒 1.0(结构性保证);色权重(可调)与投影正交叠加

## 验证实验(实施时跑,结果记入验收)

- aspirin 自对齐 colorT = 1.0
- aspirin→salicylic 投影前后 colorT 对照
- 方向性对照:同受体数不同孤对向的分子对,投影后 colorT 应降
- demo 库 top-10 前后排序对照

## 验收

- cargo test 全绿(+off 测试)、clippy 0、fmt;wasm 重建 1.9.0 双冒烟
- m5:新断言(donor 全带 h 投影/acceptor 全带 proj/sitesToEngineJson
  发 off)+ 漂移值更新;六套件全绿;390px;零 page error

## 验收结果(实施后)

- 引擎 v1.9.0:parse_sites 可选 "off"(原子位点 c=atom.c+off,α 仍取
  锚原子;ring 不变;无 off 逐位零漂移);320 测试(+1 内核)、
  clippy 0、fmt;node 冒烟:off 双侧自对齐仍 1.000000、单侧 off 交叉
  0.4493<1.0 证解析生效
- app:colorSites 邻接+键级解析提升至函数首(TDZ 教训),donor/acceptor
  位点带 proj{h|ext|nb}(逐构象重算,不烘焙绝对坐标);sitePos proj
  分支(pharmMatch/编辑器/Features 自动跟随);sitesToEngineJson 五
  调用点全接 coords,投影位点发 off(1e-6)
- **验证实验(aspirin 查询,原子中心 → 投影)**:自匹配 100.0→100.0
  (结构性保证);salicylic 55.6→39.8、benzocaine 32.6→25.5、
  paracetamol 16.5→11.6(方向不一致受罚=预期语义)、ibuprofen
  15.2→15.6(方向一致不受罚);与 flex 精修协同(可扭转对齐方向)
- m5 **44/44**(+3:donor 全@H、acceptor 全带投影、off ~1 Å 且
  hydrophobe 保持原子);六套件全绿 37/10/11/10/32/44;390px;
  零 page error;Features 球现为有指向显示(donor@H/acceptor@孤对)
