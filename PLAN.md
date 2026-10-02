# Plan: 优化器几何约束(v1.6.5→v1.7.0,引擎 + WASM additive,app UI 另立项)

## 背景

引擎收尾清单 #3:`OptimizationOptions` 目前只有 engine/coordinates/
convergence,无任何 restraint——"配体准备(docking prep)"的真实需求是
保持药效团部分刚性、只松弛其余,或锁定特定二面角/距离优化侧链。

## 引擎设计

1. **新模块 `src/optimizer/restraints.rs`**
   - `Restraint` 枚举:Distance{i,j,r0 Å,k,tol}、Angle{i,j,k,a0°,k,tol}、
     Dihedral{i,j,k,l,a0°,k,tol}(k:kcal/mol/Å² 或 /rad²,默认 10;tol:
     平底半宽,默认 0——纯谐波)
   - `RestraintSet{restraints, frozen:Vec<usize>}`;
     `energy_and_gradient(coords)->(f64, Vec<[f64;3]>)`
   - 平底语义:|x−x0|≤tol→0,否则 ½k(x−x0∓tol)²;二面角差先包到
     (−180°,180°](a0 近 ±180 时正确)
   - 梯度:距离 = r_hat 链式(metad DistanceCV 同式);二面角 = 复用
     `etkdg::dihedral_gradient_contrib(coords,i,j,k,l,dE/dφ)`(metad 先例);
     角度 = cos 链式新推导(|sinθ|<1e-8 → 跳过贡献,退化守卫文档化);
     k=0 → 项跳过
2. **`optimizer/mod.rs`**:`optimize_with_restraints(ff,coords,conv,
   Option<&RestraintSet>)`——None/空 → 原 `optimize()` 逐位不动;否则
   `RestrainedObjective` 装饰 CartesianObjective:f_and_g 加约束 E/g 后
   冻结行清零,自身 force_stats(FF+约束、排除冻结——受约束面的物理
   判据);energy() = 同一 energy_and_gradient().0(线搜索路径一致性
   同 shape color 先例)
3. **lib.rs**:`OptimizationOptions` 增 `restraints_json: String`
   (skip)+ `set_restraints(json)` setter(存原文,延迟解析);两引擎
   dispatch(MMFF94/94s、GFNFF)统一接;解析错/索引越界/internal 坐标
   组合 → `message:"Restraints error: ..."`(Parse error 先例);
   JSON schema:{"freeze":[0-based 索引],"distance":[{i,j,r0,k?,tol?}],
   "angle":[{i,j,k,a0,k?,tol?}],"dihedral":[{i,j,k,l,a0,k?,tol?}]}
4. 版本 1.7.0(新能力,additive;不设约束的调用逐位不变)

## 测试

- 单元:三类 FD 梯度(平底内外两侧、二面角 ±180 包裹、k=0 跳过、
  共线退化守卫)
- 集成:①冻结 aspirin 前四原子 → 冻结坐标逐位不变、能量有限;
  ②二面角约束驱动乙醇/丁胺 φ→a0(|φ_final−a0|<5°);③距离约束
  |d_final−r0|<0.05 Å;④错误路径(坏 JSON/越界索引/internal 组合);
  ⑤无约束全路径逐位(既有 305 项不迁改 = 零漂移)
- node 冒烟:set_restraints → optimize_from_sdf,冻结/二面角各一例,
  v1.7.0

## 验收

- cargo test 全绿(305+新增)、clippy -D warnings 0、fmt;wasm 重建
  1.7.0;六套件 m0–m5 纯回归全绿(零 app 改动)

## 验收结果(实施后)

- 引擎:src/optimizer/restraints.rs(Restraint 三类 + RestraintSet,
  平底谐波,二面角 ±180° 包裹;距离=DistanceCV 同式、二面角=复用
  etkdg::dihedral_gradient_contrib、角度=cos 链式新推导 + 共线守卫);
  optimize_with_restraints(None/空 → 原 optimize() 逐位);Restrained
  Objective 加 E/g 后冻结行清零、force_stats 为受约束面口径;
  energy() 同一携带函数(线搜索一致性)
- lib.rs:OptimizationOptions.restraints_json + set_restraints(存原文
  延迟解析,"k" 为 fc 别名);两入口解析 + "Restraints error: ..."
  早退(坏 JSON/越界/internal 组合);两引擎 dispatch 接通
- 测试:**314 全绿**(+9:FD 三类×平底内外/±180 包裹最短弧/k=0/
  共线守卫/冻结逐位/二面角 0.00°/距离 fc=1000 平衡/four 错误路径);
  clippy -D warnings 0、fmt;v1.7.0 wasm 重建,node 冒烟(二面角
  0.00°、越界索引报错);六套件纯回归全绿 37/10/11/10/32/36
- 实施修正:距离测试初版 fc=100 得 2.22 Å 为正确力平衡(MMFF 键合
  力顶住)非 bug——fc=1000 平衡点 <0.1 Å,期望收紧后通过
