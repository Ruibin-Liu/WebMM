# Plan: 柔性叠合(Schrödinger Ligand Alignment 同类,v1.7.0→v1.8.0 = Kabsch/MCS 骨架叠合 + 约束驱动 flex 精修)

## 背景

引擎三资产(ETKDG 构象、v1.5.0 刚性 shape+color 对齐、v1.7.0 几何约束)
恰好拼出 Schrödinger Ligand Alignment 的诚实近似:probe 在叠合中允许
扭转,而非仅刚体。路线:刚性系综检索(既有)→ top-N 胜者做约束驱动
柔性精修。

## 引擎(v1.8.0,additive)

1. **`Restraint::Position{i, target[xyz], fc, tol}`**(POSRES 原语):
   E=½fc·max(0,|x_i−target|−tol)²,梯度只落在原子 i;JSON key
   "position"({i,xyz,fc?,tol?},xyz 必须);FD 单测 + 集成(拉乙醇
   H8 至空间点,tol 内收敛)
2. **`kabsch_align_wasm(mobile_sdf, ref_sdf, pairs_json)`** → JSON
   {rmsd, transform:[12]}(与 applyTransformToSdf 的 m12 同构):
   Horn 四元数法(质心化互协方差 → 4×4 N 矩阵 → sym_jacobi 最大
   特征向量 → R;t = q_c − R·p_c);proper rotation only(无反射,
   文档化);pairs=[[mobile_idx,ref_idx],...] 0 基,<3 对或索引越界
   或 SDF 坏 → Err;单测:已知旋转+平移精确恢复(rmsd~1e-12)、
   噪声扰动下优于恒等、镜像集 rmsd>0、错误路径
3. 版本 1.8.0;既有无约束路径逐位不动(305+9 既有不迁改)

## app(shape 模式 flex 精修,两路径共挂)

4. 搜索栏 shapeFlexWrap(checkbox "flex top 10",title 释义,
   switchMode 与 shapeConfsWrap 同步显隐,默认关)
5. **flexRefineTop(rows, qSdf, cw)**:rows 排序后、渲染前,取
   top-10:
   a. MCS:get_mcs_as_json(MolList[qMol, tMol])(两 3D SDF 均含
      显式 H)→ SMARTS → 双侧 get_substruct_matches 首匹配 zip 成
      pairs;<6 对 → 该行跳过(诚实降级,保留刚体分)
   b. Kabsch 纠正:kabsch_align_wasm(alignedSdf, qSdf, pairs)
      → 纠正姿态(共同骨架入位;matches 恒从已对齐 SDF 出发)
   c. 位置约束:matched 目标原子 ti → 查询原子 qi 坐标(查询系,
      双方同系),fc=25 tol=0.5(文档化常量)→ optimize_from_sdf
      MMFF94s + set_restraints → flexed SDF;非 converged → 跳过
   d. 重打分:shape_align_color_wasm(qSdf, flexedSdf, ..., cw)
      全量 → 更新 score/colorT/aligned,行加 flex 徽标 + title
      (MCS n atoms · post-flex matched RMSD)
   e. 重排序(同 combo 键)+ 状态行追加 "· flex refined top 10"
6. 行渲染:flex 行名旁小徽标(样式既有 chip 先例);conf 单位
   (confs=1 路径 rows 无 conf 字段,flex 不依赖它)
7. 性能预算:top10 ×(MCS ms + Kabsch µs + MMFF 优化 50-300ms +
   全量对齐 44ms)≈ 1-3s,主线程分块 await;确定性:全确定性路径

## 测试

- 引擎:Position FD/集成、Kabsch 四单测;既有 314 不迁改
- node 冒烟:Kabsch 旋转恢复 + position 约束拉引
- m5(+4):flex 开关在位默认关;开启后 top 行有 flex 徽标(≥1 行
  badge,aspirin 自查询 MCS=全分子);自匹配 flex 后仍 100.0/200.0;
  关闭后无徽标且结果与既有逐位(零漂移)
- m0-m4 纯回归;390px;零 page error;node --check

## 验收

- cargo test 全绿(314+新增)、clippy -D warnings 0、fmt;wasm
  重建 1.8.0;六套件全绿

## 验收结果(实施后)

- 引擎(v1.8.0):Restraint::Position(POSRES 原语,平底,E=½fc·
  max(0,|x−target|−tol)²,FD 单测 + 乙醇 H8 拉点集成 0.079<tol);
  kabsch_align(原生)+kabsch_align_wasm(Horn 四元数 × sym_jacobi,
  proper-only;已知旋转 rmsd 1e-15 恢复、噪声优于恒等、镜像 >0、
  错误路径);JSON "position" 接通;319 测试(+6)、clippy 0、fmt
- app:shapeFlex 复选框 + flexRefineTop(top-10:RDKit MCS≥6 原子
  →qmol 双侧匹配→Kabsch 纠正→位置约束 fc25/tol0.5→MMFF94s→全量
  重打分+pharm 重算+flex 徽章/title(MCS n 原子·RMSD))两路径全接;
  **自匹配守卫 score≥0.9995 跳过**(flex 只会让完美叠合受损)
- 实测:aspirin 查询 flex 8/10——自匹配 21 原子 RMSD 0.00、水杨酸
  15 原子 0.32 Å、benzocaine 12 原子 2.74 Å;墙钟 +0.9s;零 page error
- **过程四 bug**:①minimal get_substruct_matches 需 get_qmol
  (SMARTS 串直传抛 "Cannot pass as a Mol");②sdfCoords 是扁平数组
  (qco[qi] 数字→"missing xyz");③badge 模板替换误伤 pharm 表第二
  处 innerHTML(无 flexBadge 声明→ReferenceError→空表)——git stash
  A/B 锁定 app 侧后定位;④m4 cancel 测试 100ms 刀刃(同步描述符相
  ~100-200ms,机器/wasm 时序微移即翻面)——输入加倍 19 行使
  mid-flight 无歧义 + 输入归位保隔离断言,三连稳
- m5 **41/41**(+5:flex 开关默认关/徽章/自匹配不降格 100.0%/关后
  零漂移/pharm 回归);六套件全绿 37/10/11/10/32/41
