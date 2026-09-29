# Plan: 形状预制缓存 —— 项枚举/自体积一次性化(v1.6.2)

## 背景与无损性论证

v1.6.1 成本模型 T ≈ 65(起点)+ 50(自体积)+ 8×25(重打分)+ 25(polish)。
其中重复计算有二:
1. **每次重打分**(overlap_full)都对查询与目标**重新枚举项规格**——
   而项剪枝判据 v/(V_e+V_j−v) ≥ EPS 只依赖原子间距离
   (cross = ½(αΣα|c|²−|Σαc|²) 平移/旋转不变)→ **刚性姿态下项结构恒同,
   枚举一次即够,位值不变(无损)**
2. **自体积 vq/vt**(位形不变量)逐次调用重算;库扫描中查询侧重复 N 次

## 设计(引擎,无 API 变化)

1. `PreparedShape { specs, atoms, self_overlap(惰性) }`:
   `prepare_shape(&atoms)`;`overlap_prepared(pa, a_atoms, pb, b_atoms)`
   (materialize 于当前姿态 + 双线性);`overlap_full` 变薄壳(测试/screen 兼容)
2. align_colored:查询/目标各 prepare 一次,重打分循环只 materialize+双线性;
   自体积惰性一次
3. **wasm 层透明缓存**:`Mutex<HashMap<u64, Arc<PreparedShape>>>`(SDF 字符串
   hash 为键,容量 512 超限清空)——库扫描中查询侧全缓存、重复检索目标侧
   复用;shape_align_wasm / shape_align_color_wasm 签名不变
4. 版本 1.6.1 → 1.6.2(纯性能,值位不变)

## 验收

1. **位值无损**:单测 prepare 路径 vs 旧 overlap_full 位值一致(同姿态
   bit-exact);全部既有测试(自对齐 1.000000、等价性 <0.02 等)不回归
2. 性能:node bench rescore_top=8 重复调用(查询缓存命中)398ms →
   预期 ≤300ms;库扫描逐条目耗时再降(实测报告)
3. 质量回归:浏览器两段式 top-20 召回 19/20 不变;sim/sub/390px/
   零 page error;294 测试 + 新增;clippy(1.98,-D warnings)、fmt
