# Plan: ECFP4+Tanimoto 引擎导出(v1.11.0)——已完成

## 现状核查(先行事实)

- **应用侧已有 ECFP4+Tanimoto**:search 模式五种指纹之一 "Morgan" =
  MinimalLib `get_morgan_fp`(r2/2048)+ JS `tanimotoBits`;m5 对
  search_refs.json **逐位全等**(55 库 × 5 指纹 × 多查询,Python RDKit
  生成参考)。此前"MinimalLib 无 Morgan"判断是探针函数名用错
  (`get_morgan_fingerprint` ≠ `get_morgan_fp`),予以更正。
- 引擎侧 `src/sascore.rs::morgan_sparse_counts` = bit-exact 展开型
  Morgan r2 标识符计数(golden 77 分子逐位全等)。ECFP4 位集 =
  {identifier % 2048 : identifier 存在}——与 RDKit
  GetMorganFingerprintAsBitVect(mol, 2, 2048) 同一折叠语义。

## 范围

### A. 引擎(Rust,additive)
1. `sascore.rs`:`pub fn ecfp4_fingerprint(g: &SaGraph, n_bits: u32) -> Vec<u32>`
   ——计数键折叠去重升序;`pub fn tanimoto(a: &[u32], b: &[u32]) -> f64`
   ——有序集交/并(全等双除法,与 JS tanimotoBits 同算术)。
2. `lib.rs` 导出:`ecfp4_fingerprint_wasm(molblock) -> Vec<u32>`(固定
   2048;复用 SA 的原始 Kekulé molblock 解析)+ `tanimoto_wasm(a, b) -> f64`。
3. cargo 测试:77 金标分子折叠位集与 golden fp 键集 %2048 全等;
   tanimoto 自身=1.0/对称/空集守卫;若干跨分子值(与金标 fp 在测试内
   折叠重算,不引新夹具)。

### B. 应用(最小接线)
1. index.html 与 platform.html 检索指纹选项 "Morgan" 标签 →
   "ECFP4 (Morgan r2)"(value 键 `morgan` 不动——零行为变更)。
2. m5 增**跨实现对拍**:查询与若干库分子上,引擎
   `ecfp4_fingerprint_wasm(get_mol(smiles).get_molblock())` 位集 ==
   MinimalLib `get_morgan_fp()` 位集,且 `tanimoto_wasm` ==
   `tanimotoBits`(双精度全等)——防 vendor 升级漂移的真回归价值,
   亦为导出的真实消费者。

### C. 版本与构建
- v1.11.0(Cargo.toml);wasm-pack 双重建;pkg 暂存
  (site/index.html + app/fpscores.bin);node 冒烟(导出可用、自相似 1.0)。

## 不做
- 不替换应用检索路径(MinimalLib 已验证且同值,换路径零收益);
  不加属性过滤/PAINS 列(LigandLab 卡片剩余项,另行立项);
  不动探索器预筛(语义变更风险)。

## 验收
1. cargo 全量 + sascore_golden(含新 ecfp4/tanimoto 测试)绿;
2. m0–m6 七套件 + 平台 Node 全绿(m5 +1 跨实现断言);
3. node 冒烟:ecfp4 位集非空、tanimoto 自身=1.0;
4. 零 page error;clippy 0;fmt。


## 实施与验收(完成)

1. **引擎**:`sascore::ecfp4_fingerprint(g, n_bits)`(展开型标识符集
   折叠 mod n,去重升序——RDKit AsBitVect 同语义)+
   `sascore::tanimoto`(有序集双指针交/并;**空并集=0.0 镜像 app
   tanimotoBits/ RDKit 约定**)。lib.rs:抽 `sa_parse_molblock` 共享
   解析(机械重构,SA 路径错误串不变);`ecfp4_fingerprint_wasm`(固定
   2048,原始 Kekulé molblock)+ `tanimoto_wasm`。v1.11.0。
2. **cargo 测试 +2**:77 金标分子折叠位集 vs golden fp 键集 %2048 全等;
   tanimoto 恒等式(空集=0/自身=1/对称/与金标折叠集重算值全等)。
3. **应用**:两页检索选项 "Morgan" → "ECFP4 (Morgan r2)"(value 不动);
   m5 +1 跨实现断言:浏览器内 `ecfp4_fingerprint_wasm(get_molblock())`
   位集 == MinimalLib `get_morgan_fp()` 位集(全部 sim 查询)+
   `tanimoto_wasm` == 集合算术(逐对双精度全等)——vendor 升级漂移的
   真回归护栏,亦为导出的真实消费者。
4. **构建**:wasm-pack v1.11.0 双重建,pkg 暂存;node 冒烟(自 T=1/
   空=0/CCO↔paracetamol 0.0833)。

## 验收数字

cargo 322 + sascore_golden **6**(+2);clippy 0;fmt;m0-m6
37/10/11/10/32/**45**/30(m5 +1 跨实现);平台 Node 34/34;零 page error。
更正记录:此前"MinimalLib 无 Morgan"为探针函数名误用
(get_morgan_fingerprint ≠ get_morgan_fp),应用侧 ECFP4+Tanimoto
检索自始存在且 m5 已逐位对拍。
