# Plan: 检索属性过滤(MW/cLogP 区间)+ PAINS 标记列 —— 已完成

## 范围(纯应用层,两页同构)

1. 检索控制行新增:MW min–max / cLogP min–max 数字输入 + "标 PAINS"
   复选 + filterStatus 提示行——全部 view-layer,经
   renderSearchResultsCurrent 重渲染即时生效,清空即恢复。
2. cLogP 惰性:CrippenClogp(MinimalLib get_descriptors)按条目缓存
   (e.clogp);MW 用库条目自带值。
3. PAINS 标记:lbdd_data A/B/C 目录(~480 SMARTS,lbddQueries 既有)
   惰性匹配(e.pains={n,first} 缓存),仅勾选时计算,硬顶前 200 显示行;
   行内 fail 徽章(flex 徽章同款),filterStatus 报 flagged 计数。
4. 不加表格列(m5/m6 断言依赖 children 索引);不做 PAINS 过滤语义
   (LigandLab 是"标记"非"筛除")。

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


## 实施与验收(完成;检索属性过滤 + PAINS)

1. **两页同构五处补丁**(controls/helpers/filter block/badge/status)
   一次锚定落盘;7 内联脚本 node --check 通过。
2. **实测**:T=0 全库 55 行,MW≥200 → 8 行("47 filtered by MW/cLogP"),
   叠加 cLogP≤2 → 0,清空 → 55 全恢复;PAINS 探针 C#CC(=O)C#C
   (pentadiyn-3-one,命中 ene_one_yne_A(1),无显式氢依赖)入库检索,
   勾选后该行 PAINS 徽章 + filterStatus "1 PAINS-flagged"。
3. **排雷两条**:①SMARTS `-[#1]` 只匹显式氢原子——CC=S(乙硫醛)不中
   thio_aldehyd_A(Python RDKit 证实;与 batch 列行为一致),探针改用
   无 [#1] 模式的 pentadiyn-3-one;②loadSearchLibraryFrom 收对象数组
   [{smiles,name}] 非文本(传文本静默得 0 库)。
4. **m5 +2**(过滤收窄+恢复;PAINS 徽章)且**段落自洁**:PAINS 探针
   后恢复 55 库(状态污染教训:下游 shape 断言硬编码 screen 55×10)。
5. **环境维护**:playwright chromium 缓存 1234→1243(旧二进制被清),
   七套件 + 基线脚本路径按 README 惯例同步升级。

## 验收数字

m0-m6 37/10/11/10/32/**47**/30(m5 +2);平台 Node 34/34;两页 7×2
脚本 check;390px 溢出 0;控件目检(两行 wrap 无重叠);零 page error。
