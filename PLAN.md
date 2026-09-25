# Plan: ETKDG L-BFGS 历史 m=8→20 → v1.2.4

## 背景

v1.2.3 把 4D 阶段统一到 lbfgs_minimize(历史 m=8),而主优化器
(optimizer/mod.rs)用 m=20。RDKit 用全内存 BFGS(全部曲率信息)。
在 33 原子(99 维)下,增大历史几乎零成本(每次迭代多 ~12 个
O(dim) 向量运算),可能减少迭代数。

## 任务

1. src/etkdg/mod.rs 两处 `const M: usize = 8` → `20`
   (lbfgs_minimize ~1862、minimize_etkdg ~5059)
2. 门禁:cargo test 275 全绿、ensemble 6/6、benchmark 230/230、
   clippy 0、fmt
3. 交错 A/B(vs v1.2.3):embed1 ibuprofen ms/构象 + 4D 收敛观察
4. 若中性或更好 → 保留 + 发布 1.2.4;若变差 → 回滚并记录

## 验收(实施后实测记录)

- `cargo test` 275/275 全绿;benchmark 230/230 逐位一致(2 个既有
  离群不变);clippy 0;fmt;wasm(node 冒烟 1.2.4);API 零变化
- **交错 A/B(load ~8,3 轮)**:
  - aspirin:24.7/29.0/25.1 → 22.9/21.3/23.5 ms/embed(**~12% 提速**,
    3/3 轮一致)
  - ibuprofen:50.9/54.1/49.9 → 47.0/46.3/50.3 ms/embed(~6%,2/3 轮
    一致,第 3 轮平)
- 结论:小而一致的收益,保留;m=20 与主优化器对齐(一致性)。
