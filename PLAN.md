# Plan: 仓库清理第二轮 — stash/诊断例程/一次性脚本/计划归档全清

## 范围

用户批准第一轮报告中的全部保留项清理。另含同类别 decisive 项:旧对拍
bin 两件套(其唯一驱动脚本 compare_rdkit.py 一并清除)。git 历史均可恢复。

## 依赖核查结论(不可删)

- `scripts/benchmark_mmff.py`(门禁)依赖 examples `bench_mmff`、`dump_types_energy` → 保留
- RDKit 2026.03 参考重生成依赖 `parity_2026` → 保留
- CI(rust.yml)仅跑 fmt/check/clippy/test,不引用任何待删文件

## 任务

1. `git stash clear`(3 个 v0.5.0 时代 stash,用户确认不可恢复操作)
2. `git rm` examples ×10:conf_parity, diag_angles, diag_embed, gff_audit,
   gff_bdump, gff_metals, gff_qdump, gff_rep, mmff_charged, test_nh2
3. `git rm` 旧对拍两件套:src/bin/compare_mmff.rs, src/bin/compare_etkdg.rs,
   scripts/compare_rdkit.py(已被 benchmark_mmff.py 与
   gen_etkdg_ref.py+validate_etkdg.py 取代)
4. `git rm` scripts ×5:diag_angle_sb, diag_compare, diag_mmff_divergence,
   diff_atom_types, diff_detail
5. `git rm` docs/plans/ ×5(2026-03-21/03-22/03-28/04-19/09-14,均已完成)
6. `rm` pkg/ 未跟踪杂项 ×3:index_old.html, caff_check.sdf, server.py
7. CODE_STATUS.md Recently Completed 顶部追加本任务条目

## 验收

- `cargo test` 256 全绿(例程减少,测试数不变)、`cargo clippy --all-targets` 0 警告
- `git stash list` 为空
- grep 无待删文件的活引用(README/CI/scripts 门禁)
- 提交后 git status 干净
