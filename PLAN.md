# Plan: 仓库清理 — 过时文档/垃圾文件/.gitignore 卫生

## 范围

评审后续清理。只处理"自证过时 / 纯垃圾 / 失效规则"三类;其余候选项列入
"明确不动"待用户决策。无 Rust/JS 代码改动,无 API 变更。

## 调查结论(清理依据)

| 项 | 依据 |
|---|---|
| `.DS_Store` | macOS 垃圾文件,2026-07-24 误入库且 tracked |
| `PROJECT_STATUS.md` | 文件头自证 outdated("Please see CODE_STATUS.md"),数据停留在 165 测试时代;README/site/app 无引用 |
| `ETKDG_MMFF_REVIEW.md` | RDKit 2025.09 时代一次性审计;CODE_STATUS 历史条目已判其 "unreliable"(标注 FIXED 实则未修);仅历史条目与旧计划引用;git 历史保留 |
| `.gitignore` 失效行 | `GENTS.md` 为拼写错(匹配不到任何文件);`PLAN.md`/`CODE_STATUS.md` 规则惰性(两文件均 tracked) |
| `CODE_STATUS.md` Current Focus | 仍为 de-rotation/Playground v1.1.1 时代描述,v1.0.0 已发布,需刷新(模板各节不动) |

## 任务

1. `git rm .DS_Store`;`.gitignore` 增加 `.DS_Store`
2. `git rm PROJECT_STATUS.md`
3. `git rm ETKDG_MMFF_REVIEW.md`
4. `.gitignore`:删除 `GENTS.md` 错拼行与 `PLAN.md`/`CODE_STATUS.md` 惰性规则
5. `CODE_STATUS.md` Current Focus 刷新为 v1.0.0 后状态;Recently Completed
   顶部按既有 `+- **标题。** …` 格式追加本任务条目

## 明确不动(待用户决策)

- 3 个 v0.5.0 时代 stash(破坏性操作)
- examples/ 下 10 个无引用诊断例程(conf_parity/diag_angles/diag_embed/
  gff_audit/gff_bdump/gff_metals/gff_qdump/gff_rep/mmff_charged/test_nh2)
  —— AGENTS.md 约定的诊断机制,删除属范围决策
- scripts/ 下 diag_*/diff_* 等一次性诊断脚本
- docs/plans/*.md(刻意的 dated 计划归档)
- pkg/ 内未跟踪本地杂项(index_old.html/caff_check.sdf/server.py,
  可能是本地 dev 文件)

## 验收

- `cargo test` 256 全绿、`cargo clippy --all-targets` 0 警告(验证不回归)
- grep 确认无 PROJECT_STATUS/ETKDG_MMFF_REVIEW 活引用
- git status 仅剩本任务变更,提交后干净
