# Plan: v1.0.0 评审修复 — README 默认迭代数同步 + CODE_STATUS 乱码

## 背景

对 v1.0.0 提交(53e05d7)的代码评审发现两处文档问题:

1. **README API 参考段陈旧(该提交自身引入的不一致)**:优化器默认
   `max_iterations` 已在 `src/lib.rs` 改为 1000,但 README 的
   `OptimizationOptions` 示例(L248)与参数默认值表(L328)仍写 `200`。
   陈旧值随构建拷贝传播到未跟踪的 `pkg/README.md` 与 `site/pkg/README.md`。
2. **CODE_STATUS.md 乱码(历史遗留,G6 提交引入)**:L23 "胍!根" 含
   U+FFFD 替换字符,应为 "胍根"(guanidino,胍基)。

纯文档修复:无 Rust/JS 代码改动,无 API 变更,测试数与 clippy 状态不变。

评审发现 3(.gitignore 的 `GENTS.md` 拼写及 PLAN/CODE_STATUS 惰性忽略规则)
明确不在本次范围。

## 任务

1. `README.md` L248 示例 `options.convergence.max_iterations = 200;` → `1000;`
2. `README.md` L328 参数表 `` `convergence.max_iterations` | `200` `` → `` `1000` ``
3. `CODE_STATUS.md` L23 "胍!根"(胍 + U+FFFD + ! + 根)→ "胍根"
4. 同步未跟踪构建副本:`cp README.md pkg/README.md site/pkg/README.md`
   (复刻 wasm-pack 构建 + site staging 的 README 拷贝语义)
5. `CODE_STATUS.md` Recently Completed 顶部按既有 `+- **标题。** …` 格式
   追加本任务条目

## 验收

- `grep -n "max_iterations" README.md`:无默认 200 残留(248/328 两处为 1000)
- `grep -c $'\xef\xbf\xbd' CODE_STATUS.md` = 0
- `git status`:仅 README.md / CODE_STATUS.md / PLAN.md 三个跟踪文件变更
- 不跑全量测试(无代码改动);256 测试 / clippy 0 状态沿用
