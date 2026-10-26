# Plan: BeautifulRust lint 分级接入(clippy.toml 激活)

## 范围

1. **Cargo.toml [lints] 分级集合**(用户已放置上游 clippy.toml,无 [lints] 则惰性):
   启用 —— rust: unsafe_op_in_unsafe_fn/rust_2018_idioms/unreachable_pub;
   clippy: all + unwrap_used/expect_used(allow)/panic_in_result_fn/
   dbg_macro/todo/string_add/wildcard_imports/redundant_clone。
   延期(注释记录理由):pedantic、as_conversions/cast_*(数值内核 2000+ 条)、
   print_stdout/stderr(examples/benches 是诊断 CLI)、missing_docs 族
   (584 项,内部 crate)、allow_attributes*(41 处存量,随改随转)。
   实测全量 4810 警告 → 分级后清零可行。
2. **src 修到零**:29 处 unwrap→expect(带理由:gfnff 内嵌参数编译期不变式/
   graph BFS 树性质/shape colored 蕴含/etkdg 环原子必键/internal_opt 缓存成
   对设置);8 处 redundant_clone 删除;unreachable_pub 降级(instr 计时模块
   →pub(super),sascore→pub(crate));6 处隐藏生命周期补 <'_>;
   sascore 两个 _public 测试专用访问器 allow(dead_code) 带理由(仅金测
   直含文件使用)。
3. **诊断文件豁免**:benches/examples/tests/sascore_golden 文件级
   #![allow(clippy::unwrap_used)] 带理由(诊断工具,panic 即报告)。
4. clippy.toml(上游逐字节一致)随 [lints] 入库成套。

## 验收
- clippy --all-targets 与 --target wasm32 双 0 警告;fmt clean;
  cargo test 328+7 全绿;WASM 重编 + CDP 全量 19/19。
