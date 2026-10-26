# Plan: Rust Review 应修项修复(4 提交,每步带回归)

## 范围(承接本轮 rust-review)

1. **parser 加固**(发现 1/2/3,单提交):
   - parser.rs:130 守卫 `len < 7` → `len < 9`(键行 [6..9] 切片);
   - 六处字节切片(parser.rs:80,116,123,145,417,438)统一 `.get(a..b)`
     + 带行内容的错误信息(UTF-8 边界 + 越界一族收敛);
   - 键行 atom1/atom2 上界校验:越界 → `Err`(替代 graph.rs:15 的静默
     丢弃;守卫旧键不动,graph 侧守卫随后续兼容清理);
   - 回归:8 字符键行 → Err;越界键 → Err;多字节字符跨列 → Err;
     合法输入全量测试不回归。
2. **graph 防御**(发现 5):环感知 BFS 父链 `parent[v].unwrap()`
   (graph.rs:524/535/541 一带,以实际代码为准)改为 None 安全路径;
   回归:含孤立原子的合法分子(如水合离子)过 `MMFFForceField::new`
   不 panic。
3. **NaN 排序 total_cmp 族**(发现 7,八处):shape/mod.rs:882,894,
   1260;optimizer/internals.rs:204;etkdg/mod.rs:1805;gfnff/mod.rs:
   758,1617,3348 —— `partial_cmp().unwrap()` → `total_cmp`;
   回归:含 NaN 输入的模块内单测不 panic。
4. **sascore 溢出**(发现 6):sa_load_table 的 `4 + n*8` 改 u64 精确
   算术,截断校验不可绕过;回归:声明 n=0xFFFFFFF0 的短表 → Err。

## 不做
- 发现 4(坐标 unwrap_or(0.0))与 8/9/10/11(utils 缓存、allow 注释、
  etkdg 守卫、optimizer 缓冲区)—— 记录在案,另行立项。
- 不动 gfnff 内嵌 JSON 初始化 unwrap 群(编译期固定数据)。

## 安全
- parser 行为变化:畸形输入从"静默容忍/panic"变为显式 Err —— 属修复
  意图;328 存量测试 + CDP 19 套件为回归线。
- total_cmp 语义:NaN 排最后而非 panic;-0.0/+0.0 全序差异对排序无影响。

## 验收
- 每提交:cargo test 全绿 + 该步新增回归测试绿。
- 终态:clippy 双目标 0 警告、fmt clean、cargo test 全量、Node 73、
  CDP 19/19(parser 触及全部 wasm 路径,必须全量)。
