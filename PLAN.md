# Plan: Workbench file:// 打开时的裸 TypeError 修复(RDKit 未加载防御 + 可操作提示)

## 背景

用户观察:页脚 `engine: WebMM` 的版本号似乎丢了——随后自行更正:
没丢,只是引擎(7MB 级 WASM,双布局探测 + 动态导入)加载更慢,
期间静态占位就是裸 `WebMM`,与 RDKit 侧的 `Loading...` 不一致,
初看像版本号缺失。另发现:双布局两条导入路径都失败时无任何处理
(webmmready 永不触发,占位永远停留且 3D 按钮静默不可用)。

## 任务(app/index.html,两处微改)

1. #engineVersion 初始占位 `WebMM` → `WebMM (loading…)`,与 RDKit
   侧行为对齐,消除"版本号丢失"的误读;
2. 引擎加载器外层 try/catch:双路径均失败时占位改为
   `WebMM (unavailable)` 并 console.error 归因(取代静默死态)。

## 验收(实施后实测记录)

- Playwright:正常 http——早期占位 `WebMM (loading…)`(route 延迟 4s
  验证),就绪后 `WebMM v1.3.1`、3D 动作启用、零 page error;引擎双路径
  阻断——占位 `WebMM (unavailable)` + console.error,页面其余部分
  (RDKit)正常;
- **门禁附带修复(阻塞性既有脆性,与本任务无关但卡 cargo test)**:
  prop_tests::gradient_finite_difference 在随机种子下失败(proptest
  将种子持久化到 proptest-regressions/prop_tests.txt 后必复现)。
  数值定论:解析梯度正确(g2[2] = −0.08324 与闭式一致);单侧有限
  差分在该构型(拉伸键 dE/dr ≈ −6241、z 路径曲率 ∂²r/∂z² = 1/r = 2)
  的截断误差 (eps/2)·|dE/dr|·(1/r) ≈ 6.2e-4,与观测差逐位吻合——
  纯测试数值方法问题。修复:改中心差分(偶阶曲率项严格抵消,残余
  O(eps²) 截断 + ~1e-6 舍入),注释记录推导;3 个失败种子重放全过;
- `cargo test` 281/281、clippy 0、fmt 干净。
