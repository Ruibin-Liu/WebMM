# Plan: 全面 Review 应修项修复(6 项,零行为变化为主)

## 范围(对应 Review 发现,逐项)

1. **README 刷新**:删除 "workflow stepper" 陈旧描述(Stepper 已移除);
   核对平台描述其余指控(task ticker/概览面板/诚实横幅均存在)后修正。
2. **版本对齐**:package.json 1.10.0 → 1.11.0(对齐 Cargo.toml,
   后者为 webmm_version() 的真相源);"main": "pkg/index.js" →
   "pkg/webmm.js"(实际构建产物名)。
3. **src/lib.rs crate 级 `#![allow(unused_variables)]` 移除**:
   移除后跑 clippy,逐个修复真实警告(下划线前缀或删除未用绑定);
   clippy --all-targets 必须回到 0 警告、cargo test 全绿。
4. **db.js 两处**:`tx()` 死代码删除;`del()` resolve 点从
   req.onsuccess 改为 tx.oncomplete + onabort reject(镜像 r2 加固后
   的 put/putBatch;删除持久性同论证)。
5. **site/ CDN 本地化**:site/index.html 与 site/playground.html 的
   jsdelivr 3Dmol 改为 vendored 引用(按 pages.yml 部署结构选路径,
   产物 app/vendor/3Dmol-min.js 已在库);部署工作流如需携带该文件
   则同步补 copy 步骤。
6. **杂项**:`git rm --cached pi_resume` + .gitignore 加 pi_resume 与
   .claude/(镜像 .opencode/ 条目)。

## 不做
- 不动观察项(CODE_STATUS 归档/lib.rs 测试迁移/vendor 体积/logHash)。
- 不改任何 Rust 逻辑(第 3 项只删 allow 属性与消警,行为零变化)。

## 安全
- db.js del() 变更触及存储层:Node 73 + CDP 全量 19 套件为验收线。
- site 改动只动 script src 与工作流拷贝,不动页面逻辑。

## 验收
- cargo clippy --all-targets 0 警告;cargo test 328 全绿。
- node tests/platform/run.js 73/73;tests/cdp/run.sh 19/19。
- 两 site 页面经本地服务器确认 3Dmol 从 vendored 路径 200 加载、零外联。
- git status 干净(除本轮改动);CODE_STATUS.md 记录修复条目。
