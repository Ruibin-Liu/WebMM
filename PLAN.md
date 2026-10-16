# Plan: 测试基建固化(README 重写 + run.sh + Upcoming 清理)

## 事实(已核实)
- 12 个 CDP 测试(m0-m6 + 5 个专项),README 只描述到 m6 旧版。
- 真实外部输入仅两个:/tmp/caff24.sdf(m1;重启即丢,RDKit 重生配方在
  旧 README)、/tmp/m5-project.json(m6 运行时导出;import_loop 与
  rdkit_race 消费,缺失时它们各自报错退出 rc=2/直接失败)。
  其余 /tmp 路径均为测试自产输出(下载/截图)。
- 环境变量:LOOP_N(import_loop 迭代数,默认 10)、STOP_ON_FAIL(首败
  停 + 诊断);scripts/build-ketcher.sh 另有 KEEP。
- 服务器:repo 根 8901 端口。
- CODE_STATUS ## Upcoming 两项均过期:site SDF 修复(site/index.html
  无 test.sdf 残留)、Pages 部署(pages.yml 在跑)。

## 改动
1. tests/cdp/README.md 重写:套件总表(名称/一句话/输入依赖)、夹具
   依赖与重生配方、环境变量、运行顺序、run.sh 用法。
2. 新增 tests/cdp/run.sh:检查/按需起服(只杀自己起的)、按依赖序跑
   全套(m0-m6 → platform Node → 专项),可选 -k 单套参数;set -e 语义
   但汇总各套结果。
3. CODE_STATUS.md:## Upcoming 清为空或注明已迁移完成;新增本轮条目。

## 验收
- run.sh 全量跑通(绿);缺 caff24 时自动重生再跑 m1;README 与实际
  文件/变量一一对应(人工核对)。
