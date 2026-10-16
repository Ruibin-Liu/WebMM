# Plan: scripts/build-ketcher.sh(Ketcher bundle 可复现构建 + 版本 pin)

## 背景
- app/vendor/ketcher/ketcher-bundle.js(31MB IIFE)是一次性 esbuild 产物,
  构建命令未入库 —— 版本升级/CVE 修复时不可复现。
- 入口封装 WebMMKetcher(mount/getMolfile/setMolecule/unmount)只存在于
  编译产物中,源码不在库 —— 需一并还原为 scripts/ketcher-entry.js。

## 范围
- 新增 scripts/build-ketcher.sh:临时目录 npm ci 精确版本
  (ketcher-react 3.18.0 / ketcher-standalone 3.18.0 / esbuild pin),
  esbuild --bundle --format=iife --minify + 守卫式 process banner,
  产出 ketcher-bundle.js + ketcher.css + LICENSE 到 app/vendor/ketcher/。
- 新增 scripts/ketcher-entry.js:WebMMKetcher 入口封装源码
  (从现有 bundle 逆向还原,API 与 platform.html 调用面一致)。
- 不改动现有 vendor 产物、platform.html、页面逻辑;默认输出到临时目录,
  --install 才写入 app/vendor/ketcher/。

## 验收
- 脚本产出 bundle 含守卫 banner 与 WebMMKetcher 全 API;
- 用重建 bundle 替换后 Ketcher 模态 mount/setMolecule/getMolfile 冒烟通过
  (Playwright,临时替换、测完还原);
- 版本标记 3.18.0 与现有产物一致。
