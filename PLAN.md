# Plan: 指纹选择器按检索方式门控(用户发现 sub 模式可选指纹)

## 缺陷
- fp 只在 sim 分支使用(Tanimoto 需要指纹);sub 是 SMARTS 子结构匹配、
  shape 走构象叠合 —— 两页 onSearchModeChange 都只对 shape 禁用
  searchFp,sub 模式下指纹下拉仍可选,误导用户。

## 改动(app/platform.html + app/index.html 各一行)
- searchFp.disabled = (mode !== 'sim')。

## 验收
- m6_platform 检索段加断言:切 sub → searchFp.disabled===true;切回 sim →
  false;shape 仍禁用。m5 同款(index 页)。全量门复跑。
