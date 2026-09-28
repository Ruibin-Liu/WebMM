# Plan: Search 库持久化(localStorage)

## 背景

Search 标签页的库目前刷新即失(内存态),用户每次都要重新粘贴/上传。
用户确认持久化有必要。约束:指纹位串(2000 条 × 5 × 2048 bit ≈ 20MB)
远超 localStorage 配额——**持久化输入(canonical SMILES + 名称,~100-200KB),
不持久化指纹缓存**;恢复时重算(2000 条约 1s,一次性)。

## 设计

1. **保存**:`loadSearchLibraryFrom` 成功后,把解析结果(canonical SMILES +
   名称,JSON)写 `localStorage['wb-searchLibrary']`(配额失败 try/catch,
   状态行提示"library loaded but too large to persist")
2. **恢复**:首次 `switchMode('search')` 时,若内存无库且存在已存库 →
   自动装载,状态行注明 "restored from your last session"
   (不恢复 textarea 内容——textarea 三模式共享)
3. **清除**:searchBar 加 "Clear" 按钮 = 清内存库 + 清 localStorage +
   收起结果面板;状态行确认
4. Demo library 同样走保存路径(它就是当前库)

## 验收(Playwright)

1. demo 装载 → 刷新页面 → 切 Search:库自动恢复(计数一致)、检索可用
2. 自定义库(粘贴 2 条)→ 刷新 → 恢复 2 条;检索正常
3. Clear → 刷新 → 不恢复,状态行正确;textarea 不受影响
4. localStorage 损坏(手写非法 JSON)→ 不抛错,按无库处理
5. 390px / 零 page error / 门禁 sanity(281 tests、clippy、fmt)

## 边界

- 不做库管理 UI(多库命名保存等);不持久化检索历史
