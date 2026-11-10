/* test_core.js — node 单测: 用真实剪枝树 + bundle 验证 phylo_core 纯函数.
 * 运行: node test_core.js  (退出码 0 = 全绿)
 */
"use strict";
const fs = require("fs");
const path = require("path");
const PC = require("./phylo_core.js");

let pass = 0, fail = 0;
function ok(cond, msg) {
  if (cond) { pass += 1; }
  else { fail += 1; console.error("  FAIL: " + msg); }
}

/* --- 1. Newick 解析: 语法特性 --- */
{
  const t1 = PC.parseNewick("(A:0.1,(B:0.2,(C:0.3,D:0.4)98/100:0.05)'clade x':0.1);");
  ok(t1.children.length === 2, "二叉根");
  ok(t1.children[0].name === "A" && Math.abs(t1.children[0].length - 0.1) < 1e-12, "tip 名与枝长");
  const clade = t1.children[1];
  ok(clade.name === "clade x", "引号标签");
  const cd = clade.children[1];
  ok(cd.support === 100 && cd.supports.join(",") === "98,100", "SH-aLRT/UFBoot 双支持值");
  ok(cd.children[0].name === "C" && cd.children[1].name === "D", "嵌套结构");
  const t2 = PC.parseNewick("((A:1,B:1):5,C:2)[&color=red];");
  ok(t2.meta === "&color=red", "根注释");
  const t3 = PC.parseNewick("((A:0.1,B:0.2));");
  ok(t3.children.length === 2, "唯一根包装提升");
  let threw = false;
  try { PC.parseNewick("(A:0.1,B:0.2"); } catch (e) { threw = true; }
  ok(threw, "不平衡树抛错");
}

/* --- 2. 布局不变量 --- */
{
  const t = PC.parseNewick("(A:0.1,(B:0.2,C:0.3)90:0.05);");
  const lay = PC.layoutPhylogram(t);
  ok(lay.nTips === 3, "tip 计数");
  ok(new Set(lay.tips.map(x => x.y)).size === 3, "y 唯一");
  const ys = lay.tips.map(x => x.y);
  ok(ys.every(v => v >= 0 && v < 3), "y 在 [0,n)");
  const inner = t.children[1];
  ok(Math.abs(inner.y - (1 + 2) / 2) < 1e-12, "内节点=子女中点");
  ok(lay.tips[1].depth === 0.05 + 0.2 && lay.tips[0].depth === 0.1, "depth 累计枝长");
}

/* --- 3. midpoint rooting(中值根) --- */
{
  // 手算例 1: ((A:1,B:1):1,(C:2,D:2):3) — A-C patristic=7, 中点 3.5 落 R0→CD 边(offset 1.5)
  const t1 = PC.parseNewick("((A:1,B:1):1,(C:2,D:2):3);");
  const tips1 = PC.indexTips(t1);
  const r1 = PC.midpointRoot(t1);
  const lay1 = PC.layoutPhylogram(r1);
  ok(r1.children && r1.children.length === 2, "mid1: 新根二叉");
  const idx1 = PC.indexTips(r1);
  ok(idx1.size === 4 && [...idx1.keys()].sort().join() === "A,B,C,D", "mid1: tip 集不变");
  ok(_upChainOk(r1), "mid1: up 链完整无环");
  ok(Math.abs(idx1.get("A").depth - 3.5) < 1e-9 && Math.abs(idx1.get("C").depth - 3.5) < 1e-9,
     `mid1: 直径两端等深 3.5 (A=${idx1.get("A").depth.toFixed(4)}, C=${idx1.get("C").depth.toFixed(4)})`);
  // 手算例 2(polytomy): ((A:1,B:2,C:3):4,D:10) — C-D=17, 中点 8.5 落 R0→D 边 offset 1.5
  const t2 = PC.parseNewick("((A:1,B:2,C:3):4,D:10);");
  const r2 = PC.midpointRoot(t2);
  PC.layoutPhylogram(r2);
  const idx2 = PC.indexTips(r2);
  ok(_upChainOk(r2), "mid2: up 链完整");
  ok(Math.abs(idx2.get("C").depth - 8.5) < 1e-9 && Math.abs(idx2.get("D").depth - 8.5) < 1e-9,
     `mid2: 直径两端等深 8.5 (C=${idx2.get("C").depth.toFixed(4)}, D=${idx2.get("D").depth.toFixed(4)})`);
  // 手算例 3(中点落在内边): ((A:1,B:4):2,(C:1,D:1):9) — 直径 B-C=16(B=6,C=10),
  // 中点 8 落 R0→CD 边(累计 6→15, offset 2); A 深=2+2+1=5, D=7+1=8
  const t3 = PC.parseNewick("((A:1,B:4):2,(C:1,D:1):9);");
  const r3 = PC.midpointRoot(t3);
  PC.layoutPhylogram(r3);
  const idx3 = PC.indexTips(r3);
  ok(_upChainOk(r3) && Math.abs(idx3.get("B").depth - 8) < 1e-9 && Math.abs(idx3.get("C").depth - 8) < 1e-9,
     `mid3: 中点在内边, 直径两端等深 8 (B=${idx3.get("B").depth.toFixed(4)}, C=${idx3.get("C").depth.toFixed(4)})`);
  // 幂等性: 再 midpoint 一次 tip 深度不变
  const r3b = PC.midpointRoot(PC.parseNewick("((A:1,B:4):2,(C:1,D:1):9);"));
  const r3c = PC.midpointRoot(r3b);
  const l3b = PC.layoutPhylogram(r3b), l3c = PC.layoutPhylogram(r3c);
  ok(l3b.tips.every(t => Math.abs(t.depth - PC.indexTips(r3c).get(t.name).depth) < 1e-9),
     "mid: 幂等(重复中值根深度不变)");
  // 真实树: Astro 骨架, 性质(直径两端等深 + tip 不变)
  const fs = require("fs");
  const nwkA = fs.readFileSync(path.join(__dirname, "data", "backbone_Astro-Poty.nwk"), "utf8");
  const treeA = PC.parseNewick(nwkA);
  const before = [...PC.indexTips(treeA).keys()].sort();
  const rA = PC.midpointRoot(treeA);
  const after = [...PC.indexTips(rA).keys()].sort();
  ok(before.length === after.length && before.every((v, i) => v === after[i]),
     `Astro: 中值根后 tip ${after.length} 不变`);
  ok(_upChainOk(rA), "Astro: up 链完整");
  // 平衡性: 新根两子树最大深度差应远小于原根
  const bal = n => Math.abs(_maxDepth(n.children[0]) - _maxDepth(n.children[1]));
  const tA0 = PC.parseNewick(nwkA);
  PC.layoutPhylogram(tA0);
  PC.layoutPhylogram(rA);
  ok(bal(rA) < bal(tA0), `Astro: 根平衡改善 (${bal(tA0).toFixed(3)} → ${bal(rA).toFixed(3)})`);

  function _upChainOk(root) {
    let okAll = true;
    (function w(n) {
      for (const c of n.children) {
        if (c.up !== n) okAll = false;
        w(c);
      }
    })(root);
    if (root.up !== null) okAll = false;
    return okAll;
  }
  function _maxDepth(n) {
    let mx = n.depth;
    (function w(x) { if (x.depth > mx) mx = x.depth; x.children.forEach(w); })(n);
    return mx;
  }
}

/* --- 3.5 布局折叠(clade collapse): 折叠集内节点成为虚拟叶 --- */
{
  const t = PC.parseNewick("((A:1,B:2):3,(C:4,D:5):6);");
  const lay0 = PC.layoutPhylogram(t);
  const tree2 = PC.parseNewick("((A:1,B:2):3,(C:4,D:5):6);");
  const inner2 = tree2.children[0];            // (A,B) 子树(注意: 用第二棵树自身的节点对象)
  const lay = PC.layoutPhylogram(tree2, new Set([inner2]));
  ok(lay.nTips === 3, `折叠后 tips=3(虚拟叶) (${lay.nTips})`);
  ok(inner2._collapsed === true, "折叠节点标记 _collapsed");
  ok(new Set(lay.tips.map(x => x.y)).size === 3, "y 连续无空洞");
  ok(lay.tips.some(x => x === inner2), "折叠节点本身占一行");
  ok(!lay0.tips.some(x => x === t.children[0]), "未折叠时内节点不在 tips");
  // _subDepth: 折叠子树 maxDepth(depth-only 传播), 供三角底宽(A 深 4, B 深 5 → max 5)
  ok(inner2._subDepth === 5, `折叠节点 _subDepth=子树最大深度 5 (${inner2._subDepth})`);
  ok(Math.abs(inner2.children[0].depth - 4) < 1e-9, "子树 depth 仍有效(A 深 4)");
  // 幂等兼容: 不传 set 行为与原版一致
  const layOld = PC.layoutPhylogram(PC.parseNewick("((A:1,B:2):3,(C:4,D:5):6);"));
  ok(layOld.nTips === 4, `无折叠集时 tips=4 向后兼容 (${layOld.nTips})`);
}

/* --- 3.6 宿主大类名称关键词推断已移除(taxonomy 图标才是宿主证据) --- */
{
  ok(typeof PC.HOST_CLASSES === "undefined", "HOST_CLASSES 已移除");
  ok(typeof PC.HOST_LABEL === "undefined", "HOST_LABEL 已移除");
  ok(typeof PC.hostClassForTip === "undefined", "hostClassForTip 已移除");
  for (const name of ["Mammal_some_virus", "Cartilaginous_DTHH2075", "Bat_coronavirus", "Plant_virus"]) {
    const meta = PC.parseTipMeta(name);
    ok(!Object.prototype.hasOwnProperty.call(meta, "host_class"), `${name}: parseTipMeta 不得推断宿主大类`);
    ok(meta.source !== "env", `${name}: 不得因名称关键词标记 env`);
  }
}

/* --- 4. 真实数据: 两个超群 --- */
for (const sg of ["Astro-Poty", "Picorna-Calici"]) {
  const nwk = fs.readFileSync(path.join(__dirname, "data", `backbone_${sg}.nwk`), "utf8");
  const bundle = JSON.parse(fs.readFileSync(path.join(__dirname, "data", `bundle_virid_${sg}.json`), "utf8"));
  const tree = PC.parseNewick(nwk);
  const idx = PC.indexTips(tree);
  ok(idx.size === bundle.n_backbone, `${sg}: 骨架 tip=${idx.size} 与 bundle 一致`);
  ok(PC.allNodes(tree).length >= idx.size, `${sg}: 节点数>=tip数`);

  // 骨架不含 query
  const qids = bundle.placements.map(p => p.query_id);
  ok(qids.every(q => !idx.has(q)), `${sg}: query 均不在骨架树`);

  // 插入
  const ins = PC.insertQueries(tree, bundle.placements);
  ok(ins.inserted.length === bundle.n_query && ins.conflicts.length === 0,
     `${sg}: 插入 ${ins.inserted.length}/${bundle.n_query}, 冲突 ${ins.conflicts.length}`);
  const lay = PC.layoutPhylogram(tree);
  ok(lay.nTips === bundle.n_backbone + bundle.n_query, `${sg}: 插入后 tip=${lay.nTips}`);

  // query 与 nn1 成 sister
  let sisterOk = 0;
  for (const q of ins.inserted) {
    const mid = q.up;
    okSib: for (const c of mid.children) if (c !== q && c.name === q.placement.neighbors[0][0]) { sisterOk += 1; break okSib; }
  }
  ok(sisterOk === bundle.n_query, `${sg}: 全部 query 与 nn1 成 sister (${sisterOk})`);

  // 深度守恒: 插入后 nn1 tip 的深度变化 ≤ 其原枝长(分叉近似不炸深度)
  const idx2 = PC.indexTips(tree);
  let depthOk = true;
  const lay0 = PC.layoutPhylogram(PC.parseNewick(nwk));
  for (const t0 of lay0.tips) {
    const t2n = idx2.get(t0.name);
    if (t2n && Math.abs(t2n.depth - t0.depth) > 0.5001) { depthOk = false; }
  }
  ok(depthOk, `${sg}: 插入后骨架 tip 深度漂移均 <0.5(sister 分叉近似)`);

  // 元数据
  const metas = lay0.tips.map(t => PC.parseTipMeta(t.name));
  ok(metas.every(m => m.source !== undefined), `${sg}: 元数据解析全覆盖`);
  const env = metas.filter(m => m.source === "env").length;
  const refseq = metas.filter(m => m.source === "refseq").length;
  console.log(`  [${sg}] meta 分布: refseq=${refseq} env=${env} ` +
              `group=${metas.filter(m => m.source === "group").length} ` +
              `genbank=${metas.filter(m => m.source === "genbank").length}`);
}

console.log(`\n${pass} passed, ${fail} failed`);
process.exit(fail ? 1 : 0);
