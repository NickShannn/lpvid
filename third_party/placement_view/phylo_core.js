/* phylo_core.js — 可测试的纯逻辑核心: Newick 解析 / phylogram 布局 / query 插入 / 元数据解析 / 邻域折叠
 * 构建时被内联进 tree_explorer.html; 同时可被 node test_core.js 直接 require 测试.
 * 约定: 无 DOM 依赖, 全部函数通过 module.exports 或 window.PhyloCore 暴露.
 */
(function (global) {
"use strict";

/* ---------- Newick 解析 ----------
 * 支持: 引号标签, 内部节点 support(单值 "98" 或 "98.5/100" 形式), [&key=val] 注释,
 *       多歧树(polytomy), 结尾可选分号.
 * 节点字段: {name, length, support, supports[], meta{}, children[], up, depth, y, x}
 */
function parseNewick(text) {
  const s = String(text || "").trim();
  let i = 0;

  function label() {
    let out = "";
    if (s[i] === "'") {                      // 引号标签(内部 '' 转义)
      i += 1;
      while (i < s.length) {
        if (s[i] === "'") {
          if (s[i + 1] === "'") { out += "'"; i += 2; continue; }
          i += 1; break;
        }
        out += s[i]; i += 1;
      }
      return out.trim();
    }
    while (i < s.length && ",():;[".indexOf(s[i]) === -1) { out += s[i]; i += 1; }
    return out.trim();
  }

  function comment() {                        // [&... ] / [...] → meta 字符串
    if (s[i] !== "[") return null;
    const start = i + 1;
    while (i < s.length && s[i] !== "]") i += 1;
    const body = s.slice(start, i);
    i += 1;
    return body;
  }

  function branchLen() {
    if (s[i] !== ":") return 0;
    i += 1;
    let num = "";
    while (i < s.length && ",();[".indexOf(s[i]) === -1) { num += s[i]; i += 1; }
    return Math.max(0, Number(num) || 0);
  }

  function node() {
    const n = { name: "", length: 0, support: null, supports: null, meta: null,
                children: [], up: null, depth: 0, x: 0, y: 0 };
    if (s[i] === "(") {
      i += 1;
      n.children.push(node());
      while (s[i] === ",") { i += 1; n.children.push(node()); }
      if (s[i] !== ")") throw new Error("unbalanced tree at " + i);
      i += 1;
      const raw = label();
      if (raw) {
        const nums = raw.split("/").map(Number).filter(v => !isNaN(v));
        if (n.children.length && nums.length) { n.support = nums[nums.length - 1]; n.supports = nums; }
        else if (raw) n.name = raw;
      }
      const c1 = comment();
      if (c1) n.meta = c1;
      n.length = branchLen();
      const c2 = comment();
      if (c2) n.meta = (n.meta ? n.meta + " " : "") + c2;
      n.children.forEach(c => { c.up = n; });   // up 在解析期维护, 不依赖后续布局
      return n;
    }
    n.name = label();
    n.meta = comment();
    n.length = branchLen();
    return n;
  }

  const root = node();
  while (i < s.length && (s[i] === ";" || s[i] === " ")) i += 1;
  if (i < s.length) throw new Error("trailing content at " + i);
  // 唯一根且无名义子树 ((A,B)) → 提升
  while (!root.name && root.children.length === 1 && !root.length) {
    root.children[0].up = null;
    return root.children[0];
  }
  return root;
}

/* ---------- phylogram 布局 ----------
 * 叶子按先序定 y(0..n-1); 内节点 y=子女中点; x=depth=累计枝长.
 * collapsedSet(可选): 集合内的内节点折叠为"虚拟叶"(clade collapse)—— 子树不参与布局,
 *   节点本身占一行并标记 _collapsed; 不传时行为与原版完全一致.
 * 返回 {tips[], nTips, maxWidth} 并原地写 node.y/node.depth(node._collapsed).
 */
function layoutPhylogram(root, collapsedSet) {
  const tips = [];
  let maxWidth = 0;
  (function walk(n, depth) {
    n.depth = depth + n.length;
    if (n.depth > maxWidth) maxWidth = n.depth;
    const fold = !!(collapsedSet && collapsedSet.has(n) && n.children.length);
    n._collapsed = fold;
    if (fold) {
      // 折叠: 本节点占一行; 子树做 depth-only 传播(depth 有效供三角宽度, y 全部压到本行)
      let sd = n.depth;
      n.y = tips.length; tips.push(n);
      (function dw(m) {
        m.y = n.y;
        for (const c of m.children) {
          c.depth = m.depth + c.length;
          c._collapsed = false;
          if (c.depth > sd) sd = c.depth;
          dw(c);
        }
      })(n);
      n._subDepth = sd;
      return;
    }
    if (!n.children.length) { n.y = tips.length; tips.push(n); return; }
    for (const c of n.children) walk(c, n.depth);
    n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2;
  })(root, 0);
  return { tips, nTips: tips.length, maxWidth: maxWidth || 1 };
}

/* ---------- tip 索引 / 全树遍历 ---------- */
function indexTips(root) {
  const map = new Map();
  (function walk(n) {
    if (!n.children.length && n.name) map.set(n.name, n);
    n.children.forEach(walk);
  })(root);
  return map;
}
function allNodes(root) {
  const out = [];
  (function walk(n) { out.push(n); n.children.forEach(walk); })(root);
  return out;
}

/* ---------- query 插入(EPA-sister 近似) ----------
 * 对每条 placement: 找 nn1 tip, 将其与 query 组成新 clade:
 *   分叉点取 nn1 原枝中点 → nn1.length/2（纯绘图近似；LucaPhylo 不输出 query 枝长，不对外作为指标展示）。
 * 同一 nn1 的多条 query 并列为 polytomy(children[0] 恒为 nn1), 不层层嵌套.
 * 返回 {inserted, conflicts[]}: conflicts = nn1 不在骨架树上的 query.
 */
function insertQueries(root, placements) {
  const tips = indexTips(root);
  const inserted = [], conflicts = [];
  for (const p of placements) {
    const nn1 = p.neighbors && p.neighbors[0] && tips.get(p.neighbors[0][0]);
    if (!nn1 || !nn1.up) { conflicts.push(p.query_id); continue; }
    const q = { name: p.query_id, length: nn1.length / 2,
                support: null, supports: null, meta: null, children: [],
                up: null, isQuery: true, placement: p };
    const parent = nn1.up;
    if (parent.isQueryAnchor) {              // 同锚点已放过 query → 并列
      parent.children.push(q);
      q.up = parent;
    } else {
      const half = nn1.length / 2;
      const mid = { name: "", length: half, support: null, supports: null, meta: null,
                    children: [nn1, q], up: parent, isQueryAnchor: true };
      nn1.length = half; nn1.up = mid; q.up = mid;
      parent.children[parent.children.indexOf(nn1)] = mid;
    }
    inserted.push(q);
  }
  return { inserted, conflicts };
}

/* ---------- tip 元数据启发式(RdRP_tree_230328 命名规约) ----------
 * 形如:  NC_016896_Astrovirus_wild_boar          → RefSeq 参考
 *        YP_009345903.1_Qinghai_..._astrovirus_2  → GenBank 蛋白(注释序列)
 *        UWK09061.1_Gerbil_astrovirus_...         → 蛋白 acc + 病毒描述
 *        tombus_... / virga_... / birna-like_...  → 类群标签前缀
 * 宿主标签只来自结构化 taxonomy host_marks，绝不按名称关键词推断。
 */
function parseTipMeta(name) {
  const first = name.split("_")[0];
  const meta = { source: "other", group_tag: null, raw: name };
  if (/^N[CGYZ]_/.test(name)) { meta.source = "refseq"; return meta; }
  // 蛋白 accession 跨下划线: YP_009345903.1_ / NP_041234 / WP_... / UMO75842.1_ / BHTH16136_ / GDSWC170753_
  if (/^[A-Z]{2,3}_\d+(\.\d+)?_/.test(name) || /^[A-Z]{2,8}\d{3,}(\.\d+)?_/.test(name)) {
    meta.source = "genbank"; return meta;
  }
  // 无 accession 的描述性前缀 → 类群标签
  meta.group_tag = first.toLowerCase();
  meta.source = "group";
  return meta;
}

/* ---------- 中值根(midpoint rooting) ----------
 * 两次最远 tip 搜索定树直径 A-B, 在 A→B 路径中点插入新根:
 *   1) farthest(): 无向 DFS(边权: child.length = child→parent);
 *   2) rerootNode(p): 翻转 p 的祖先链使 p 成根;
 *   3) 切割边 (c,p) 上按 offset 拆两段, 新根 R.children=[c,p].
 * 返回新根; 原 root 变成普通内部节点.
 */
function farthest(start) {
  let best = null, bestD = -1, bestPath = null;
  (function dfs(n, from, d, path) {
    if (!n.children.length && d > bestD) { bestD = d; best = n; bestPath = path.slice(); }
    const nbs = n.children.concat(n.up ? [n.up] : []);
    for (const m of nbs) {
      if (m === from) continue;
      const w = (m.up === n) ? m.length : n.length;   // child 方向用 m.length, parent 方向用 n.length
      path.push(n);
      dfs(m, n, d + w, path);
      path.pop();
    }
  })(start, null, 0, []);
  return { tip: best, dist: bestD, path: bestPath };
}

function rerootNode(r) {
  // 使 r 成为根: 沿 r 的祖先链翻转 up/children(非链子树不动)
  let prev = null, prevEdge = 0, cur = r;
  while (cur) {
    const next = cur.up, edgeUp = cur.length;
    if (next) { const i = next.children.indexOf(cur); if (i >= 0) next.children.splice(i, 1); }
    if (prev) { prev.children.push(cur); cur.up = prev; cur.length = prevEdge; }
    else cur.up = null;
    prev = cur; prevEdge = edgeUp;
    cur = next;
  }
  return r;
}

function midpointRoot(root) {
  const first = (function anyTip(n) { return n.children.length ? anyTip(n.children[0]) : n; })(root);
  const A = farthest(first).tip;
  const B = farthest(A);
  const m = B.dist / 2;
  // 沿 A→B 路径累计找切割边 (c=A 侧, p=B 侧), offset 从 c 端量
  // farthest().path = [start, ..., 终点的父] (含起点不含终点)
  const path = B.path.concat([B]);
  let cum = 0, c = null, p = null, L = 0, offset = 0;
  for (let i = 0; i + 1 < path.length; i++) {
    const u = path[i], v = path[i + 1];
    const w = (v.up === u) ? v.length : u.length;
    if (cum + w >= m - 1e-12) { c = u; p = v; L = w; offset = m - cum; break; }
    cum += w;
  }
  if (!p) { c = path[path.length - 2]; p = path[path.length - 1]; L = 1; offset = 0.5; }
  rerootNode(p);                             // p 成根; 此时 c 必为 p 的孩子(两个方向统一)
  offset = Math.min(L, Math.max(0, offset));
  const i = p.children.indexOf(c);
  if (i >= 0) p.children.splice(i, 1);
  c.length = offset;
  const R = { name: "", length: 0, support: null, supports: null, meta: null,
              children: [c, p], up: null, depth: 0, x: 0, y: 0, isMidRoot: true };
  c.up = R; p.up = R; p.length = L - offset;
  return R;
}

const PhyloCore = { parseNewick, layoutPhylogram, indexTips, allNodes,
                    insertQueries, parseTipMeta,
                    midpointRoot, rerootNode, farthest };
if (typeof module !== "undefined" && module.exports) module.exports = PhyloCore;
global.PhyloCore = PhyloCore;
})(typeof window !== "undefined" ? window : globalThis);
