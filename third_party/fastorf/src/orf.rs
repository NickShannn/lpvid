// 核心: 6 框 ORF 扫描与嵌套过滤。语义口径集中在本文件, 便于对拍校准 (M4)。
use crate::genetic_code::{codon_index, Code, INVALID_CODON};

pub struct Orf {
    pub orf_num: usize, // 序列内编号(1-based, 按输出顺序)
    pub strand: u8,     // b'+' | b'-'
    pub frame: u32,     // 该链上的框 1..=3 (反向链指 revcomp 序列上的框)
    pub start: usize,   // 1-based 正向坐标, 闭区间
    pub end: usize,     // 同上; circular 环绕 ORF 用线性化坐标, end 可能 > seq_len
    pub nt_len: usize,  // 链方向 5'->3' 长度, 含终止密码子
    pub aa_len: usize,  // 蛋白残基数(不含终止)
    pub has_stop: bool,
    pub edge5: bool, // 起点贴搜索窗口头 (官方 partial 语义)
    pub edge3: bool, // 终点贴搜索窗口尾
}

pub struct OrfParams {
    pub min_len: usize,    // nt, 含终止密码子口径
    pub nested_drop: bool,
    pub frames: [bool; 6], // 下标 0-2 = 正向框 1-3, 3-5 = 反向框 1-3
    pub circular: bool,
    pub seq_len: usize,
}

fn encode_base(b: u8) -> u8 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' | b'U' => 3,
        _ => 4, // N 及全部 IUPAC 模糊碱基: 不构成起始/终止, 翻译为 X
    }
}

fn comp_enc(e: u8) -> u8 {
    match e {
        0 => 3,
        1 => 2,
        2 => 1,
        3 => 0,
        _ => 4,
    }
}

// 一条链上的候选 ORF (链上 0-based 坐标)
struct RawOrf {
    q_start: usize,
    q_end: usize, // 排他末端; 有 stop 时为 stop 末碱基之后
    has_stop: bool,
    edge5: bool,
    edge3: bool,
}

// 官方 0.4.3 实测语义 (M4 校准):
// - 同框内每个 stop 区间只报第一个起始密码子的 ORF (后续 start 为同框嵌套, 不报)
// - min_len 按不含终止密码子的 nt 长度
fn scan_chain(enc: &[u8], off: usize, code: &Code, min_len: usize, chain_len: usize) -> Vec<RawOrf> {
    let mut out = Vec::new();
    let mut pending: Option<usize> = None; // 只保留区间内第一个 start
    let mut i = off;
    while i + 3 <= enc.len() {
        let c = codon_index(enc[i], enc[i + 1], enc[i + 2]);
        if c != INVALID_CODON && (code.stop_mask >> c) & 1 == 1 {
            if let Some(qs) = pending.take() {
                if i - qs >= min_len {
                    // 官方语义: 3' partial = 无终止密码子; 有 stop 即便贴窗口尾也非 partial (seqM 实测)
                    out.push(RawOrf { q_start: qs, q_end: i + 3, has_stop: true, edge5: qs == off, edge3: false });
                }
            }
        } else if c != INVALID_CODON && (code.start_mask >> c) & 1 == 1 {
            if pending.is_none() {
                pending = Some(i);
            }
        }
        i += 3;
    }
    // 链末端: 3' partial (有起始无终止)。官方口径: 3' 端截到最后一个完整密码子,
    // 不含末尾 1-2 nt 零头。
    if let Some(qs) = pending {
        let last_codon_end = off + 3 * ((chain_len - off) / 3); // <= chain_len 的最大框对齐末端
        if last_codon_end > qs && last_codon_end - qs >= min_len {
            out.push(RawOrf { q_start: qs, q_end: last_codon_end, has_stop: false, edge5: qs == off, edge3: true });
        }
    }
    out
}

pub fn find_orfs(seq: &[u8], code: &Code, p: &OrfParams) -> Vec<Orf> {
    // 官方行为: 序列(片段)长度 < min_len 时整条跳过
    if seq.len() < p.min_len {
        return Vec::new();
    }
    let mut orfs: Vec<Orf> = Vec::new();
    let enc: Vec<u8> = seq.iter().map(|&b| encode_base(b)).collect();
    let mut work = enc.clone();
    if p.circular {
        // 环状: seq ++ seq[..len-1], 允许跨原点密码子; ORF 全长仍限 <= seq_len
        let n = enc.len();
        work.extend_from_slice(&enc[..n - 1]);
    }
    let rev: Vec<u8> = work.iter().rev().map(|&e| comp_enc(e)).collect();
    for f in 0..6usize {
        if !p.frames[f] {
            continue;
        }
        let plus = f < 3;
        let chain: &[u8] = if plus { &work } else { &rev };
        let frame = (f % 3) as u32 + 1;
        let raws = scan_chain(chain, f % 3, code, p.min_len, chain.len());
        for r in raws {
            // 环状: 不输出 3' partial (环无末端); ORF 全长 <= 序列长;
            // 且起点必须落在原序列范围内 (q_start >= n 的是 work 后半的相位重叠产物,
            // n 非 3 倍数时三个框会绕环互访, 此规则给出无重复且覆盖全环的 ORF 集)
            if p.circular
                && (!r.has_stop || r.q_start >= p.seq_len || r.q_end - r.q_start > p.seq_len)
            {
                continue;
            }
            let nt_len = r.q_end - r.q_start;
            let (start, end) = if plus {
                (r.q_start + 1, r.q_end) // 1-based 闭区间
            } else {
                // rev 链上 0-based q 对应 work 正向 0-based l-1-q
                let l = rev.len();
                let a = l - r.q_end;     // 正向左端 (work 坐标)
                let b = l - 1 - r.q_start; // 正向右端
                if p.circular {
                    // 环绕序列: 折回原坐标 [0, n); 跨原点 ORF 用 end>n 的线性化表示
                    let n = p.seq_len;
                    if b < n {
                        (a + 1, b + 1)
                    } else if a >= n {
                        (a - n + 1, b - n + 1)
                    } else {
                        (a + 1, b + 1) // a < n <= b: 跨原点
                    }
                } else {
                    (a + 1, b + 1)
                }
            };
            orfs.push(Orf {
                orf_num: 0,
                strand: if plus { b'+' } else { b'-' },
                frame,
                start,
                end,
                nt_len,
                aa_len: nt_len / 3 - if r.has_stop { 1 } else { 0 },
                has_stop: r.has_stop,
                edge5: r.edge5,
                edge3: r.edge3,
            });
        }
    }
    if p.nested_drop {
        drop_nested(&mut orfs);
    }
    for (i, o) in orfs.iter_mut().enumerate() {
        o.orf_num = i + 1;
    }
    orfs
}

/// ORF 的核酸序列 (链方向 5'->3', 字母 ACGTN, 含终止密码子)。
/// 环绕 ORF (end > seq_len) 按线性化坐标取模拼接。
pub fn cds_of(seq: &[u8], o: &Orf) -> Vec<u8> {
    let len = seq.len();
    let (a, b) = (o.start - 1, o.end); // 0-based 半开
    if o.strand == b'+' {
        if b <= len {
            seq[a..b].to_vec()
        } else {
            let mut v = seq[a..].to_vec();
            v.extend_from_slice(&seq[..b - len]);
            v
        }
    } else {
        // minus: 链方向 = 正向坐标递减取互补
        let rev_range = |lo: usize, hi: usize| -> Vec<u8> {
            seq[lo..hi]
                .iter()
                .rev()
                .map(|&b| match b {
                    b'A' => b'T',
                    b'T' | b'U' => b'A',
                    b'G' => b'C',
                    b'C' => b'G',
                    _ => b'N',
                })
                .collect()
        };
        if b <= len {
            rev_range(a, b)
        } else {
            // 环绕: 环上链序 = revcomp(seq[0..b-len]) 接 revcomp(seq[a..len])
            let mut v = rev_range(0, b - len);
            v.extend(rev_range(a, len));
            v
        }
    }
}

/// ORF 翻译为蛋白 (起始密码子按惯例译 M, 待 M4 校准; 终止不入; 含 N 密码子译 X)。
pub fn prot_of(seq: &[u8], o: &Orf, code: &Code) -> Vec<u8> {
    let cds = cds_of(seq, o);
    let n_codons = cds.len() / 3;
    let mut prot = Vec::with_capacity(n_codons);
    for k in 0..n_codons {
        let c = codon_index(
            encode_base(cds[k * 3]),
            encode_base(cds[k * 3 + 1]),
            encode_base(cds[k * 3 + 2]),
        );
        if k == n_codons - 1 && o.has_stop && c != INVALID_CODON && (code.stop_mask >> c) & 1 == 1 {
            break; // 终止密码子不入蛋白
        }
        let aa = if c == INVALID_CODON {
            b'X'
        } else {
            code.aa[c]
        };
        // 官方口径(219万条全量对拍确认): 起始密码子译 M; 但起点贴窗口头(5' partial,
        // edge5)时保留原义 (TTG/CTG->L 等), 因为无法确认其为真实起始
        prot.push(if k == 0 && !o.edge5 { b'M' } else { aa});
    }
    prot
}

/// 嵌套删除: 坐标区间完全包含即删, 不分框/链 (官方 help 字面语义; M4 对拍校准)。
/// 栈算法: 按 (start asc, end desc) 排序, 被任一已保留区间完全包含者删。
fn drop_nested(orfs: &mut Vec<Orf>) {
    let mut order: Vec<usize> = (0..orfs.len()).collect();
    order.sort_by_key(|&i| (orfs[i].start, std::cmp::Reverse(orfs[i].end)));
    let mut keep = vec![false; orfs.len()];
    let mut stack: Vec<usize> = Vec::new(); // 已保留区间的索引
    for &i in &order {
        while let Some(&top) = stack.last() {
            if orfs[top].end < orfs[i].start {
                stack.pop(); // 不可能再包含任何后续区间
            } else {
                break;
            }
        }
        let nested = stack
            .last()
            .map(|&top| orfs[i].end <= orfs[top].end)
            .unwrap_or(false);
        if nested {
            keep[i] = false;
        } else {
            keep[i] = true;
            stack.push(i);
        }
    }
    let mut out = Vec::with_capacity(keep.iter().filter(|k| **k).count());
    for (i, o) in orfs.drain(..).enumerate() {
        if keep[i] {
            out.push(o);
        }
    }
    *orfs = out;
}
