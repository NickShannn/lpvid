// 输出格式: ORF 蛋白 FASTA / CDS FASTA / Text ASN.1 / feature table / TSV 映射表 / GFF3。
// 全部函数写进 &mut Vec<u8>, 格式化在各线程内完成, 主线程只做拼接 (性能关键)。
// ASN.1 与 feature table 为近似格式, 待 M4 拿官方样本逐字节对齐。
use crate::genetic_code::Code;
use crate::orf::{cds_of, prot_of, Orf};
use std::io::{Result as IoResult, Write};

pub struct SeqView<'a> {
    pub seq_id: &'a str,  // 规范化后的序列 ID(不含 _ORFn), 已保证仅含安全字符
    pub orig_id: &'a str, // 原始 FASTA ID
    pub desc: &'a str,
    pub orfs: &'a [Orf],
    pub code_id: u32,
    pub code: &'a Code,
    pub seq: &'a [u8], // 待输出序列(片段裁剪后), 用于按需重建 CDS/蛋白
}

// 官方 0.4.3 语义: partial = ORF 端点紧贴搜索窗口边界 (-b/-e 或序列全长), 与有无 stop 无关
fn status(o: &Orf) -> &'static str {
    match (o.edge5, o.edge3) {
        (false, false) => "complete",
        (true, false) => "5partial",
        (false, true) => "3partial",
        (true, true) => "partial",
    }
}

fn write_wrapped(w: &mut dyn Write, sv: &SeqView, o: &Orf, desc_rest: &str, seq: &[u8]) -> IoResult<()> {
    w.write_all(b">")?;
    write!(w, "{}_ORF{}", sv.seq_id, o.orf_num)?;
    if !desc_rest.is_empty() {
        write!(w, " {desc_rest}")?;
    }
    writeln!(w)?;
    for chunk in seq.chunks(60) {
        w.write_all(chunk)?;
        writeln!(w)?;
    }
    Ok(())
}

pub fn write_orf_fasta(w: &mut dyn Write, sv: &SeqView) -> IoResult<()> {
    for o in sv.orfs {
        let desc = format!(
            "[{}] Nucleotide {}..{} strand {} frame {} ({}, {} aa)",
            sv.orig_id, o.start, o.end, o.strand as char, o.frame, status(o), o.aa_len
        );
        write_wrapped(w, sv, o, &desc, &prot_of(sv.seq, o, sv.code))?;
    }
    Ok(())
}

pub fn write_cds_fasta(w: &mut dyn Write, sv: &SeqView) -> IoResult<()> {
    for o in sv.orfs {
        let desc = format!(
            "[{}] Nucleotide {}..{} strand {} frame {} ({}, {} nt)",
            sv.orig_id, o.start, o.end, o.strand as char, o.frame, status(o), o.nt_len
        );
        write_wrapped(w, sv, o, &desc, &cds_of(sv.seq, o))?;
    }
    Ok(())
}

pub fn write_tsv(w: &mut dyn Write, sv: &SeqView) -> IoResult<()> {
    for o in sv.orfs {
        let (startc, stopc) = codons_at(sv.seq, o);
        w.write_all(sv.seq_id.as_bytes())?;
        write!(w, "_ORF{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            o.orf_num,
            tab_safe(sv.orig_id),
            tab_safe(sv.desc),
            o.strand as char,
            o.frame,
            o.start,
            o.end,
            o.nt_len,
            o.aa_len,
        )?;
        w.write_all(&startc)?;
        if o.has_stop {
            w.write_all(b"\t")?;
            w.write_all(&stopc)?;
        } else {
            w.write_all(b"\t")?;
        }
        writeln!(w, "\t{}", status(o))?;
    }
    Ok(())
}

// 取 ORF 首/末密码子 (TSV 用; 栈数组零分配)。3' partial 无 stop, 返回全 0 占位。
// 索引按模 seq.len() 取, 兼容 circular 环绕 ORF。
fn codons_at(seq: &[u8], o: &Orf) -> ([u8; 3], [u8; 3]) {
    let comp = |b: u8| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        _ => b'N',
    };
    let n = seq.len();
    let at = |i: usize| {
        // circular 环绕 ORF 的索引按模 n 取 (与 cds_of 的环绕拼接一致)
        [seq[i % n], seq[(i + 1) % n], seq[(i + 2) % n]]
    };
    // minus 链: 取正向切片后反转+互补
    let rcat = |i: usize| {
        let r = at(i);
        [comp(r[2]), comp(r[1]), comp(r[0])]
    };
    if o.strand == b'+' {
        let s = at(o.start - 1);
        let e = if o.has_stop { at(o.end - 3) } else { [0, 0, 0] };
        (s, e)
    } else {
        // minus 链 5' 起始密码子位于正向 end 处, 占 [end-3, end)
        let s = rcat(o.end - 3);
        let e = if o.has_stop { rcat(o.start - 1) } else { [0, 0, 0] };
        (s, e)
    }
}

// TSV 单元格安全: 制表符/换行替换为空格 (记录于 README)
fn tab_safe(s: &str) -> String {
    s.chars()
        .map(|c| if c == '\t' || c == '\n' || c == '\r' { ' ' } else { c })
        .collect()
}

pub fn write_gff3(w: &mut dyn Write, sv: &SeqView) -> IoResult<()> {
    // seq_id 与 orf_id 均为白名单字符, 无需 percent-encoding; orig_id 需要
    for o in sv.orfs {
        w.write_all(sv.seq_id.as_bytes())?;
        write!(w, "\tfastorf\tORF\t{}\t{}\t.\t{}\t0\tID={}_ORF{};orig_id={}\n",
            o.start, o.end, o.strand as char, sv.seq_id, o.orf_num, pct(sv.orig_id))?;
    }
    Ok(())
}

// GFF3 属性值只允许安全字符, 其余 percent-encoding
fn pct(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for b in s.bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'.' | b'-' | b'_' | b':' | b'@' => {
                out.push(b as char)
            }
            _ => out.push_str(&format!("%{b:02X}")),
        }
    }
    out
}

// ---- 以下两种格式为近似实现, M4 对拍后校准 ----

pub fn write_asn1(w: &mut dyn Write, sv: &SeqView) -> IoResult<()> {
    for o in sv.orfs {
        writeln!(w, "Seq-annot ::= {{")?;
        writeln!(w, "  type feat-table,")?;
        writeln!(w, "  data ftable {{")?;
        writeln!(w, "    {{")?;
        writeln!(w, "      data cdregion {{")?;
        writeln!(w, "        code {{ id {} }},", sv.code_id)?;
        writeln!(w, "        frame {},", (o.frame - 1) % 3)?;
        writeln!(w, "        partial-product {},", !o.has_stop)?;
        writeln!(w, "      }},")?;
        // NCBI int-location 为 0-based 闭区间
        writeln!(
            w,
            "      location int {{ from {}, to {}, strand {} }},",
            o.start - 1,
            o.end - 1,
            if o.strand == b'+' { "plus" } else { "minus" }
        )?;
        write!(w, "      comment \"{}_ORF{}\",\n", sv.seq_id, o.orf_num)?;
        writeln!(w, "    }}")?;
        writeln!(w, "  }}")?;
        writeln!(w, "}}")?;
    }
    Ok(())
}

pub fn write_ftable(w: &mut dyn Write, sv: &SeqView) -> IoResult<()> {
    write!(w, ">Feature {}\n", sv.seq_id)?;
    for o in sv.orfs {
        if o.strand == b'+' {
            writeln!(w, "{}\t{}\tORF", o.start, o.end)?;
        } else {
            writeln!(w, "{}\t{}\tORF", o.end, o.start)?; // 反向: 起<止 交换
        }
        writeln!(w, "\t\t\tcodon_start\t1")?;
        if o.edge5 || o.edge3 {
            writeln!(w, "\t\t\tpartial\t{}", status(o))?;
        }
        write!(w, "\t\t\tproduct\t{}_ORF{}\n", sv.seq_id, o.orf_num)?;
    }
    Ok(())
}
