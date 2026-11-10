// FASTA 读入(.gz 自动识别)与 ID 规范化。
use std::fs::File;
use std::io::{BufRead, BufReader, Read};

pub struct FastaRec {
    pub id: String,   // header 首个空白前的 token
    pub desc: String, // header 其余部分(可为空)
    pub seq: Vec<u8>, // 已大写; 保留全部字符(含 N 等模糊碱基), 编码在 orf.rs 做
}

pub fn read_fasta(path: &str) -> Result<Vec<FastaRec>, String> {
    let mut raw = Vec::new();
    if path == "-" {
        std::io::stdin()
            .read_to_end(&mut raw)
            .map_err(|e| format!("读 stdin 失败: {e}"))?;
    } else {
        File::open(path)
            .map_err(|e| format!("打开 {path} 失败: {e}"))?
            .read_to_end(&mut raw)
            .map_err(|e| format!("读 {path} 失败: {e}"))?;
    }
    // gzip 魔数 1f 8b: 解压到内存再解析 (输入本就整体读入, 无额外代价)
    let text: Vec<u8> = if raw.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::new();
        flate2::read::GzDecoder::new(&raw[..])
            .read_to_end(&mut out)
            .map_err(|e| format!("解压 {path} 失败: {e}"))?;
        out
    } else {
        raw
    };
    parse(BufReader::new(&text[..]), path)
}

fn parse(mut r: impl BufRead, path: &str) -> Result<Vec<FastaRec>, String> {
    let mut recs = Vec::new();
    let mut cur: Option<FastaRec> = None;
    let mut line = Vec::new();
    let mut lineno = 0usize;
    loop {
        line.clear();
        let n = r
            .read_until(b'\n', &mut line)
            .map_err(|e| format!("读 {path} 第 {} 行失败: {e}", lineno + 1))?;
        if n == 0 {
            break;
        }
        lineno += 1;
        let l = trim_eol(&line);
        if l.is_empty() {
            continue;
        }
        if l[0] == b'>' {
            if let Some(c) = cur.take() {
                recs.push(c);
            }
            let h = String::from_utf8_lossy(&l[1..]).to_string();
            let (id, desc) = match h.find(|c: char| c.is_whitespace()) {
                Some(i) => (h[..i].to_string(), h[i..].trim().to_string()),
                None => (h.trim().to_string(), String::new()),
            };
            if id.is_empty() {
                return Err(format!("{path} 第 {lineno} 行: 空 FASTA ID"));
            }
            cur = Some(FastaRec {
                id,
                desc,
                seq: Vec::new(),
            });
        } else {
            let c = cur
                .as_mut()
                .ok_or_else(|| format!("{path} 第 {lineno} 行: 序列出现在首个 header 之前"))?;
            c.seq.extend(l.iter().filter(|b| !b.is_ascii_whitespace()).map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(c) = cur.take() {
        recs.push(c);
    }
    if recs.is_empty() {
        return Err(format!("{path}: 未发现任何 FASTA 记录"));
    }
    for (i, r) in recs.iter().enumerate() {
        if r.seq.is_empty() {
            return Err(format!("{path}: 第 {} 条序列 '{}' 为空", i + 1, r.id));
        }
    }
    Ok(recs)
}

fn trim_eol(l: &[u8]) -> &[u8] {
    let mut end = l.len();
    while end > 0 && (l[end - 1] == b'\n' || l[end - 1] == b'\r') {
        end -= 1;
    }
    &l[..end]
}

/// keep 模式: 白名单 [A-Za-z0-9._-] 之外替换为 '_', 清洗后为空则回退 fallback。
pub fn sanitize_id(raw: &str, fallback: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    let mut last_underscore = false;
    for c in raw.chars() {
        if c.is_ascii_alphanumeric() || c == '.' || c == '_' || c == '-' {
            out.push(c);
            last_underscore = false;
        } else if !last_underscore {
            out.push('_'); // 连续特殊字符只折叠为一个下划线
            last_underscore = true;
        }
    }
    let trimmed = out.trim_matches('_');
    if trimmed.is_empty() {
        fallback.to_string()
    } else {
        trimmed.to_string()
    }
}
