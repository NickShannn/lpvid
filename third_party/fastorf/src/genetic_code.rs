// 遗传密码表查询: 从生成的 TABLES 构建 64 项翻译表 + 起始/终止密码子位掩码。
// 密码子索引 = b1*16 + b2*4 + b3, 碱基编码 A=0,C=1,G=2,T=3; 含模糊碱基的密码子为 255。
use crate::genetic_code_data::{CodeTable, TABLES};

pub const INVALID_CODON: usize = 255;

pub struct Code {
    pub aa: [u8; 64],    // 氨基酸单字母或 b'*'
    pub start_mask: u64, // bit i = 密码子 i 可作为起始
    pub stop_mask: u64,  // bit i = 密码子 i 为终止
}

pub fn codon_index(b1: u8, b2: u8, b3: u8) -> usize {
    // 入参已是 0-3 编码; 任一为 4(N) 则无效
    if b1 > 3 || b2 > 3 || b3 > 3 {
        INVALID_CODON
    } else {
        ((b1 as usize) << 4) | ((b2 as usize) << 2) | b3 as usize
    }
}

pub fn lookup_table(id: u32) -> Result<&'static CodeTable, String> {
    TABLES
        .iter()
        .find(|t| t.id == id)
        .ok_or_else(|| format!("未知遗传密码表 id {id}; 可用: {}", ids_csv()))
}

fn ids_csv() -> String {
    TABLES.iter().map(|t| t.id.to_string()).collect::<Vec<_>>().join(",")
}

/// start_mode: 0=仅ATG, 1=ATG+表内替代起始密码子, 2=任意 sense 密码子。
/// extra_starts/extra_stops: 显式密码子列表(如 "ATG,GTG"), 覆盖对应预设。
pub fn build(
    id: u32,
    start_mode: u8,
    extra_starts: Option<&[String]>,
    extra_stops: Option<&[String]>,
) -> Result<Code, String> {
    let t = lookup_table(id)?;
    let mut aa = [0u8; 64];
    aa.copy_from_slice(t.aa.as_bytes());
    let mut stop_mask = 0u64;
    for i in 0..64 {
        if aa[i] == b'*' {
            stop_mask |= 1u64 << i;
        }
    }
    let mut start_mask = match start_mode {
        0 => {
            let i = codon_index(0, 3, 2); // ATG
            1u64 << i
        }
        1 => {
            let mut m = 0u64;
            for i in 0..64 {
                if t.starts.as_bytes()[i] == b'M' {
                    m |= 1u64 << i;
                }
            }
            m |= 1u64 << codon_index(0, 3, 2); // ATG 恒为起始
            m
        }
        2 => !stop_mask, // 任意非终止密码子
        _ => return Err(format!("start 模式必须是 0/1/2, 收到 {start_mode}")),
    };
    if let Some(list) = extra_starts {
        start_mask = 0;
        for s in list {
            let i = parse_codon(s)?;
            start_mask |= 1u64 << i;
        }
    }
    if let Some(list) = extra_stops {
        stop_mask = 0;
        for s in list {
            let i = parse_codon(s)?;
            stop_mask |= 1u64 << i;
        }
    }
    Ok(Code {
        aa,
        start_mask,
        stop_mask,
    })
}

fn parse_codon(s: &str) -> Result<usize, String> {
    let s = s.trim().to_ascii_uppercase();
    let b = s.as_bytes();
    if b.len() != 3 {
        return Err(format!("密码子长度应为 3: '{s}'"));
    }
    let enc = |c: u8| match c {
        b'A' => Ok(0u8),
        b'C' => Ok(1),
        b'G' => Ok(2),
        b'T' | b'U' => Ok(3),
        _ => Err(format!("密码子含非 ACGT(U) 字符: '{s}'")),
    };
    let (e1, e2, e3) = (enc(b[0])?, enc(b[1])?, enc(b[2])?);
    Ok(codon_index(e1, e2, e3))
}
