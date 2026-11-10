// fastorf: 多线程 ORF 查找 (NCBI ORFfinder 的重写)。用法见 README。
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod fasta;
mod genetic_code;
mod genetic_code_data;
mod orf;
mod output;

use orf::OrfParams;
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::process::ExitCode;

const USAGE: &str = "\
fastorf - 多线程开放阅读框 (ORF) 查找

用法: fastorf --in <FASTA> [选项]

输入:
  --in FILE             输入 FASTA (支持 .gz, '-' 为 stdin)
  --begin N             片段起始 (1-based, 默认 1)
  --end N               片段终止 (0 = 到末尾, 默认 0)
  --circular            环状序列, 允许跨原点 ORF

搜索:
  --gc N                遗传密码表 1-6,9-16,21-33 (默认 1)
  --start 0|1|2         0=仅 ATG; 1=ATG+该表替代起始密码子; 2=任意 sense (默认 1)
  --start-codons LIST   显式起始密码子, 逗号分隔 (覆盖 --start)
  --stop-codons LIST    显式终止密码子, 逗号分隔 (覆盖表内终止)
  --ml N                最小 ORF 长度 nt, 含终止密码子 (默认 75; <30 警告但照用)
  --nested all|drop     drop = 删除完全嵌套于另一 ORF 内的 (默认 all)
  --strand both|plus|minus
  --frames LIST         指定框, 如 1,3,4-6 (1-3 正向, 4-6 反向; 覆盖 --strand)

输出:
  --out FILE            主输出 (默认 stdout)
  --outfmt 0|1|2|3      0=ORF 蛋白 FASTA(默认) 1=CDS FASTA 2=Text ASN.1 3=feature table
  --faa FILE            蛋白 FASTA 独立输出 (与 --outfmt 1 的 --out 等价, 可同轮组合)
  --fnt FILE            CDS 核酸 FASTA 独立输出
  --tsv FILE            ID/坐标映射表 (TSV)
  --gff3 FILE           GFF3
  --id-mode numeric|keep  numeric=重编号(默认); keep=清洗保留原 ID
  --prefix STR          numeric 模式前缀 (默认 SEQ)

其他:
  --threads N           并行线程数 (默认全部核)
  -V/--version, -h/--help";

struct Args {
    input: Option<String>,
    out: Option<String>,
    outfmt: u8,
    faa: Option<String>,
    fnt: Option<String>,
    tsv: Option<String>,
    gff3: Option<String>,
    gc: u32,
    start_mode: u8,
    start_codons: Option<Vec<String>>,
    stop_codons: Option<Vec<String>>,
    ml: usize,
    nested_drop: bool,
    strand: String,
    frames: Option<Vec<usize>>,
    begin: usize,
    end: usize,
    circular: bool,
    id_mode_keep: bool,
    prefix: String,
    threads: Option<usize>,
}

impl Default for Args {
    fn default() -> Self {
        Args {
            input: None,
            out: None,
            outfmt: 0,
            faa: None,
            fnt: None,
            tsv: None,
            gff3: None,
            gc: 1,
            start_mode: 1,
            start_codons: None,
            stop_codons: None,
            ml: 75,
            nested_drop: false,
            strand: "both".into(),
            frames: None,
            begin: 1,
            end: 0,
            circular: false,
            id_mode_keep: false,
            prefix: "SEQ".into(),
            threads: None,
        }
    }
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args::default();
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    let need = |name: &str, i: &mut usize, argv: &Vec<String>| -> Result<String, String> {
        *i += 1;
        argv.get(*i)
            .cloned()
            .ok_or_else(|| format!("参数 {name} 缺少值"))
    };
    while i < argv.len() {
        let arg = argv[i].clone();
        match arg.as_str() {
            "-h" | "--help" => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            "-V" | "--version" => {
                println!("fastorf {}", env!("CARGO_PKG_VERSION"));
                std::process::exit(0);
            }
            "--in" => a.input = Some(need(&arg, &mut i, &argv)?),
            "--out" => a.out = Some(need(&arg, &mut i, &argv)?),
            "--outfmt" => a.outfmt = need(&arg, &mut i, &argv)?.parse().map_err(|_| "--outfmt 需为 0-3".to_string())?,
            "--faa" => a.faa = Some(need(&arg, &mut i, &argv)?),
            "--fnt" => a.fnt = Some(need(&arg, &mut i, &argv)?),
            "--tsv" => a.tsv = Some(need(&arg, &mut i, &argv)?),
            "--gff3" => a.gff3 = Some(need(&arg, &mut i, &argv)?),
            "--gc" => a.gc = need(&arg, &mut i, &argv)?.parse().map_err(|_| "--gc 需为整数".to_string())?,
            "--start" => a.start_mode = need(&arg, &mut i, &argv)?.parse().map_err(|_| "--start 需为 0/1/2".to_string())?,
            "--start-codons" => a.start_codons = Some(need(&arg, &mut i, &argv)?.split(',').map(String::from).collect()),
            "--stop-codons" => a.stop_codons = Some(need(&arg, &mut i, &argv)?.split(',').map(String::from).collect()),
            "--ml" => a.ml = need(&arg, &mut i, &argv)?.parse().map_err(|_| "--ml 需为非负整数".to_string())?,
            "--nested" => {
                let v = need(&arg, &mut i, &argv)?;
                a.nested_drop = match v.as_str() {
                    "drop" => true,
                    "all" => false,
                    _ => return Err("--nested 需为 all|drop".into()),
                };
            }
            "--strand" => {
                a.strand = need(&arg, &mut i, &argv)?;
                if !matches!(a.strand.as_str(), "both" | "plus" | "minus") {
                    return Err("--strand 需为 both|plus|minus".into());
                }
            }
            "--frames" => a.frames = Some(parse_frames(&need(&arg, &mut i, &argv)?)?),
            "--begin" => a.begin = need(&arg, &mut i, &argv)?.parse().map_err(|_| "--begin 需为正整数".to_string())?,
            "--end" => a.end = need(&arg, &mut i, &argv)?.parse().map_err(|_| "--end 需为非负整数".to_string())?,
            "--circular" => a.circular = true,
            "--id-mode" => {
                let v = need(&arg, &mut i, &argv)?;
                a.id_mode_keep = match v.as_str() {
                    "keep" => true,
                    "numeric" => false,
                    _ => return Err("--id-mode 需为 numeric|keep".into()),
                };
            }
            "--prefix" => a.prefix = need(&arg, &mut i, &argv)?,
            "--threads" => a.threads = Some(need(&arg, &mut i, &argv)?.parse().map_err(|_| "--threads 需为正整数".to_string())?),
            _ => return Err(format!("未知参数: {arg} (--help 查看用法)")),
        }
        i += 1;
    }
    if a.input.is_none() {
        return Err("必须提供 --in <FASTA>".into());
    }
    if a.outfmt > 3 {
        return Err("--outfmt 需为 0-3".into());
    }
    Ok(a)
}

// "1,3,4-6" -> [0,2,3,4,5] (0-based 框下标)
fn parse_frames(s: &str) -> Result<Vec<usize>, String> {
    let mut out = Vec::new();
    for part in s.split(',') {
        let part = part.trim();
        let range: Vec<&str> = part.splitn(2, '-').collect();
        let (lo, hi) = match range.len() {
            1 => (part.parse::<usize>(), part.parse::<usize>()),
            2 => (
                range[0].trim().parse::<usize>(),
                range[1].trim().parse::<usize>(),
            ),
            _ => unreachable!(),
        };
        let (lo, hi) = (lo.map_err(|_| format!("--frames 非法: {part}"))?, hi.map_err(|_| format!("--frames 非法: {part}"))?);
        for f in lo..=hi {
            if !(1..=6).contains(&f) {
                return Err(format!("--frames 框号需在 1-6: {f}"));
            }
            out.push(f - 1);
        }
    }
    if out.is_empty() {
        return Err("--frames 为空".into());
    }
    Ok(out)
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("fastorf: 错误: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), String> {
    let a = parse_args()?;
    if a.ml < 30 {
        eprintln!("fastorf: 警告: --ml {} 低于官方下限 30, 按给定值执行 (官方会静默改为 30)", a.ml);
    }
    if a.begin < 1 {
        return Err("--begin 需 >= 1".into());
    }
    let code = genetic_code::build(
        a.gc,
        a.start_mode,
        a.start_codons.as_deref(),
        a.stop_codons.as_deref(),
    )?;

    let recs = fasta::read_fasta(a.input.as_deref().unwrap())?;

    // 串行分配规范化 ID, 保证全局唯一且确定
    let n = recs.len();
    let seq_ids: Vec<String> = if a.id_mode_keep {
        let mut used = std::collections::HashSet::new();
        recs.iter()
            .enumerate()
            .map(|(i, r)| {
                let fallback = format!("{}{:06}", a.prefix, i + 1);
                let mut s = fasta::sanitize_id(&r.id, &fallback);
                if !used.insert(s.clone()) {
                    let base = s;
                    let mut k = 2;
                    loop {
                        s = format!("{base}_{k}");
                        if used.insert(s.clone()) {
                            break;
                        }
                        k += 1;
                    }
                }
                s
            })
            .collect()
    } else {
        (0..n).map(|i| format!("{}{:06}", a.prefix, i + 1)).collect()
    };

    let mut frames = [false; 6];
    match &a.frames {
        Some(list) => {
            for &f in list {
                frames[f] = true;
            }
        }
        None => match a.strand.as_str() {
            "plus" => {
                frames[0] = true;
                frames[1] = true;
                frames[2] = true;
            }
            "minus" => {
                frames[3] = true;
                frames[4] = true;
                frames[5] = true;
            }
            _ => frames = [true; 6],
        },
    }

    if let Some(t) = a.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .map_err(|e| format!("线程池初始化失败: {e}"))?;
    }

    // 每条序列: 并行完成扫描与全部格式化, 串行只拼接字节 -> 结果与线程数无关
    let want_tsv = a.tsv.is_some();
    let want_gff3 = a.gff3.is_some();
    let want_faa = a.faa.is_some();
    let want_fnt = a.fnt.is_some();
    let outfmt = a.outfmt;
    let results: Vec<SeqOutput> = recs
        .par_iter()
        .enumerate()
        .map(|(i, r)| {
            let frag = fragment(&r.seq, a.begin, a.end)?;
            let orfs = orf::find_orfs(
                &frag,
                &code,
                &OrfParams {
                    min_len: a.ml,
                    nested_drop: a.nested_drop,
                    frames,
                    circular: a.circular,
                    seq_len: frag.len(),
                },
            );
            // 报告坐标 = 全长序列坐标 (片段偏移加回, 口径 M4 对拍校准)
            let shift = a.begin - 1;
            let mut orfs = orfs;
            if shift > 0 {
                for o in orfs.iter_mut() {
                    o.start += shift;
                    o.end += shift;
                }
            }
            let sv = output::SeqView {
                seq_id: &seq_ids[i],
                orig_id: &r.id,
                desc: &r.desc,
                orfs: &orfs,
                code_id: a.gc,
                code: &code,
                seq: &r.seq, // 坐标基于全长, 序列也用全长
            };
            let mut main_buf = Vec::new();
            match outfmt {
                0 => output::write_orf_fasta(&mut main_buf, &sv),
                1 => output::write_cds_fasta(&mut main_buf, &sv),
                2 => output::write_asn1(&mut main_buf, &sv),
                _ => output::write_ftable(&mut main_buf, &sv),
            }
            .map_err(|e| format!("格式化输出失败: {e}"))?;
            let faa_buf = if want_faa {
                let mut b = Vec::new();
                output::write_orf_fasta(&mut b, &sv).map_err(|e| e.to_string())?;
                b
            } else {
                Vec::new()
            };
            let fnt_buf = if want_fnt {
                let mut b = Vec::new();
                output::write_cds_fasta(&mut b, &sv).map_err(|e| e.to_string())?;
                b
            } else {
                Vec::new()
            };
            let tsv_buf = if want_tsv {
                let mut b = Vec::new();
                output::write_tsv(&mut b, &sv).map_err(|e| e.to_string())?;
                b
            } else {
                Vec::new()
            };
            let gff3_buf = if want_gff3 {
                let mut b = Vec::new();
                output::write_gff3(&mut b, &sv).map_err(|e| e.to_string())?;
                b
            } else {
                Vec::new()
            };
            Ok(SeqOutput {
                main: main_buf,
                faa: faa_buf,
                fnt: fnt_buf,
                tsv: tsv_buf,
                gff3: gff3_buf,
                orf_count: orfs.len(),
            })
        })
        .collect::<Result<Vec<_>, String>>()?;

    write_outputs(&a, &results)?;
    let total: usize = results.iter().map(|r| r.orf_count).sum();
    let gc_name = genetic_code::lookup_table(a.gc)
        .map(|t| t.name)
        .unwrap_or("?");
    eprintln!(
        "fastorf: {} 条序列, {} 个 ORF (gc={} [{}], start={}, ml={}, nested={}, strand={})",
        recs.len(), total, a.gc, gc_name, a.start_mode, a.ml,
        if a.nested_drop { "drop" } else { "all" }, a.strand
    );
    Ok(())
}

fn fragment(seq: &[u8], begin: usize, end: usize) -> Result<Vec<u8>, String> {
    let b = begin - 1;
    if b >= seq.len() {
        return Err(format!("--begin {begin} 超出序列长度 {}", seq.len()));
    }
    let e = if end == 0 || end > seq.len() { seq.len() } else { end };
    if e < b {
        return Err("--end 小于 --begin".into());
    }
    Ok(seq[b..e].to_vec())
}

// 每条序列的最终输出字节 (并行阶段生成完毕)
struct SeqOutput {
    main: Vec<u8>,
    faa: Vec<u8>,
    fnt: Vec<u8>,
    tsv: Vec<u8>,
    gff3: Vec<u8>,
    orf_count: usize,
}

fn write_outputs(a: &Args, results: &[SeqOutput]) -> Result<(), String> {
    // 多文件输出时并行写 (文件间无顺序依赖, 文件内保持输入顺序)
    let mut errs: Vec<String> = Vec::new();
    std::thread::scope(|s| {
        let mut handles = Vec::new();
        if let Some(p) = &a.out {
            handles.push(s.spawn(move || write_main_file(p, results)));
        }
        if let Some(p) = &a.faa {
            handles.push(s.spawn(move || write_plain(p, &results, |r| &r.faa)));
        }
        if let Some(p) = &a.fnt {
            handles.push(s.spawn(move || write_plain(p, &results, |r| &r.fnt)));
        }
        if let Some(p) = &a.tsv {
            handles.push(s.spawn(move || write_tsv_file(p, a, results)));
        }
        if let Some(p) = &a.gff3 {
            handles.push(s.spawn(move || write_gff3_file(p, results)));
        }
        for h in handles {
            match h.join() {
                Ok(Ok(())) => {}
                Ok(Err(e)) => errs.push(e),
                Err(_) => errs.push("输出线程 panic".into()),
            }
        }
    });
    if a.out.is_none() {
        let mut w = BufWriter::with_capacity(4 << 20, std::io::stdout().lock());
        for r in results {
            w.write_all(&r.main).map_err(|e| format!("写 stdout 失败: {e}"))?;
        }
        w.flush().map_err(|e| format!("flush stdout 失败: {e}"))?;
    }
    if errs.is_empty() {
        Ok(())
    } else {
        Err(errs.join("; "))
    }
}

fn write_main_file(p: &str, results: &[SeqOutput]) -> Result<(), String> {
    write_plain(p, results, |r| &r.main)
}

fn write_plain(p: &str, results: &[SeqOutput], pick: fn(&SeqOutput) -> &[u8]) -> Result<(), String> {
    let mut w = file_writer(p)?;
    for r in results {
        w.write_all(pick(r)).map_err(|e| format!("写 {p} 失败: {e}"))?;
    }
    w.flush().map_err(|e| format!("flush {p} 失败: {e}"))
}

fn write_tsv_file(p: &str, a: &Args, results: &[SeqOutput]) -> Result<(), String> {
    let mut w = file_writer(p)?;
    writeln!(&mut w, "# fastorf {}\tgc={}\tstart={}\tml={}\tnested={}\tcircular={}",
        env!("CARGO_PKG_VERSION"), a.gc, a.start_mode, a.ml,
        if a.nested_drop { "drop" } else { "all" }, a.circular)
        .map_err(|e| e.to_string())?;
    writeln!(&mut w, "orf_id\torig_id\tdesc\tstrand\tframe\tstart\tend\tnt_len\taa_len\tstart_codon\tstop_codon\tstatus")
        .map_err(|e| e.to_string())?;
    for r in results {
        w.write_all(&r.tsv).map_err(|e| format!("写 {p} 失败: {e}"))?;
    }
    w.flush().map_err(|e| format!("flush {p} 失败: {e}"))
}

fn write_gff3_file(p: &str, results: &[SeqOutput]) -> Result<(), String> {
    let mut w = file_writer(p)?;
    writeln!(&mut w, "##gff-version 3").map_err(|e| e.to_string())?;
    for r in results {
        w.write_all(&r.gff3).map_err(|e| format!("写 {p} 失败: {e}"))?;
    }
    w.flush().map_err(|e| format!("flush {p} 失败: {e}"))
}
fn file_writer(p: &str) -> Result<Box<dyn Write>, String> {
    let f = std::fs::File::create(p).map_err(|e| format!("创建 {p} 失败: {e}"))?;
    Ok(Box::new(BufWriter::with_capacity(4 << 20, f)))
}
