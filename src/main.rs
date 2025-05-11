// src/main.rs
use anyhow::{Context, Result};
use clap::{ArgAction, Parser};
use csv::ReaderBuilder;
use log::{debug, info};
use rust_htslib::bcf::{self, header::HeaderRecord, Header, Read, Record, Writer};
use rust_htslib::bcf::record::GenotypeAllele;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{self, Read as _};

/// Command‑line options ------------------------------------------------------
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Opts {
    /// VCF input (default: stdin)
    #[arg(value_name = "VCF", default_value = "-")]
    input: String,

    /// Output VCF (default: stdout)
    #[arg(short, long, default_value = "-")]
    output: String,

    /// Tab‑delimited 2‑col file: <sample> <group>
    #[arg(short, long)]
    labels: String,

    /// Comma‑sep tags to add, or ‘all’
    #[arg(short, long, default_value = "all")]
    tags: String,

    /// Abort if labels contain samples absent from VCF
    #[arg(long, action = ArgAction::SetTrue)]
    strict: bool,

    /// Verbose logging
    #[arg(long, action = ArgAction::SetTrue)]
    debug: bool,
}

/// Per‑variant statistics ----------------------------------------------------
#[derive(Default, Debug, Clone)]
struct AfStats {
    ac: [u32; 2],     // allele counts (REF, ALT)
    an: u32,          // allele number
    n_hemi: u32,
    n_homref: u32,
    n_het: u32,
    n_homalt: u32,
    n_miss: u32,
    af: f64,
    maf: f64,
    mac: u32,
    exc_het: f64,
    hwe: f64,
}

/// Hardy‑Weinberg exact test + excess‑het flag (简化版，与 truvari 近似)
fn calc_hwe(ac0: u32, ac1: u32, n_het: u32) -> (f64, f64) {
    // 参考 Wigginton et al. 2005 (PMID:15789306)  
    // 这里给出一个数值稳定的近似实现
    let n_homref = (ac0 - n_het) / 2;
    let n_homalt = (ac1 - n_het) / 2;
    let n = n_homref + n_het + n_homalt;

    if n == 0 {
        return (0.0, 0.0);
    }

    let exp_het = (2.0 * (ac0 as f64) * (ac1 as f64)) / ((2 * n) as f64);
    let chi_sq = ((n_het as f64 - exp_het).powi(2)) / exp_het.max(1.0);
    let p = (-0.5 * chi_sq).exp();      // 近似 p‑value
    let exc_het = if p < 1e-6 { 0.0 } else { 1.0 }; // 1=good, 0=bad
    (exc_het, p)
}

/// 计算所有统计量（单群体，等价于 Python 的 calc_af）
fn calc_af(genotypes: &[Option<[Option<u8>; 2]>]) -> AfStats {
    let mut st = AfStats::default();

    for g in genotypes {
        match g {
            None => st.n_miss += 1,
            Some(alleles) => {
                // haploid?
                if alleles[1].is_none() {
                    match alleles[0] {
                        None => st.n_miss += 1,
                        Some(a) => {
                            st.an += 1;
                            st.ac[a as usize] += 1;
                            st.n_hemi += 1;
                        }
                    }
                } else {
                    // diploid
                    let a0 = alleles[0];
                    let a1 = alleles[1];
                    for a in [a0, a1] {
                        if let Some(x) = a {
                            st.an += 1;
                            st.ac[x as usize] += 1;
                        }
                    }
                    match (a0, a1) {
                        (None, None) => st.n_miss += 1,
                        (Some(x), Some(y)) if x == y && x == 0 => st.n_homref += 1,
                        (Some(x), Some(y)) if x == y && x == 1 => st.n_homalt += 1,
                        (Some(_), Some(_)) => st.n_het += 1,
                        _ => st.n_hemi += 1, // one missing
                    }
                }
            }
        }
    }

    if st.an > 0 {
        st.af = st.ac[1] as f64 / st.an as f64;
        st.mac = st.ac[0].min(st.ac[1]);
        st.maf = st.mac as f64 / st.an as f64;
        let (exc_het, hwe) = calc_hwe(st.ac[0], st.ac[1], st.n_het);
        st.exc_het = exc_het;
        st.hwe = hwe;
    }
    st
}

/// ------------- 主程序 ------------------------------------------------------
fn main() -> Result<()> {
    let opts = Opts::parse();

    // logging
    if opts.debug {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    }

    // ------------------- 读取 label 文件 ------------------------------------
    #[derive(serde::Deserialize)]
    struct Rec {
        sample: String,
        group: String,
    }
    let mut rdr = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(&opts.labels)
        .with_context(|| "opening --labels")?;
    let mut group_map: HashMap<String, Vec<String>> = HashMap::new();
    for rec in rdr.deserialize::<Rec>() {
        let rec = rec?;
        group_map.entry(rec.group).or_default().push(rec.sample);
    }
    let groups: Vec<String> = group_map.keys().cloned().collect();
    info!("Loaded {} groups from {}", groups.len(), opts.labels);

    // ------------------- 打开 VCF -------------------------------------------
    let mut reader: Box<dyn Read> = if opts.input == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(&opts.input).with_context(|| "opening VCF")?)
    };
    let mut bcf = bcf::Reader::from_stream(&mut reader).with_context(|| "reading VCF")?;
    let samples: Vec<String> = bcf.header().samples().iter().map(|s| String::from_utf8_lossy(s).into()).collect();

    // sanity check
    if opts.strict {
        for (s, _) in group_map.values().flat_map(|v| v.iter()) {
            if !samples.contains(s) {
                anyhow::bail!("Sample `{s}` in --labels but not in VCF header (strict mode)");
            }
        }
    }

    // -------------- 构建每组 sample mask ------------------------------------
    let mut masks: HashMap<String, Vec<bool>> = HashMap::new();
    for (grp, list) in &group_map {
        let set: HashSet<&String> = list.iter().collect();
        masks.insert(
            grp.clone(),
            samples.iter().map(|s| set.contains(s)).collect(),
        );
    }

    // ----------------- 克隆并扩充 header ------------------------------------
    let mut out_hdr = Header::from_template(bcf.header());
    let all_tags = [
        "AF", "MAF", "ExcHet", "HWE", "MAC", "AC", "AN", "N_HEMI", "N_MISS",
        "N_HOMREF", "N_HET", "N_HOMALT",
    ];
    let want_tags: Vec<&str> = if opts.tags == "all" {
        all_tags.to_vec()
    } else {
        opts.tags.split(',').collect()
    };

    let mut add_info_line = |id: &str, num: &str, typ: &str, desc: &str| {
        let line = format!("##INFO=<ID={id},Number={num},Type={typ},Description=\"{desc}\">");
        out_hdr.push_record(HeaderRecord::from_str(&line));
    };

    for grp in &groups {
        let count = group_map[grp].len();
        for t in &want_tags {
            match *t {
                "AC" | "MAC" => add_info_line(
                    &format!("{t}_{grp}"),
                    "A",
                    if *t == "AC" || *t == "MAC" { "Integer" } else { "Float" },
                    &format!("{t} on {count} {grp} samples"),
                ),
                _ => add_info_line(
                    &format!("{t}_{grp}"),
                    "1",
                    if ["AF", "MAF", "ExcHet", "HWE"].contains(t) {
                        "Float"
                    } else {
                        "Integer"
                    },
                    &format!("{t} on {count} {grp} samples"),
                ),
            }
        }
    }

    // ----------------- 打开输出 ---------------------------------------------
    let mut writer = if opts.output == "-" {
        Writer::from_stdout(&out_hdr, true, bcf::Format::VCF)?
    } else {
        Writer::from_path(&opts.output, &out_hdr, true, bcf::Format::VCF)?
    };

    // ----------------- 逐 variant 处理 --------------------------------------
    for rec_result in bcf.records() {
        let mut rec = rec_result?;
        let gts = rec.genotypes()?;

        // 为每个 group 收集其 sample genotype
        for grp in &groups {
            let mask = &masks[grp];
            let mut gt_vec: Vec<Option<[Option<u8>; 2]>> = Vec::new();
            for (samp_idx, alleles) in gts.iter().enumerate() {
                if !mask[samp_idx] {
                    continue;
                }
                // extract (max 2) alleles
                let mut pair = [None, None];
                for (i, a) in alleles.iter().take(2).enumerate() {
                    pair[i] = match a.index() {
                        Some(idx) if idx >= 0 => Some(idx as u8),
                        _ => None,
                    };
                }
                if pair.iter().all(|x| x.is_none()) {
                    gt_vec.push(None);
                } else {
                    gt_vec.push(Some(pair));
                }
            }
            let stats = calc_af(&gt_vec);

            // ALT allele 仅写 AC/MAC 的 ALT 值，与 Python 逻辑保持一致
            for tag in &want_tags {
                let full = format!("{tag}_{grp}");
                match *tag {
                    "AC" => rec.push_info_integer(&full, &[stats.ac[1] as i32])?,
                    "MAC" => rec.push_info_integer(&full, &[stats.mac as i32])?,
                    "AN" => rec.push_info_integer(&full, &[stats.an as i32])?,
                    "N_HEMI" => rec.push_info_integer(&full, &[stats.n_hemi as i32])?,
                    "N_MISS" => rec.push_info_integer(&full, &[stats.n_miss as i32])?,
                    "N_HOMREF" => rec.push_info_integer(&full, &[stats.n_homref as i32])?,
                    "N_HET" => rec.push_info_integer(&full, &[stats.n_het as i32])?,
                    "N_HOMALT" => rec.push_info_integer(&full, &[stats.n_homalt as i32])?,
                    "AF" => rec.push_info_float(&full, &[stats.af as f32])?,
                    "MAF" => rec.push_info_float(&full, &[stats.maf as f32])?,
                    "ExcHet" => rec.push_info_float(&full, &[stats.exc_het as f32])?,
                    "HWE" => rec.push_info_float(&full, &[stats.hwe as f32])?,
                    _ => {}
                }
            }
        }
        writer.write(&rec)?;
    }

    info!("Finished grpaf‑rs");
    Ok(())
}
