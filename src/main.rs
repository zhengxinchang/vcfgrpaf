// src/main.rs
use anyhow::{Context, Result};
use clap::Parser;
use csv::ReaderBuilder;
use rust_htslib::{
    bcf::{self, Header, HeaderRecord, Read, Reader as BcfReader, Writer},
};

use std::{

    collections::{HashMap, HashSet}
};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "This program is the rust re-implimentation of grpaf.py from truvari"
)]
struct Opts {
    /// VCF input
    #[arg(value_name = "VCF")]
    input: String,

    /// Output VCF
    #[arg(short, long)]
    output: String,

    /// Tab‑delimited 2‑col file: <sample> <group>
    #[arg(short, long)]
    labels: String,
}

/// Per‑variant statistics ----------------------------------------------------
#[derive(Default, Debug, Clone)]
struct AfStats {
    ac: [u32; 2], // allele counts (REF, ALT)
    an: u32,      // allele number
    n_hemi: u32,
    n_homref: u32,
    n_het: u32,
    n_homalt: u32,
    n_miss: u32,
    af: f64,
    maf: f64,
    mac: u32,
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
    }
    st
}

fn main() -> Result<()> {
    let opts = Opts::parse();

    // read labels
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
    println!("Loaded {} groups from {}", groups.len(), opts.labels);

    // open input VCF

    let mut bcf: BcfReader = BcfReader::from_path(&opts.input).expect("Error opening file.");
    let mut headerview: bcf::header::HeaderView = bcf.header().clone();

    let samples: Vec<String> = headerview
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).into())
        .collect();

    // make masks for each group
    let mut masks: HashMap<String, Vec<bool>> = HashMap::new();
    for (grp, list) in &group_map {
        let set: HashSet<&String> = list.iter().collect();
        masks.insert(
            grp.clone(),
            samples.iter().map(|s| set.contains(s)).collect(),
        );
    }

    // inject new header lines

    let want_tags = [
        "AF", "MAF", "MAC", "AC", "AN", "N_HEMI", "N_MISS", "N_HOMREF", "N_HET", "N_HOMALT",
    ];
    let all_tags = [
        "ExcHet_",
        "HWE_",
        "AF_",
        "MAF_",
        "MAC_",
        "AC_",
        "AN_",
        "N_HEMI_",
        "N_MISS_",
        "N_HOMREF_",
        "N_HET_",
        "N_HOMALT_",
    ];

    let mut out_hdr = Header::from_template(bcf.header());
    

    let mut add_info_line = |id: &str, num: &str, typ: &str, desc: &str| {
        let line = format!("##INFO=<ID={id},Number={num},Type={typ},Description=\"{desc}\">");
        out_hdr.push_record(line.as_bytes());
    };

    for grp in &groups {
        let count = group_map[grp].len();
        for t in &want_tags {
            match *t {
                "AC" | "MAC" => add_info_line(
                    &format!("{t}_{grp}"),
                    "A",
                    if *t == "AC" || *t == "MAC" {
                        "Integer"
                    } else {
                        "Float"
                    },
                    &format!("{t} on {count} {grp} samples"),
                ),
                _ => add_info_line(
                    &format!("{t}_{grp}"),
                    "1",
                    if ["AF", "MAF"].contains(t) {
                        "Float"
                    } else {
                        "Integer"
                    },
                    &format!("{t} on {count} {grp} samples"),
                ),
            }
        }
    }

    // generates all tags that start with all_tags in to a vec
    let mut all_tags_combination: Vec<String> = Vec::new();
    for (_, values) in headerview
        .header_records()
        .iter()
        .filter_map(|record| match record {
            HeaderRecord::Info { key, values } => Some((key, values)),
            _ => None,
        })
    {
        let empty = String::from("NOT_FOUND");
        let id = values.get("ID").unwrap_or(&empty);
        if all_tags.iter().any(|x| {
            id.starts_with(x)
        }) {
            all_tags_combination.push(id.to_string());
        }
    }

    eprintln!("Found related tags in input VCF:  {:?}", all_tags_combination);

    // open output file
    let mut writer = if opts.output == "-" {
        Writer::from_stdout(&out_hdr, true, bcf::Format::Vcf)?
    } else {
        Writer::from_path(&opts.output, &out_hdr, true, bcf::Format::Vcf)?
    };

    // process each record
    let mut processed = 0;
    for rec_result in bcf.records() {
        processed += 1;
        if processed % 10000 == 0 {
            println!("Processed {processed} variants");
        }
        let mut rec = rec_result?;

        // remove all_tags if present

        let gt_vec_map: HashMap<String, Vec<Option<[Option<u8>; 2]>>> = {
            let gts = rec.genotypes()?;
            let mut map = HashMap::new();
            for grp in &groups {
                let mask = &masks[grp];
                let mut gt_vec: Vec<Option<[Option<u8>; 2]>> = Vec::new();
                for samp_idx in 0..headerview.samples().len() {
                    let alleles = gts.get(samp_idx);
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
                map.insert(grp.clone(), gt_vec);
            }
            map
        };

        for tag in all_tags_combination.iter() {
            let full = format!("{tag}");
            rec.push_info_string(&full.as_bytes(), &[])?;
        }

        for grp in &groups {
            let stats = calc_af(&gt_vec_map[grp]);

            for tag in &want_tags {
                let full = format!("{tag}_{grp}");

                rec.push_info_string(&full.as_bytes(), &[])?;

                match *tag {
                    "AC" => rec.push_info_integer(&full.as_bytes(), &[stats.ac[1] as i32])?,
                    "MAC" => rec.push_info_integer(&full.as_bytes(), &[stats.mac as i32])?,
                    "AN" => rec.push_info_integer(&full.as_bytes(), &[stats.an as i32])?,
                    "N_HEMI" => rec.push_info_integer(&full.as_bytes(), &[stats.n_hemi as i32])?,
                    "N_MISS" => rec.push_info_integer(&full.as_bytes(), &[stats.n_miss as i32])?,
                    "N_HOMREF" => {
                        rec.push_info_integer(&full.as_bytes(), &[stats.n_homref as i32])?
                    }
                    "N_HET" => rec.push_info_integer(&full.as_bytes(), &[stats.n_het as i32])?,
                    "N_HOMALT" => {
                        rec.push_info_integer(&full.as_bytes(), &[stats.n_homalt as i32])?
                    }
                    "AF" => rec.push_info_float(&full.as_bytes(), &[stats.af as f32])?,
                    "MAF" => rec.push_info_float(&full.as_bytes(), &[stats.maf as f32])?,
                    _ => {}
                }
            }
        }
        writer.write(&rec)?;
    }

    println!("Finished vcfgrpaf");
    Ok(())
}
