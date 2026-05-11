use crate::args::CompareArgs;
use anyhow::{Context, Result};
use rinfuse_io::fc_intermediates::{FcCandidateWithSource, FcOutputDir};
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct NormalizedPair(String, String);

impl NormalizedPair {
    fn new(a: &str, b: &str) -> Self {
        let mut genes = [a.to_string(), b.to_string()];
        genes.sort();
        Self(genes[0].clone(), genes[1].clone())
    }
}

pub fn run(args: CompareArgs) -> Result<()> {
    let fc_candidates = load_candidates(&args.fc)?;
    let rs_candidates = load_candidates(&args.rs)?;

    let mut fc_map = HashMap::new();
    for cand in fc_candidates {
        let pair = NormalizedPair::new(&cand.gene_5p, &cand.gene_3p);
        fc_map.insert(pair, cand);
    }

    let mut rs_map = HashMap::new();
    for cand in rs_candidates {
        let pair = NormalizedPair::new(&cand.gene_5p, &cand.gene_3p);
        rs_map.insert(pair, cand);
    }

    let mut all_pairs: HashSet<NormalizedPair> = fc_map.keys().cloned().collect();
    all_pairs.extend(rs_map.keys().cloned());

    let mut sorted_pairs: Vec<NormalizedPair> = all_pairs.into_iter().collect();
    sorted_pairs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    // 1. Write compare.tsv
    let file = fs::File::create(&args.out)
        .with_context(|| format!("failed to create output file {}", args.out.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "gene_a\tgene_b\tstatus\tfc_spanning\trs_spanning\tfc_split\trs_split\tfc_source\trs_source"
    )?;

    let mut shared = 0;
    let mut only_fc = 0;
    let mut only_rs = 0;

    for pair in &sorted_pairs {
        let fc_opt = fc_map.get(pair);
        let rs_opt = rs_map.get(pair);

        let status = match (fc_opt.is_some(), rs_opt.is_some()) {
            (true, true) => {
                shared += 1;
                "both"
            }
            (true, false) => {
                only_fc += 1;
                "only_fc"
            }
            (false, true) => {
                only_rs += 1;
                "only_rs"
            }
            _ => unreachable!(),
        };

        let fc_spanning = fc_opt
            .and_then(|c| c.spanning_pairs.map(|v| v.to_string()))
            .unwrap_or_else(|| "-".to_string());
        let rs_spanning = rs_opt
            .and_then(|c| c.spanning_pairs.map(|v| v.to_string()))
            .unwrap_or_else(|| "-".to_string());

        let fc_split = fc_opt
            .and_then(|c| c.split_reads.map(|v| v.to_string()))
            .unwrap_or_else(|| "-".to_string());
        let rs_split = rs_opt
            .and_then(|c| c.split_reads.map(|v| v.to_string()))
            .unwrap_or_else(|| "-".to_string());

        let fc_source = fc_opt
            .map(|c| c.source.display().to_string())
            .unwrap_or_else(|| "-".to_string());
        let rs_source = rs_opt
            .map(|c| c.source.display().to_string())
            .unwrap_or_else(|| "-".to_string());

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            pair.0,
            pair.1,
            status,
            fc_spanning,
            rs_spanning,
            fc_split,
            rs_split,
            fc_source,
            rs_source
        )?;
    }
    writer.flush()?;

    // 2. Write compare_summary.md
    let summary_path = args.out.with_extension("md");
    if let Some(parent) = summary_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    let mut sw = BufWriter::new(fs::File::create(&summary_path)?);

    writeln!(sw, "# rinfuse-fc Comparison Summary\n")?;
    writeln!(sw, "## Stats")?;
    writeln!(sw, "- **Total FC Candidates**: {}", fc_map.len())?;
    writeln!(sw, "- **Total RS Candidates**: {}", rs_map.len())?;
    writeln!(sw, "- **Shared**: {}", shared)?;
    writeln!(sw, "- **Only FC**: {}", only_fc)?;
    writeln!(sw, "- **Only RS**: {}", only_rs)?;

    let focus_genes: HashSet<String> = args
        .focus_gene
        .iter()
        .flat_map(|s| s.split(',').map(|g| g.trim().to_string()))
        .filter(|s| !s.is_empty())
        .collect();

    if !focus_genes.is_empty() {
        writeln!(sw, "\n## Focus Genes ({})", args.focus_gene.join(", "))?;
        writeln!(
            sw,
            "| Gene A | Gene B | Status | FC Spanning | RS Spanning |"
        )?;
        writeln!(sw, "|---|---|---|---|---|")?;
        for pair in &sorted_pairs {
            if focus_genes.contains(&pair.0) || focus_genes.contains(&pair.1) {
                let fc_opt = fc_map.get(pair);
                let rs_opt = rs_map.get(pair);
                let status = match (fc_opt.is_some(), rs_opt.is_some()) {
                    (true, true) => "both",
                    (true, false) => "only_fc",
                    (false, true) => "only_rs",
                    _ => unreachable!(),
                };
                let fc_spanning = fc_opt
                    .and_then(|c| c.spanning_pairs.map(|v| v.to_string()))
                    .unwrap_or_else(|| "-".to_string());
                let rs_spanning = rs_opt
                    .and_then(|c| c.spanning_pairs.map(|v| v.to_string()))
                    .unwrap_or_else(|| "-".to_string());
                writeln!(
                    sw,
                    "| {} | {} | {} | {} | {} |",
                    pair.0, pair.1, status, fc_spanning, rs_spanning
                )?;
            }
        }
    }

    sw.flush()?;

    eprintln!(
        "compare complete: fc={} rs={} shared={} only_fc={} only_rs={}",
        fc_map.len(),
        rs_map.len(),
        shared,
        only_fc,
        only_rs
    );

    Ok(())
}

fn load_candidates(path: &Path) -> Result<Vec<FcCandidateWithSource>> {
    if path.is_file() {
        let out = FcOutputDir {
            root: path
                .parent()
                .unwrap_or_else(|| Path::new("."))
                .to_path_buf(),
            ..Default::default()
        };
        out.parse_candidate_file(path)
    } else if path.is_dir() {
        let fc_dir = FcOutputDir::detect(path, true, 3)?;
        fc_dir.parse_all_candidates()
    } else {
        anyhow::bail!("Path does not exist: {}", path.display())
    }
}

// I need to add parse_all_candidates_from_file to FcOutputDir in fc_intermediates.rs
// to make this clean.
