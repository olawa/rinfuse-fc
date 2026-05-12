use crate::args::ValidateSampleArgs;
use anyhow::Result;
use rinfuse_core::FusionCandidateLite;
use rinfuse_io::fc_intermediates::{FcCandidateWithSource, FcOutputDir};
use rinfuse_io::{ReadRegistry, ReadRegistryBuildOptions};
use serde_json::json;
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct NormalizedPair(String, String);

impl NormalizedPair {
    fn new(a: &str, b: &str) -> Self {
        let mut genes = [a.to_uppercase(), b.to_uppercase()];
        genes.sort();
        Self(genes[0].clone(), genes[1].clone())
    }
}

pub fn run(args: ValidateSampleArgs) -> Result<()> {
    if !args.out.exists() {
        fs::create_dir_all(&args.out)?;
    }

    // 1. Load FusionCatcher Candidates
    let fc_dir = FcOutputDir::detect(&args.fc_out, true, 3)?;
    let fc_candidates = fc_dir.parse_all_candidates()?;
    let fc_read_tokens = fc_dir.collect_read_id_tokens()?;

    // 2. Load STAR Candidates
    let star_candidates = load_star_candidates(&args.star_candidates)?;

    // 3. Extract reads if requested
    let mut extracted_reads_count = 0;
    if !args.reads.is_empty() {
        let mut registry = ReadRegistry::from_tokens(&fc_read_tokens);
        let opts = ReadRegistryBuildOptions::default();
        registry.collect_from_fastq_paths(&args.reads, &opts)?;
        registry.write_fastq(&args.out.join("fusioncatcher_supporting_reads.fq"))?;
        extracted_reads_count = registry.found_record_count();
    }

    // 4. Compare
    let mut fc_map: HashMap<NormalizedPair, Vec<&FcCandidateWithSource>> = HashMap::new();
    for cand in &fc_candidates {
        let pair = NormalizedPair::new(&cand.gene_5p, &cand.gene_3p);
        fc_map.entry(pair).or_default().push(cand);
    }

    let mut star_map: HashMap<NormalizedPair, Vec<&FusionCandidateLite>> = HashMap::new();
    for cand in &star_candidates {
        let pair = NormalizedPair::new(&cand.gene_a, &cand.gene_b);
        star_map.entry(pair).or_default().push(cand);
    }

    let mut all_pairs: HashSet<NormalizedPair> = fc_map.keys().cloned().collect();
    all_pairs.extend(star_map.keys().cloned());
    let mut sorted_pairs: Vec<NormalizedPair> = all_pairs.into_iter().collect();
    sorted_pairs.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

    let mut shared_count = 0;
    let mut only_fc_count = 0;
    let mut only_star_count = 0;

    let focus_genes: HashSet<String> = args
        .focus_gene
        .iter()
        .flat_map(|s| s.split(','))
        .map(|s| s.trim().to_uppercase())
        .filter(|s| !s.is_empty())
        .collect();

    let mut missing_from_star = Vec::new();

    // Write TSVs
    let mut w_fc_star = BufWriter::new(fs::File::create(args.out.join("fc_vs_star.tsv"))?);
    let mut w_focus = if !focus_genes.is_empty() {
        Some(BufWriter::new(fs::File::create(
            args.out.join("focus_fc_vs_star.tsv"),
        )?))
    } else {
        None
    };

    let mut w_missing = BufWriter::new(fs::File::create(args.out.join("missing_from_star.tsv"))?);
    let mut w_recovered = BufWriter::new(fs::File::create(args.out.join("recovered_by_star.tsv"))?);

    let header = "gene_a\tgene_b\tstatus\tfc_spanning\tstar_unique_reads\tfc_source";
    writeln!(w_fc_star, "{}", header)?;
    if let Some(w) = &mut w_focus {
        writeln!(w, "{}", header)?;
    }
    writeln!(w_missing, "{}", header)?;
    writeln!(w_recovered, "{}", header)?;

    for pair in &sorted_pairs {
        let has_fc = fc_map.contains_key(pair);
        let has_star = star_map.contains_key(pair);

        let status = match (has_fc, has_star) {
            (true, true) => {
                shared_count += 1;
                "shared"
            }
            (true, false) => {
                only_fc_count += 1;
                "only_fc"
            }
            (false, true) => {
                only_star_count += 1;
                "only_star"
            }
            _ => unreachable!(),
        };

        let is_focus = focus_genes.contains(&pair.0) || focus_genes.contains(&pair.1);
        if is_focus && has_fc && !has_star {
            missing_from_star.push(pair.clone());
        }

        let fc_cand = fc_map.get(pair).and_then(|v| v.first());
        let star_cand = star_map.get(pair).and_then(|v| v.first());

        let fc_spanning = fc_cand
            .and_then(|c| c.spanning_pairs.map(|v| v.to_string()))
            .unwrap_or_else(|| "-".to_string());
        let star_reads = star_cand
            .map(|c| c.unique_read_count.to_string())
            .unwrap_or_else(|| "-".to_string());
        let fc_source = fc_cand
            .map(|c| c.source.display().to_string())
            .unwrap_or_else(|| "-".to_string());

        let row = format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            pair.0, pair.1, status, fc_spanning, star_reads, fc_source
        );

        writeln!(w_fc_star, "{}", row)?;
        if is_focus {
            if let Some(w) = &mut w_focus {
                writeln!(w, "{}", row)?;
            }
        }

        if has_fc && !has_star {
            writeln!(w_missing, "{}", row)?;
        } else if !has_fc && has_star {
            writeln!(w_recovered, "{}", row)?;
        }
    }

    // 5. Write FC Candidates
    let mut w_fc = BufWriter::new(fs::File::create(args.out.join("fc_candidates.tsv"))?);
    writeln!(w_fc, "gene_5p\tgene_3p\tspanning\tsplit\tsource")?;
    for cand in &fc_candidates {
        writeln!(
            w_fc,
            "{}\t{}\t{}\t{}\t{}",
            cand.gene_5p,
            cand.gene_3p,
            cand.spanning_pairs
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string()),
            cand.split_reads
                .map(|v| v.to_string())
                .unwrap_or_else(|| "-".to_string()),
            cand.source.display()
        )?;
    }

    // 6. Write Markdown Summary
    let mut w_sum = BufWriter::new(fs::File::create(
        args.out.join("sample_validation_summary.md"),
    )?);
    writeln!(w_sum, "# Sample Validation Summary\n")?;
    writeln!(w_sum, "- **FusionCatcher Total**: {}", fc_map.len())?;
    writeln!(w_sum, "- **STAR Candidates Total**: {}", star_map.len())?;
    writeln!(w_sum, "- **Shared**: {}", shared_count)?;
    writeln!(w_sum, "- **Only FusionCatcher**: {}", only_fc_count)?;
    writeln!(w_sum, "- **Only STAR (rinfuse)**: {}", only_star_count)?;

    if !focus_genes.is_empty() {
        writeln!(w_sum, "\n## Focus Genes ({})", args.focus_gene.join(", "))?;
        if missing_from_star.is_empty() {
            writeln!(w_sum, "No focus candidates were missed by STAR.")?;
        } else {
            writeln!(w_sum, "### Focus candidates missing from STAR:")?;
            for p in &missing_from_star {
                writeln!(w_sum, "- {} -- {}", p.0, p.1)?;
            }
        }
    }

    // 7. Write Manifest
    let manifest_path = args.out.join("manifest.json");
    let mut w_man = BufWriter::new(fs::File::create(&manifest_path)?);
    let manifest = json!({
        "fc_out": args.fc_out,
        "star_candidates": args.star_candidates,
        "reads": args.reads,
        "extracted_reads_count": extracted_reads_count,
        "counts": {
            "fc_candidates": fc_map.len(),
            "star_candidates": star_map.len(),
            "shared": shared_count,
            "only_fc": only_fc_count,
            "only_star": only_star_count
        }
    });
    serde_json::to_writer_pretty(&mut w_man, &manifest)?;

    eprintln!(
        "Validation complete: shared={} only_fc={} only_star={} -> {}",
        shared_count,
        only_fc_count,
        only_star_count,
        args.out.display()
    );

    Ok(())
}

fn load_star_candidates(path: &Path) -> Result<Vec<FusionCandidateLite>> {
    let mut candidates = Vec::new();

    // Check if it's a directory
    let actual_path = if path.is_dir() {
        let jsonl = path.join("star_candidates.jsonl");
        if jsonl.exists() {
            jsonl
        } else {
            anyhow::bail!(
                "star_candidates.jsonl not found in directory {}",
                path.display()
            );
        }
    } else {
        path.to_path_buf()
    };

    let file = fs::File::open(&actual_path)?;
    let reader = BufReader::new(file);

    // Try parsing as JSONL
    if actual_path.extension().and_then(|s| s.to_str()) == Some("jsonl") {
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let cand: FusionCandidateLite = serde_json::from_str(&line)?;
            candidates.push(cand);
        }
    } else {
        // Assume TSV. Let's just require JSONL for simplicity since aggregate-star writes it.
        anyhow::bail!("Please provide the JSONL star candidates file (star_candidates.jsonl) or the directory containing it.");
    }

    Ok(candidates)
}
