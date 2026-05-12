use crate::args::ValidateCohortArgs;
use anyhow::{Context, Result};
use std::collections::HashSet;
use std::fs;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Debug)]
struct SampleMetrics {
    sample_id: String,
    total_fc: u64,
    total_star: u64,
    shared: u64,
    only_fc: u64,
    only_star: u64,
    focus_missing_count: u64,
    focus_shared_count: u64,
}

#[derive(Debug)]
struct CandidateRow {
    sample_id: String,
    gene_a: String,
    gene_b: String,
    status: String,
    fc_spanning: String,
    star_unique_reads: String,
    fc_source: String,
}

pub fn run(args: ValidateCohortArgs) -> Result<()> {
    if !args.out.exists() {
        fs::create_dir_all(&args.out)?;
    }

    let focus_genes: HashSet<String> = args
        .focus_gene
        .iter()
        .flat_map(|s| s.split(','))
        .map(|s| s.trim().to_uppercase())
        .filter(|s| !s.is_empty())
        .collect();

    let mut sample_metrics_list = Vec::new();
    let mut all_missing = Vec::new();
    let mut all_recovered = Vec::new();
    let mut focus_missing = Vec::new();

    for dir in &args.validation_dir {
        let sample_id = dir
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();

        let manifest_path = dir.join("manifest.json");
        let fc_vs_star_path = dir.join("fc_vs_star.tsv");

        if !manifest_path.exists() || !fc_vs_star_path.exists() {
            eprintln!(
                "Warning: skipping {} because manifest.json or fc_vs_star.tsv is missing.",
                dir.display()
            );
            continue;
        }

        // Parse manifest
        let manifest_content = fs::read_to_string(&manifest_path)?;
        let manifest: serde_json::Value = serde_json::from_str(&manifest_content)?;

        let counts = manifest.get("counts").context("manifest missing counts")?;
        let total_fc = counts
            .get("fc_candidates")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        let total_star = counts
            .get("star_candidates")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);
        let shared = counts.get("shared").and_then(|v| v.as_u64()).unwrap_or(0);
        let only_fc = counts.get("only_fc").and_then(|v| v.as_u64()).unwrap_or(0);
        let only_star = counts
            .get("only_star")
            .and_then(|v| v.as_u64())
            .unwrap_or(0);

        let mut focus_missing_count = 0;
        let mut focus_shared_count = 0;

        // Parse fc_vs_star.tsv
        let file = fs::File::open(&fc_vs_star_path)?;
        let reader = BufReader::new(file);
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if i == 0 || line.trim().is_empty() {
                continue; // skip header or empty
            }

            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 6 {
                continue;
            }

            let gene_a = parts[0].to_string();
            let gene_b = parts[1].to_string();
            let status = parts[2].to_string();

            let row = CandidateRow {
                sample_id: sample_id.clone(),
                gene_a: gene_a.clone(),
                gene_b: gene_b.clone(),
                status: status.clone(),
                fc_spanning: parts[3].to_string(),
                star_unique_reads: parts[4].to_string(),
                fc_source: parts[5].to_string(),
            };

            let is_focus = focus_genes.contains(&gene_a.to_uppercase())
                || focus_genes.contains(&gene_b.to_uppercase());

            if status == "only_fc" {
                all_missing.push(row);
                if is_focus {
                    focus_missing_count += 1;
                }
            } else if status == "only_star" {
                all_recovered.push(row);
            } else if status == "shared" && is_focus {
                focus_shared_count += 1;
            }
        }

        // Re-read for focus_missing (we had to clone or parse row, just adding to a list is easier)
        // Actually, we can just push clones of row into the respective lists
        let reader = BufReader::new(fs::File::open(&fc_vs_star_path)?);
        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            if i == 0 || line.trim().is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 6 {
                continue;
            }

            let gene_a = parts[0].to_string();
            let gene_b = parts[1].to_string();
            let status = parts[2].to_string();
            let is_focus = focus_genes.contains(&gene_a.to_uppercase())
                || focus_genes.contains(&gene_b.to_uppercase());

            if is_focus && status == "only_fc" {
                focus_missing.push(CandidateRow {
                    sample_id: sample_id.clone(),
                    gene_a,
                    gene_b,
                    status,
                    fc_spanning: parts[3].to_string(),
                    star_unique_reads: parts[4].to_string(),
                    fc_source: parts[5].to_string(),
                });
            }
        }

        sample_metrics_list.push(SampleMetrics {
            sample_id,
            total_fc,
            total_star,
            shared,
            only_fc,
            only_star,
            focus_missing_count,
            focus_shared_count,
        });
    }

    // Sort metrics by sample ID for stable output
    sample_metrics_list.sort_by(|a, b| a.sample_id.cmp(&b.sample_id));

    // Write cohort_summary.tsv
    let mut w_summary = BufWriter::new(fs::File::create(args.out.join("cohort_summary.tsv"))?);
    writeln!(w_summary, "sample_id\ttotal_fc\ttotal_star\tshared\tonly_fc\tonly_star\tfocus_missing_count\tfocus_shared_count")?;
    for m in &sample_metrics_list {
        writeln!(
            w_summary,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            m.sample_id,
            m.total_fc,
            m.total_star,
            m.shared,
            m.only_fc,
            m.only_star,
            m.focus_missing_count,
            m.focus_shared_count
        )?;
    }

    // Write lists
    write_rows(&args.out.join("all_missing_from_star.tsv"), &all_missing)?;
    write_rows(&args.out.join("all_recovered_by_star.tsv"), &all_recovered)?;
    write_rows(&args.out.join("focus_missing.tsv"), &focus_missing)?;

    // Write Markdown
    let mut w_md = BufWriter::new(fs::File::create(args.out.join("cohort_summary.md"))?);
    writeln!(w_md, "# Cohort Validation Summary\n")?;
    writeln!(
        w_md,
        "Total Samples Processed: {}",
        sample_metrics_list.len()
    )?;

    let total_missing: u64 = sample_metrics_list
        .iter()
        .map(|m| m.focus_missing_count)
        .sum();
    let total_shared: u64 = sample_metrics_list
        .iter()
        .map(|m| m.focus_shared_count)
        .sum();

    if !focus_genes.is_empty() {
        writeln!(
            w_md,
            "\n## Focus Genes Performance ({})",
            args.focus_gene.join(", ")
        )?;
        writeln!(
            w_md,
            "- **Total Focus Candidates Shared**: {}",
            total_shared
        )?;
        writeln!(
            w_md,
            "- **Total Focus Candidates Missing (FC Only)**: {}",
            total_missing
        )?;
    }

    writeln!(w_md, "\n## Sample Metrics Table")?;
    writeln!(w_md, "| Sample ID | Shared | Only FC | Only STAR |")?;
    writeln!(w_md, "|---|---|---|---|")?;
    for m in &sample_metrics_list {
        writeln!(
            w_md,
            "| {} | {} | {} | {} |",
            m.sample_id, m.shared, m.only_fc, m.only_star
        )?;
    }

    eprintln!(
        "Cohort validation complete: {} samples -> {}",
        sample_metrics_list.len(),
        args.out.display()
    );

    Ok(())
}

fn write_rows(path: &Path, rows: &[CandidateRow]) -> Result<()> {
    let mut w = BufWriter::new(fs::File::create(path)?);
    writeln!(
        w,
        "sample_id\tgene_a\tgene_b\tstatus\tfc_spanning\tstar_unique_reads\tfc_source"
    )?;
    for r in rows {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.sample_id,
            r.gene_a,
            r.gene_b,
            r.status,
            r.fc_spanning,
            r.star_unique_reads,
            r.fc_source
        )?;
    }
    Ok(())
}
