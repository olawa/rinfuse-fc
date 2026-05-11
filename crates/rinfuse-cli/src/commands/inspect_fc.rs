use crate::args::InspectFcArgs;
use anyhow::{Context, Result};
use rinfuse_io::fc_intermediates::{FcCandidateWithSource, FcOutputDir};
use rinfuse_io::{ReadRegistry, ReadRegistryBuildOptions};
use serde_json::json;
use std::collections::HashSet;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn run(args: InspectFcArgs) -> Result<()> {
    if !args.out.exists() {
        fs::create_dir_all(&args.out)
            .with_context(|| format!("failed to create output directory {}", args.out.display()))?;
    }

    let fc_dir = FcOutputDir::detect(&args.fc_out, args.recursive, args.max_depth)?;
    let all_candidates = fc_dir.parse_all_candidates()?;
    let read_tokens = fc_dir.collect_read_id_tokens()?;

    // 1. Focus gene filtering
    let focus_genes: HashSet<String> = args
        .focus_gene
        .iter()
        .flat_map(|s| s.split(',').map(|g| g.trim().to_string()))
        .filter(|s| !s.is_empty())
        .collect();

    let (focused_candidates, focus_read_tokens) = if !focus_genes.is_empty() {
        let focused: Vec<FcCandidateWithSource> = all_candidates
            .iter()
            .filter(|c| focus_genes.contains(&c.gene_5p) || focus_genes.contains(&c.gene_3p))
            .cloned()
            .collect();

        // For focus genes, we'll extract reads that are mentioned in files but for now
        // we'll just keep all read_tokens if no more granular filtering is available.
        // In a real scenario, we might want to only keep tokens associated with focused candidates.
        (focused, read_tokens.clone())
    } else {
        (Vec::new(), Vec::new())
    };

    // 2. Write candidates.tsv
    write_candidates_tsv(&args.out.join("candidates.tsv"), &all_candidates)?;

    if !focused_candidates.is_empty() {
        write_candidates_tsv(&args.out.join("focus_candidates.tsv"), &focused_candidates)?;
    }

    // 3. Extract reads via ReadRegistry
    let mut main_registry = ReadRegistry::from_tokens(&read_tokens);
    let opts = ReadRegistryBuildOptions::default();
    if !args.reads.is_empty() {
        main_registry.collect_from_fastq_paths(&args.reads, &opts)?;
    }
    main_registry.write_fastq(&args.out.join("supporting_reads.fq"))?;
    main_registry.write_missing(&args.out.join("missing_read_ids.txt"))?;

    let mut focus_found_count = 0;
    if !focus_read_tokens.is_empty() {
        let mut focus_registry = ReadRegistry::from_tokens(&focus_read_tokens);
        if !args.reads.is_empty() {
            focus_registry.collect_from_fastq_paths(&args.reads, &opts)?;
        }
        focus_registry.write_fastq(&args.out.join("focus_supporting_reads.fq"))?;
        focus_found_count = focus_registry.found_record_count();

        let mut focus_id_writer =
            BufWriter::new(fs::File::create(args.out.join("focus_read_ids.txt"))?);
        for token in &focus_read_tokens {
            writeln!(focus_id_writer, "{}", token)?;
        }
    }

    // 4. Write evidence.jsonl (one line per candidate)
    let evidence_path = args.out.join("evidence.jsonl");
    let mut ev_writer = BufWriter::new(fs::File::create(&evidence_path)?);
    for cand in &all_candidates {
        let line = json!({
            "gene_5p": cand.gene_5p,
            "gene_3p": cand.gene_3p,
            "source_file": cand.source,
            "evidence_type": "fusioncatcher_candidate",
        });
        writeln!(ev_writer, "{}", line)?;
    }
    ev_writer.flush()?;

    // 5. Write manifest.json
    let manifest = json!({
        "root": args.fc_out,
        "input_reads": args.reads,
        "scanned_files": fc_dir.scanned_files,
        "candidate_reports": fc_dir.candidate_reports,
        "supporting_read_files": fc_dir.supporting_read_files,
        "warnings": fc_dir.warnings,
        "counts": {
            "total_candidates": all_candidates.len(),
            "focused_candidates": focused_candidates.len(),
            "read_id_tokens": read_tokens.len(),
            "unique_bases": main_registry.requested_base_count(),
            "found_records": main_registry.found_record_count(),
            "missing_bases": main_registry.missing_base_count(),
            "focus_found_records": focus_found_count,
        }
    });

    let mut manifest_writer = BufWriter::new(fs::File::create(args.out.join("manifest.json"))?);
    serde_json::to_writer_pretty(&mut manifest_writer, &manifest)?;
    manifest_writer.flush()?;

    // 6. Write inspect_summary.md
    write_summary_md(
        &args.out.join("inspect_summary.md"),
        &fc_dir,
        &all_candidates,
        &focused_candidates,
        &main_registry,
    )?;

    eprintln!(
        "inspect-fc complete: candidates={} tokens={} found_reads={} missing={}",
        all_candidates.len(),
        read_tokens.len(),
        main_registry.found_record_count(),
        main_registry.missing_base_count()
    );

    Ok(())
}

fn write_candidates_tsv(path: &Path, candidates: &[FcCandidateWithSource]) -> Result<()> {
    let mut writer = BufWriter::new(fs::File::create(path)?);
    writeln!(writer, "Gene_5p\tGene_3p\tSource\tRaw_Line")?;
    for cand in candidates {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}",
            cand.gene_5p,
            cand.gene_3p,
            cand.source.display(),
            cand.raw_fields.join("\t")
        )?;
    }
    writer.flush()?;
    Ok(())
}

fn write_summary_md(
    path: &Path,
    fc_dir: &FcOutputDir,
    all_candidates: &[FcCandidateWithSource],
    focused_candidates: &[FcCandidateWithSource],
    registry: &ReadRegistry,
) -> Result<()> {
    let mut writer = BufWriter::new(fs::File::create(path)?);
    writeln!(writer, "# rinfuse-fc Inspection Summary\n")?;
    writeln!(writer, "## Discovery")?;
    writeln!(writer, "- **Root**: `{}`", fc_dir.root.display())?;
    writeln!(
        writer,
        "- **Scanned Files**: {}",
        fc_dir.scanned_files.len()
    )?;
    writeln!(
        writer,
        "- **Candidate Reports**: {}",
        fc_dir.candidate_reports.len()
    )?;
    for report in &fc_dir.candidate_reports {
        writeln!(writer, "  - `{}`", report.display())?;
    }
    writeln!(
        writer,
        "- **Supporting Read Files**: {}",
        fc_dir.supporting_read_files.len()
    )?;
    for file in &fc_dir.supporting_read_files {
        writeln!(writer, "  - `{}`", file.display())?;
    }

    writeln!(writer, "\n## Results")?;
    writeln!(writer, "- **Total Candidates**: {}", all_candidates.len())?;
    if !focused_candidates.is_empty() {
        writeln!(
            writer,
            "- **Focused Candidates**: {}",
            focused_candidates.len()
        )?;
    }
    writeln!(
        writer,
        "- **Unique Read Bases**: {}",
        registry.requested_base_count()
    )?;
    writeln!(
        writer,
        "- **Extracted FASTQ Records**: {}",
        registry.found_record_count()
    )?;
    writeln!(
        writer,
        "- **Missing IDs**: {}",
        registry.missing_base_count()
    )?;

    if !focused_candidates.is_empty() {
        writeln!(writer, "\n## Focus Genes")?;
        writeln!(writer, "| Gene 5' | Gene 3' | Source |")?;
        writeln!(writer, "|---|---|---|")?;
        for cand in focused_candidates {
            writeln!(
                writer,
                "| {} | {} | {} |",
                cand.gene_5p,
                cand.gene_3p,
                cand.source.display()
            )?;
        }
    }

    if !fc_dir.warnings.is_empty() {
        writeln!(writer, "\n## Warnings")?;
        for warning in &fc_dir.warnings {
            writeln!(writer, "- {}", warning)?;
        }
    }

    writer.flush()?;
    Ok(())
}
