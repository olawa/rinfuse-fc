use crate::args::AggregateStarArgs;
use anyhow::Result;
use rinfuse_core::FusionCandidateLite;
use rinfuse_fc::steps::aggregate::aggregate_star_junctions;
use std::fs;
use std::io::BufWriter;
use std::path::Path;

pub fn run(args: AggregateStarArgs) -> Result<()> {
    if !args.out.exists() {
        fs::create_dir_all(&args.out)?;
    }

    let candidates = aggregate_star_junctions(&args.junctions, &args.genes)?;

    // Handle focus genes
    let focus_genes: Vec<String> = args
        .focus_gene
        .iter()
        .flat_map(|g| g.split(','))
        .map(|g| g.trim().to_uppercase())
        .collect();

    if !focus_genes.is_empty() {
        let mut focus_candidates = Vec::new();
        for c in &candidates {
            if focus_genes.contains(&c.gene_a.to_uppercase())
                || focus_genes.contains(&c.gene_b.to_uppercase())
            {
                focus_candidates.push(c.clone());
            }
        }

        write_outputs(&args.out, "focus_star", &focus_candidates)?;

        eprintln!(
            "Aggregated {} focus candidates (filtered by {:?}) -> {}",
            focus_candidates.len(),
            focus_genes,
            args.out.display()
        );
    }

    write_outputs(&args.out, "star", &candidates)?;

    eprintln!(
        "Aggregated {} total candidates -> {}",
        candidates.len(),
        args.out.display()
    );

    Ok(())
}

pub fn write_outputs(
    out_dir: &Path,
    prefix: &str,
    candidates: &[FusionCandidateLite],
) -> Result<()> {
    // Write JSONL
    let jsonl_path = out_dir.join(format!("{}_candidates.jsonl", prefix));
    let mut w_jsonl = BufWriter::new(fs::File::create(&jsonl_path)?);
    for c in candidates {
        let line = serde_json::to_string(c)?;
        use std::io::Write;
        writeln!(w_jsonl, "{}", line)?;
    }

    // Write TSV
    let tsv_path = out_dir.join(format!("{}_candidates.tsv", prefix));
    let mut w_tsv = BufWriter::new(fs::File::create(&tsv_path)?);
    use std::io::Write;
    writeln!(
        w_tsv,
        "gene_a\tgene_b\tgene_id_a\tgene_id_b\tchrom_a\tchrom_b\tsupport_junction_count\tunique_read_count\tmax_overhang\tjunction_types\texample_reads"
    )?;
    for c in candidates {
        let jtypes = c
            .junction_types
            .iter()
            .map(|t| t.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let reads = c.example_reads.join(",");
        writeln!(
            w_tsv,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            c.gene_a,
            c.gene_b,
            c.gene_id_a,
            c.gene_id_b,
            c.chrom_a,
            c.chrom_b,
            c.support_junction_count,
            c.unique_read_count,
            c.max_overhang,
            jtypes,
            reads
        )?;
    }

    // Write Summary JSON
    let summary_path = out_dir.join(format!("{}_candidate_summary.json", prefix));
    let mut w_sum = BufWriter::new(fs::File::create(&summary_path)?);
    let summary = serde_json::json!({
        "total_candidates": candidates.len(),
    });
    serde_json::to_writer_pretty(&mut w_sum, &summary)?;

    Ok(())
}
