use anyhow::{Context, Result};
use rinfuse_core::{FusionCandidateLite, GeneAnnotationIndex, StarChimericJunction};
use rinfuse_io::annotation::parse_gene_intervals;
use rinfuse_io::star::parse_chimeric_junctions;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub fn aggregate_star_junctions(
    junctions_path: &Path,
    genes_path: &Path,
) -> Result<Vec<FusionCandidateLite>> {
    let (intervals, _) =
        parse_gene_intervals(genes_path).context("Failed to parse gene intervals")?;
    let index = GeneAnnotationIndex::new(intervals);

    let junctions = load_junctions(junctions_path).context("Failed to load STAR junctions")?;

    let mut candidate_map: HashMap<(String, String), CandidateAggregator> = HashMap::new();

    for j in junctions {
        let hit1 = index
            .lookup(&j.seg1.chrom, j.seg1.genomic_pos)
            .unwrap_or_else(|| FusionCandidateLite::unknown_hit(&j.seg1.chrom));

        let hit2 = index
            .lookup(&j.seg2.chrom, j.seg2.genomic_pos)
            .unwrap_or_else(|| FusionCandidateLite::unknown_hit(&j.seg2.chrom));

        // Create an order-insensitive key based on gene symbols
        let mut genes = [
            (
                hit1.gene_symbol.clone(),
                hit1.gene_id.clone(),
                j.seg1.chrom.clone(),
            ),
            (
                hit2.gene_symbol.clone(),
                hit2.gene_id.clone(),
                j.seg2.chrom.clone(),
            ),
        ];
        genes.sort_by(|a, b| a.0.cmp(&b.0)); // Sort by symbol

        let key = (genes[0].0.clone(), genes[1].0.clone());

        let agg = candidate_map
            .entry(key)
            .or_insert_with(|| CandidateAggregator {
                gene_a: genes[0].0.clone(),
                gene_b: genes[1].0.clone(),
                gene_id_a: genes[0].1.clone(),
                gene_id_b: genes[1].1.clone(),
                chrom_a: genes[0].2.clone(),
                chrom_b: genes[1].2.clone(),
                junctions: Vec::new(),
            });

        agg.junctions.push(j);
    }

    let mut candidates: Vec<FusionCandidateLite> = candidate_map
        .into_values()
        .map(|agg| agg.into_candidate())
        .collect();

    // Sort by unique read count descending for stable output
    candidates.sort_by(|a, b| b.unique_read_count.cmp(&a.unique_read_count));

    Ok(candidates)
}

fn load_junctions(path: &Path) -> Result<Vec<StarChimericJunction>> {
    let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");
    if ext == "jsonl" {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let mut junctions = Vec::new();
        for line in reader.lines() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let j: StarChimericJunction = serde_json::from_str(&line)?;
            junctions.push(j);
        }
        Ok(junctions)
    } else {
        // Assume raw Chimeric.out.junction
        let (junctions, _) = parse_chimeric_junctions(path)?;
        Ok(junctions)
    }
}

struct CandidateAggregator {
    gene_a: String,
    gene_b: String,
    gene_id_a: String,
    gene_id_b: String,
    chrom_a: String,
    chrom_b: String,
    junctions: Vec<StarChimericJunction>,
}

impl CandidateAggregator {
    fn into_candidate(self) -> FusionCandidateLite {
        let mut unique_reads = HashSet::new();
        let mut max_overhang = 0;
        let mut junction_types = HashSet::new();
        let mut example_reads = Vec::new();

        for j in &self.junctions {
            unique_reads.insert(j.read_name.clone());
            max_overhang = max_overhang.max(j.max_overhang);
            junction_types.insert(j.junction_type);

            if example_reads.len() < 5 {
                example_reads.push(j.read_name.clone());
            }
        }

        let mut junction_types: Vec<i32> = junction_types.into_iter().collect();
        junction_types.sort_unstable();

        FusionCandidateLite {
            gene_a: self.gene_a,
            gene_b: self.gene_b,
            gene_id_a: self.gene_id_a,
            gene_id_b: self.gene_id_b,
            chrom_a: self.chrom_a,
            chrom_b: self.chrom_b,
            support_junction_count: self.junctions.len() as u32,
            unique_read_count: unique_reads.len() as u32,
            max_overhang,
            junction_types,
            example_reads,
            source: "STAR".to_string(),
        }
    }
}
