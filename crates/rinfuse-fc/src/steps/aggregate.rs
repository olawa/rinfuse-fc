use anyhow::{Context, Result};
use rinfuse_core::{
    FusionCandidateLite, GeneAnnotationIndex, OrientedFusionPair, StarChimericJunction,
};
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

    let mut candidate_map: HashMap<OrientedFusionPair, CandidateAggregator> = HashMap::new();

    for j in junctions {
        let seg1_lookup = j.pos1_lookup_0based().unwrap_or(0);
        let seg2_lookup = j.pos2_lookup_0based().unwrap_or(0);
        let hit1 = index
            .lookup(&j.seg1.chrom, seg1_lookup)
            .unwrap_or_else(|| FusionCandidateLite::unknown_hit(&j.seg1.chrom));

        let hit2 = index
            .lookup(&j.seg2.chrom, seg2_lookup)
            .unwrap_or_else(|| FusionCandidateLite::unknown_hit(&j.seg2.chrom));

        let key = OrientedFusionPair::new(&hit1.gene_symbol, &hit2.gene_symbol);

        let agg = candidate_map
            .entry(key.clone())
            .or_insert_with(|| CandidateAggregator {
                pair: key,
                gene_id_5p: hit1.gene_id.clone(),
                gene_id_3p: hit2.gene_id.clone(),
                chrom_5p: j.seg1.chrom.clone(),
                chrom_3p: j.seg2.chrom.clone(),
                junctions: Vec::new(),
            });

        agg.junctions.push(j);
    }

    let mut candidates: Vec<FusionCandidateLite> = candidate_map
        .into_values()
        .map(|agg| agg.into_candidate())
        .collect();

    // Sort by unique read count descending for stable output
    candidates.sort_by(|a, b| {
        b.unique_read_count
            .cmp(&a.unique_read_count)
            .then(a.gene_5p.cmp(&b.gene_5p))
            .then(a.gene_3p.cmp(&b.gene_3p))
    });

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
    pair: OrientedFusionPair,
    gene_id_5p: String,
    gene_id_3p: String,
    chrom_5p: String,
    chrom_3p: String,
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
            max_overhang = max_overhang.max(j.max_overhang.unwrap_or(0));
            junction_types.insert(j.junction_type);

            if example_reads.len() < 5 {
                example_reads.push(j.read_name.clone());
            }
        }

        let mut junction_types: Vec<i32> = junction_types.into_iter().collect();
        junction_types.sort_unstable();

        let unordered = self.pair.unordered();
        let gene_5p = self.pair.gene_5p;
        let gene_3p = self.pair.gene_3p;
        let gene_id_5p = self.gene_id_5p;
        let gene_id_3p = self.gene_id_3p;
        let chrom_5p = self.chrom_5p;
        let chrom_3p = self.chrom_3p;

        FusionCandidateLite {
            gene_5p,
            gene_3p,
            unordered_gene_a: unordered.gene_a,
            unordered_gene_b: unordered.gene_b,
            gene_id_5p,
            gene_id_3p,
            chrom_5p,
            chrom_3p,
            support_junction_count: self.junctions.len() as u32,
            unique_read_count: unique_reads.len() as u32,
            max_overhang,
            junction_types,
            example_reads,
            source: "STAR".to_string(),
        }
    }
}
