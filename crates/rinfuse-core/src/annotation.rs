use crate::evidence::Strand;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct OrientedFusionPair {
    pub gene_5p: String,
    pub gene_3p: String,
}

impl OrientedFusionPair {
    pub fn new(gene_5p: &str, gene_3p: &str) -> Self {
        Self {
            gene_5p: gene_5p.to_string(),
            gene_3p: gene_3p.to_string(),
        }
    }

    pub fn uppercased(gene_5p: &str, gene_3p: &str) -> Self {
        Self {
            gene_5p: gene_5p.to_uppercase(),
            gene_3p: gene_3p.to_uppercase(),
        }
    }

    pub fn unordered(&self) -> UnorderedGenePair {
        UnorderedGenePair::new(&self.gene_5p, &self.gene_3p)
    }

    pub fn contains_gene(&self, gene: &str) -> bool {
        self.gene_5p == gene || self.gene_3p == gene
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct UnorderedGenePair {
    pub gene_a: String,
    pub gene_b: String,
}

impl UnorderedGenePair {
    pub fn new(a: &str, b: &str) -> Self {
        let mut genes = [a.to_string(), b.to_string()];
        genes.sort();
        Self {
            gene_a: genes[0].clone(),
            gene_b: genes[1].clone(),
        }
    }
}

/// A simple genomic interval representing a gene.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneInterval {
    pub chrom: String,
    pub start_0based: u64,
    pub end_0based: u64,
    pub strand: Strand,
    pub gene_id: String,
    pub gene_symbol: String,
}

/// A hit against a gene annotation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneHit {
    pub gene_symbol: String,
    pub gene_id: String,
}

/// Simple index for gene intervals, separated by chromosome.
#[derive(Debug, Default, Clone)]
pub struct GeneAnnotationIndex {
    // Mapping from chromosome to a list of intervals.
    // For MVP, we'll keep them sorted by start position.
    pub chrom_intervals: HashMap<String, Vec<GeneInterval>>,
}

impl GeneAnnotationIndex {
    pub fn new(mut intervals: Vec<GeneInterval>) -> Self {
        let mut chrom_intervals: HashMap<String, Vec<GeneInterval>> = HashMap::new();
        for interval in intervals.drain(..) {
            chrom_intervals
                .entry(interval.chrom.clone())
                .or_default()
                .push(interval);
        }

        // Sort intervals within each chromosome by start position
        for intervals in chrom_intervals.values_mut() {
            intervals.sort_by_key(|i| i.start_0based);
        }

        Self { chrom_intervals }
    }

    /// Look up a gene by chromosome and position.
    /// Returns the first overlapping gene found.
    pub fn lookup(&self, chrom: &str, pos: u64) -> Option<GeneHit> {
        if let Some(intervals) = self.chrom_intervals.get(chrom) {
            // Linear search for MVP, can optimize with binary search or interval tree later
            for interval in intervals {
                if pos >= interval.start_0based && pos < interval.end_0based {
                    return Some(GeneHit {
                        gene_symbol: interval.gene_symbol.clone(),
                        gene_id: interval.gene_id.clone(),
                    });
                }
                // Since it's sorted by start, if interval start > pos, we can stop
                if interval.start_0based > pos {
                    break;
                }
            }
        }
        None
    }
}

/// A preliminary fusion candidate aggregated from evidence.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FusionCandidateLite {
    #[serde(alias = "gene_a")]
    pub gene_5p: String,
    #[serde(alias = "gene_b")]
    pub gene_3p: String,
    pub unordered_gene_a: String,
    pub unordered_gene_b: String,
    #[serde(alias = "gene_id_a")]
    pub gene_id_5p: String,
    #[serde(alias = "gene_id_b")]
    pub gene_id_3p: String,
    #[serde(alias = "chrom_a")]
    pub chrom_5p: String,
    #[serde(alias = "chrom_b")]
    pub chrom_3p: String,
    pub support_junction_count: u32,
    pub unique_read_count: u32,
    pub max_overhang: u32,
    pub junction_types: Vec<i32>,
    pub example_reads: Vec<String>,
    pub source: String,
}

impl FusionCandidateLite {
    /// Create an unknown hit fallback
    pub fn unknown_hit(chrom: &str) -> GeneHit {
        GeneHit {
            gene_symbol: format!("UNKNOWN_{}", chrom),
            gene_id: "UNKNOWN".to_string(),
        }
    }
}
