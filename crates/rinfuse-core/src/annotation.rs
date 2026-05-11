use crate::evidence::Strand;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

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
    pub gene_a: String,
    pub gene_b: String,
    pub gene_id_a: String,
    pub gene_id_b: String,
    pub chrom_a: String,
    pub chrom_b: String,
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
