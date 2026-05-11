use anyhow::{Context, Result};
use rinfuse_core::normalize_read_name;
use std::fs;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

/// A parsed row from a FusionCatcher candidate report with source info and counts.
#[derive(Debug, Clone)]
pub struct FcCandidateWithSource {
    pub gene_5p: String,
    pub gene_3p: String,
    pub spanning_pairs: Option<u32>,
    pub split_reads: Option<u32>,
    /// All raw TSV fields from the line.
    pub raw_fields: Vec<String>,
    pub source: PathBuf,
}

/// Detected layout of a FusionCatcher output directory.
#[derive(Debug, Default)]
pub struct FcOutputDir {
    pub root: PathBuf,
    /// Detected candidate/final-report files.
    pub candidate_reports: Vec<PathBuf>,
    /// All detected supporting-read-ID files.
    pub supporting_read_files: Vec<PathBuf>,
    /// All files encountered during scan.
    pub scanned_files: Vec<PathBuf>,
    /// Warnings encountered during scanning.
    pub warnings: Vec<String>,
}

const CANDIDATE_REPORT_NAMES: &[&str] = &[
    "final-list_candidate-fusion-genes.txt",
    "final-list_candidate-fusion-genes.tsv",
    "final-list_candidate-fusion-genes.csv",
];

fn is_candidate_report_name(name: &str) -> bool {
    if CANDIDATE_REPORT_NAMES.contains(&name) {
        return true;
    }
    let lower = name.to_ascii_lowercase();
    lower.starts_with("candidate_fusion-genes") && lower.ends_with(".txt")
}

fn is_supporting_read_file_name(name: &str) -> bool {
    let n = name.to_ascii_lowercase();
    (n.contains("supporting") && n.contains("reads"))
        || (n.contains("candidate_fusion") && n.contains("reads"))
        || (n.contains("read_ids") && n.ends_with(".txt"))
        || (n.contains("supporting") && n.contains("paired"))
}

impl FcOutputDir {
    /// Scan `dir` for FusionCatcher output files.
    pub fn detect(dir: &Path, recursive: bool, max_depth: usize) -> Result<Self> {
        let mut out = Self {
            root: dir.to_path_buf(),
            ..Default::default()
        };
        out.scan_dir(dir, 0, recursive, max_depth)?;
        out.candidate_reports.sort();
        out.supporting_read_files.sort();
        out.scanned_files.sort();
        Ok(out)
    }

    fn scan_dir(
        &mut self,
        dir: &Path,
        depth: usize,
        recursive: bool,
        max_depth: usize,
    ) -> Result<()> {
        if depth > max_depth {
            return Ok(());
        }

        let entries = match fs::read_dir(dir) {
            Ok(e) => e,
            Err(e) => {
                self.warnings
                    .push(format!("cannot read directory {}: {}", dir.display(), e));
                return Ok(());
            }
        };

        for entry in entries {
            let entry = entry?;
            let path = entry.path();
            let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");

            if path.is_dir() {
                if recursive {
                    self.scan_dir(&path, depth + 1, recursive, max_depth)?;
                }
            } else if path.is_file() {
                self.scanned_files.push(path.clone());
                if is_candidate_report_name(name) {
                    self.candidate_reports.push(path.clone());
                } else if is_supporting_read_file_name(name) {
                    self.supporting_read_files.push(path.clone());
                }
            }
        }
        Ok(())
    }

    /// Parse all detected candidate reports into rows. Returns empty vec if none found.
    pub fn parse_all_candidates(&self) -> Result<Vec<FcCandidateWithSource>> {
        let mut all_candidates = Vec::new();
        for path in &self.candidate_reports {
            let candidates = self.parse_candidate_file(path)?;
            all_candidates.extend(candidates);
        }
        Ok(all_candidates)
    }

    /// Parse a single candidate report file.
    pub fn parse_candidate_file(&self, path: &Path) -> Result<Vec<FcCandidateWithSource>> {
        let file = fs::File::open(path)
            .with_context(|| format!("cannot open candidate report {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut candidates = Vec::new();
        let mut is_rinfuse_format = false;

        for (i, line) in reader.lines().enumerate() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let fields: Vec<String> = trimmed.split('\t').map(|s| s.to_string()).collect();

            // Detect format from header (usually line 0)
            if i == 0
                && fields.len() >= 4
                && fields[0] == "Gene_5p"
                && fields[1] == "Gene_3p"
                && fields[2] == "Source"
            {
                is_rinfuse_format = true;
                continue;
            }

            // Skip header rows (col 0 starts with "Gene" or "gene")
            if fields
                .first()
                .map(|f| f.starts_with("Gene") || f.starts_with("gene") || f.starts_with('#'))
                .unwrap_or(false)
            {
                continue;
            }
            if fields.len() < 2 {
                continue;
            }

            let gene_5p = fields[0].clone();
            let gene_3p = fields[1].clone();
            let mut spanning_pairs = None;
            let mut split_reads = None;
            let raw_fields = fields.clone();
            let mut source = path.to_path_buf();

            if is_rinfuse_format && fields.len() >= 4 {
                // Rinfuse format: Gene_5p, Gene_3p, Source, Raw_Line...
                // The original fields are in the 4th column onwards if they were joined,
                // but wait, my write_candidates_tsv joined them with tabs.
                // So fields[3..] are the original fields.
                source = PathBuf::from(&fields[2]);
                if fields.len() >= 6 {
                    spanning_pairs = fields[5].parse::<u32>().ok();
                }
                if fields.len() >= 7 {
                    split_reads = fields[6].parse::<u32>().ok();
                }
            } else {
                // Original FC format
                spanning_pairs = fields.get(2).and_then(|s| s.parse::<u32>().ok());
                split_reads = fields.get(3).and_then(|s| s.parse::<u32>().ok());
            }

            candidates.push(FcCandidateWithSource {
                gene_5p,
                gene_3p,
                spanning_pairs,
                split_reads,
                raw_fields,
                source,
            });
        }

        Ok(candidates)
    }

    /// Collect normalized read ID base-names from all supporting read files.
    pub fn collect_read_id_tokens(&self) -> Result<Vec<String>> {
        let mut tokens: Vec<String> = Vec::new();
        for path in &self.supporting_read_files {
            collect_tokens_from_file(path, &mut tokens)?;
        }
        Ok(tokens)
    }
}

/// Scan a single file for read-ID-like tokens.
fn collect_tokens_from_file(path: &Path, out: &mut Vec<String>) -> Result<()> {
    let file = fs::File::open(path)
        .with_context(|| format!("cannot open supporting read file {}", path.display()))?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        for field in trimmed.split('\t') {
            for token in field.split(',') {
                let token = token.trim();
                if looks_like_read_id(token) {
                    let normalized = normalize_read_name(token);
                    if !normalized.base.is_empty() {
                        out.push(normalized.base.clone());
                    }
                }
            }
        }
    }
    Ok(())
}

/// Heuristic: does this token look like a sequencing read ID?
fn looks_like_read_id(token: &str) -> bool {
    if token.len() < 2 {
        return false;
    }
    // Pure integer -> count/coordinate, not read ID
    if token.chars().all(|c| c.is_ascii_digit()) {
        return false;
    }
    // Chromosome names
    if token.starts_with("chr") || token.starts_with("CHR") {
        return false;
    }
    // Ensembl identifiers
    if token.starts_with("ENSG")
        || token.starts_with("ENST")
        || token.starts_with("ENSE")
        || token.starts_with("ENSR")
    {
        return false;
    }
    // Known header field name prefixes
    const HEADER_PREFIXES: &[&str] = &[
        "Gene_",
        "Fusion_",
        "Spanning_",
        "Longest_",
        "Exon_",
        "Predicted_",
        "Breakpoint",
    ];
    for prefix in HEADER_PREFIXES {
        if token.starts_with(prefix) {
            return false;
        }
    }
    // Genomic coordinate: starts with digit + contains colon -> skip
    if token
        .chars()
        .next()
        .map(|c| c.is_ascii_digit())
        .unwrap_or(false)
        && token.contains(':')
    {
        return false;
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn looks_like_read_id_rejects_integers() {
        assert!(!looks_like_read_id("10"));
        assert!(!looks_like_read_id("0"));
    }

    #[test]
    fn looks_like_read_id_rejects_chromosomes() {
        assert!(!looks_like_read_id("chr1"));
        assert!(!looks_like_read_id("CHR22"));
    }

    #[test]
    fn looks_like_read_id_rejects_ensembl() {
        assert!(!looks_like_read_id("ENSG00000001"));
        assert!(!looks_like_read_id("ENST00000002"));
    }

    #[test]
    fn looks_like_read_id_rejects_header_fields() {
        assert!(!looks_like_read_id("Gene_1_symbol"));
        assert!(!looks_like_read_id("Spanning_pairs"));
    }

    #[test]
    fn looks_like_read_id_accepts_read_names() {
        assert!(looks_like_read_id("READ_A/1"));
        assert!(looks_like_read_id("SRR1234567.1/1"));
        assert!(looks_like_read_id("READ_B"));
    }

    #[test]
    fn is_candidate_report_name_works() {
        assert!(is_candidate_report_name(
            "final-list_candidate-fusion-genes.txt"
        ));
        assert!(is_candidate_report_name("candidate_fusion-genes_all.txt"));
        assert!(!is_candidate_report_name("supporting_reads.txt"));
    }

    #[test]
    fn is_supporting_read_file_name_works() {
        assert!(is_supporting_read_file_name(
            "supporting_reads_in_pairs.txt"
        ));
        assert!(is_supporting_read_file_name(
            "candidate_fusion_reads_ids.txt"
        ));
        assert!(!is_supporting_read_file_name(
            "final-list_candidate-fusion-genes.txt"
        ));
    }
}
