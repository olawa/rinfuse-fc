use anyhow::Result;
use rinfuse_core::{GeneInterval, Strand};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Default)]
pub struct AnnotationParseReport {
    pub total_lines: usize,
    pub parsed_ok: usize,
    pub skipped: usize,
    pub warnings: Vec<String>,
}

pub fn parse_gene_intervals(path: &Path) -> Result<(Vec<GeneInterval>, AnnotationParseReport)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut intervals = Vec::new();
    let mut report = AnnotationParseReport::default();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line_number = i + 1;
        report.total_lines += 1;

        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            report.skipped += 1;
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 6 {
            // Check if it's the header row
            if i == 0 && fields.first().is_some_and(|s| s.starts_with("chrom")) {
                report.skipped += 1;
                continue;
            }
            report.warnings.push(format!(
                "line {}: expected >=6 fields, got {}",
                line_number,
                fields.len()
            ));
            continue;
        }

        // Header check: chrom  start_0based  end_0based  strand  gene_id  gene_symbol
        if i == 0 && fields[0] == "chrom" {
            report.skipped += 1;
            continue;
        }

        let chrom = fields[0].to_string();

        let start_0based = match fields[1].parse::<u64>() {
            Ok(v) => v,
            Err(_) => {
                report.warnings.push(format!(
                    "line {}: invalid start '{}'",
                    line_number, fields[1]
                ));
                continue;
            }
        };

        let end_0based = match fields[2].parse::<u64>() {
            Ok(v) => v,
            Err(_) => {
                report
                    .warnings
                    .push(format!("line {}: invalid end '{}'", line_number, fields[2]));
                continue;
            }
        };

        if start_0based >= end_0based {
            report.warnings.push(format!(
                "line {}: start {} >= end {}",
                line_number, start_0based, end_0based
            ));
            continue;
        }

        let strand = match fields[3] {
            "+" => Strand::Plus,
            "-" => Strand::Minus,
            "." => Strand::Unknown,
            _ => {
                report.warnings.push(format!(
                    "line {}: invalid strand '{}'",
                    line_number, fields[3]
                ));
                continue;
            }
        };

        let gene_id = fields[4].to_string();
        let gene_symbol = fields[5].to_string();

        intervals.push(GeneInterval {
            chrom,
            start_0based,
            end_0based,
            strand,
            gene_id,
            gene_symbol,
        });
        report.parsed_ok += 1;
    }

    Ok((intervals, report))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_parse_gene_intervals() {
        let temp = tempdir().unwrap();
        let path = temp.path().join("genes.tsv");
        let content = "chrom\tstart_0based\tend_0based\tstrand\tgene_id\tgene_symbol\n\
chr1\t1000\t2000\t+\tENSG000001\tGENE_A\n\
# a comment\n\
chr2\t5000\t6000\t-\tENSG000002\tGENE_B\n\
chr3\tbad\t2000\t+\tE3\tG3\n\
chr4\t2000\t1000\t+\tE4\tG4\n\
chr5\t100\t200\tx\tE5\tG5\n";

        fs::write(&path, content).unwrap();

        let (intervals, report) = parse_gene_intervals(&path).unwrap();

        assert_eq!(intervals.len(), 2);
        assert_eq!(intervals[0].gene_symbol, "GENE_A");
        assert_eq!(intervals[1].gene_symbol, "GENE_B");

        // 7 lines total: 1 header, 2 valid, 1 comment, 3 invalid
        assert_eq!(report.total_lines, 7);
        assert_eq!(report.parsed_ok, 2);
        assert_eq!(report.skipped, 2); // header + comment
        assert_eq!(report.warnings.len(), 3);
        assert!(report.warnings[0].contains("invalid start"));
        assert!(report.warnings[1].contains("start 2000 >= end 1000"));
        assert!(report.warnings[2].contains("invalid strand"));
    }
}
