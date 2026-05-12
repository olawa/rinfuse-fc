use anyhow::Result;
use rinfuse_core::{ChimericSegment, EvidenceSource, StarChimericJunction, Strand};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Summary of a Chimeric.out.junction parse run.
#[derive(Debug, Default, Serialize, Deserialize)]
pub struct StarParseReport {
    pub total_lines: usize,
    pub parsed_ok: usize,
    pub skipped_empty: usize,
    pub parse_warnings: Vec<ParseWarning>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ParseWarning {
    pub line_number: usize,
    pub reason: String,
    pub raw_line: String,
}

/// Parse a STAR `Chimeric.out.junction` file.
///
/// Returns all successfully parsed junctions plus a report with warnings.
pub fn parse_chimeric_junctions(
    path: &Path,
) -> Result<(Vec<StarChimericJunction>, StarParseReport)> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);

    let mut junctions = Vec::new();
    let mut report = StarParseReport::default();

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line_number = i + 1;
        report.total_lines += 1;

        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            report.skipped_empty += 1;
            continue;
        }

        let fields: Vec<&str> = trimmed.split('\t').collect();

        match parse_junction_row(line_number, &fields, trimmed) {
            Ok(junction) => {
                report.parsed_ok += 1;
                junctions.push(junction);
            }
            Err(reason) => {
                report.parse_warnings.push(ParseWarning {
                    line_number,
                    reason,
                    raw_line: trimmed.to_string(),
                });
            }
        }
    }

    Ok((junctions, report))
}

/// STAR Chimeric.out.junction columns (0-indexed):
///  0  chrom1
///  1  pos1
///  2  strand1
///  3  chrom2
///  4  pos2
///  5  strand2
///  6  junction_type  (-1 not canonical, 0 not determined, 1 GT/AG, 2 CT/AC)
///  7  repeat_left
///  8  repeat_right
///  9  read_name
/// 10  start1  (position in read of segment 1)
/// 11  cigar1
/// 12  start2  (position in read of segment 2)
/// 13  cigar2
/// 14  num_chimeric_reads
/// 15  max_overhang
fn parse_junction_row(
    line_number: usize,
    fields: &[&str],
    _raw: &str,
) -> std::result::Result<StarChimericJunction, String> {
    if fields.len() < 14 {
        return Err(format!(
            "line {}: expected >=14 fields, got {}",
            line_number,
            fields.len()
        ));
    }

    let chrom1 = fields[0].to_string();
    let pos1 = fields[1]
        .parse::<u64>()
        .map_err(|_| format!("line {}: invalid pos1 '{}'", line_number, fields[1]))?;
    let strand1 = parse_strand(fields[2], line_number, "strand1")?;

    let chrom2 = fields[3].to_string();
    let pos2 = fields[4]
        .parse::<u64>()
        .map_err(|_| format!("line {}: invalid pos2 '{}'", line_number, fields[4]))?;
    let strand2 = parse_strand(fields[5], line_number, "strand2")?;

    let junction_type = fields[6].parse::<i32>().map_err(|_| {
        format!(
            "line {}: invalid junction_type '{}'",
            line_number, fields[6]
        )
    })?;

    let repeat_left = fields[7]
        .parse::<u32>()
        .map_err(|_| format!("line {}: invalid repeat_left '{}'", line_number, fields[7]))?;
    let repeat_right = fields[8]
        .parse::<u32>()
        .map_err(|_| format!("line {}: invalid repeat_right '{}'", line_number, fields[8]))?;

    let read_name = fields[9].to_string();

    let start1 = fields[10]
        .parse::<u64>()
        .map_err(|_| format!("line {}: invalid start1 '{}'", line_number, fields[10]))?;
    let cigar1 = fields[11].to_string();

    let start2 = fields[12]
        .parse::<u64>()
        .map_err(|_| format!("line {}: invalid start2 '{}'", line_number, fields[12]))?;
    let cigar2 = fields[13].to_string();

    let num_chimeric_reads = fields
        .get(14)
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0);
    let max_overhang = fields
        .get(15)
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0);

    let raw_fields = fields.iter().map(|s| s.to_string()).collect();

    Ok(StarChimericJunction {
        seg1: ChimericSegment {
            chrom: chrom1,
            genomic_pos: pos1,
            strand: strand1,
            read_start: start1,
            cigar: cigar1,
        },
        seg2: ChimericSegment {
            chrom: chrom2,
            genomic_pos: pos2,
            strand: strand2,
            read_start: start2,
            cigar: cigar2,
        },
        junction_type,
        repeat_left,
        repeat_right,
        read_name,
        num_chimeric_reads,
        max_overhang,
        source: EvidenceSource::Star,
        raw_fields,
    })
}

fn parse_strand(s: &str, line: usize, field: &str) -> std::result::Result<Strand, String> {
    match s {
        "+" => Ok(Strand::Plus),
        "-" => Ok(Strand::Minus),
        "." => Ok(Strand::Unknown),
        other => Err(format!("line {}: invalid {}: '{}'", line, field, other)),
    }
}

/// Write junctions as newline-delimited JSON (one JSON object per line).
pub fn write_junctions_jsonl(path: &Path, junctions: &[StarChimericJunction]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    for j in junctions {
        let line = serde_json::to_string(j)?;
        writeln!(w, "{}", line)?;
    }
    Ok(())
}

/// Write junctions as a TSV file with a header row.
pub fn write_junctions_tsv(path: &Path, junctions: &[StarChimericJunction]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "{}", StarChimericJunction::tsv_header())?;
    for j in junctions {
        writeln!(w, "{}", j.tsv_row())?;
    }
    Ok(())
}

/// Write a parse summary as JSON.
pub fn write_parse_summary(path: &Path, report: &StarParseReport) -> Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    serde_json::to_writer_pretty(&mut w, report)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURE_LINE: &str =
        "chr2\t33459168\t+\tchr17\t37880992\t-\t1\t0\t0\tREAD_1/1\t20\t20M30S\t50\t30M20S\t5\t40";

    #[test]
    fn parses_valid_row() {
        let fields: Vec<&str> = FIXTURE_LINE.split('\t').collect();
        let result = parse_junction_row(1, &fields, FIXTURE_LINE);
        assert!(result.is_ok(), "{:?}", result);
        let j = result.unwrap();
        assert_eq!(j.seg1.chrom, "chr2");
        assert_eq!(j.seg1.genomic_pos, 33_459_168);
        assert_eq!(j.seg1.strand, Strand::Plus);
        assert_eq!(j.seg2.chrom, "chr17");
        assert_eq!(j.seg2.strand, Strand::Minus);
        assert_eq!(j.junction_type, 1);
        assert_eq!(j.num_chimeric_reads, 5);
        assert_eq!(j.max_overhang, 40);
        assert_eq!(j.read_name, "READ_1/1");
    }

    #[test]
    fn short_row_produces_error_not_panic() {
        let fields = vec!["chr2", "123"];
        let result = parse_junction_row(2, &fields, "chr2\t123");
        assert!(result.is_err());
        let msg = result.unwrap_err();
        assert!(msg.contains("expected >=14 fields"));
    }

    #[test]
    fn invalid_pos_produces_error() {
        let mut parts: Vec<&str> = FIXTURE_LINE.split('\t').collect();
        parts[1] = "NOT_A_NUMBER";
        let result = parse_junction_row(3, &parts, "");
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("invalid pos1"));
    }

    #[test]
    fn missing_optional_cols_defaults_to_zero() {
        // Only 14 fields — no num_chimeric_reads or max_overhang
        let short: Vec<&str> = FIXTURE_LINE.split('\t').take(14).collect();
        let result = parse_junction_row(4, &short, "");
        assert!(result.is_ok());
        let j = result.unwrap();
        assert_eq!(j.num_chimeric_reads, 0);
        assert_eq!(j.max_overhang, 0);
    }

    #[test]
    fn tsv_row_roundtrip_header_columns() {
        let fields: Vec<&str> = FIXTURE_LINE.split('\t').collect();
        let j = parse_junction_row(1, &fields, FIXTURE_LINE).unwrap();
        let row = j.tsv_row();
        let header = StarChimericJunction::tsv_header();
        let hcols: Vec<&str> = header.split('\t').collect();
        let rcols: Vec<&str> = row.split('\t').collect();
        assert_eq!(hcols.len(), rcols.len(), "header/row column count mismatch");
    }
}
