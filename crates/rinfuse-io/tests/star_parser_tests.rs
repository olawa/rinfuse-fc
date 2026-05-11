use rinfuse_io::star::{
    parse_chimeric_junctions, write_junctions_jsonl, write_junctions_tsv, write_parse_summary,
};
use std::fs;
use std::path::PathBuf;
use tempfile::tempdir;

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/Chimeric.out.junction")
}

#[test]
fn parses_fixture_file_three_valid_rows() {
    let path = fixture_path();
    let (junctions, report) = parse_chimeric_junctions(&path).unwrap();

    // fixture has 3 data rows + 1 comment + 1 blank -> 3 parsed
    assert_eq!(report.parsed_ok, 3, "expected 3 parsed junctions");
    assert!(
        report.parse_warnings.is_empty(),
        "unexpected warnings: {:?}",
        report.parse_warnings
    );
    assert_eq!(report.skipped_empty, 1); // 1 comment line

    let first = &junctions[0];
    assert_eq!(first.seg1.chrom, "chr2");
    assert_eq!(first.seg2.chrom, "chr17");
    assert_eq!(first.num_chimeric_reads, 5);
    assert_eq!(first.max_overhang, 40);
}

#[test]
fn malformed_row_produces_warning_not_panic() {
    let temp = tempdir().unwrap();
    let bad_file = temp.path().join("bad.junction");
    fs::write(&bad_file, "chr1\tBAD_POS\t+\tchr2\t100\t-\t1\t0\t0\tREAD\t0\t50M\t0\t50M\n")
        .unwrap();

    let (junctions, report) = parse_chimeric_junctions(&bad_file).unwrap();

    assert_eq!(junctions.len(), 0);
    assert_eq!(report.parse_warnings.len(), 1);
    assert!(report.parse_warnings[0].reason.contains("invalid pos1"));
}

#[test]
fn write_jsonl_produces_valid_lines() {
    let temp = tempdir().unwrap();
    let path = fixture_path();
    let (junctions, _) = parse_chimeric_junctions(&path).unwrap();

    let out = temp.path().join("out.jsonl");
    write_junctions_jsonl(&out, &junctions).unwrap();

    let content = fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert_eq!(lines.len(), 3);
    // Each line should be valid JSON
    for line in &lines {
        let parsed: serde_json::Value = serde_json::from_str(line)
            .unwrap_or_else(|e| panic!("invalid JSON: {e}\nline: {line}"));
        assert!(parsed.get("seg1").is_some());
    }
}

#[test]
fn write_tsv_has_header_and_rows() {
    let temp = tempdir().unwrap();
    let path = fixture_path();
    let (junctions, _) = parse_chimeric_junctions(&path).unwrap();

    let out = temp.path().join("out.tsv");
    write_junctions_tsv(&out, &junctions).unwrap();

    let content = fs::read_to_string(&out).unwrap();
    let lines: Vec<&str> = content.lines().collect();
    assert!(
        lines[0].starts_with("chrom1\t"),
        "first line should be header"
    );
    assert_eq!(lines.len(), 4); // header + 3 data rows
}

#[test]
fn write_parse_summary_is_valid_json() {
    let temp = tempdir().unwrap();
    let path = fixture_path();
    let (_, report) = parse_chimeric_junctions(&path).unwrap();

    let out = temp.path().join("summary.json");
    write_parse_summary(&out, &report).unwrap();

    let content = fs::read_to_string(&out).unwrap();
    let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();
    assert_eq!(parsed["parsed_ok"], 3);
    assert_eq!(parsed["skipped_empty"], 1);
}
