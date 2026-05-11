#[cfg(test)]
mod compare_tests {
    use std::fs;
    use std::path::PathBuf;
    use tempfile::tempdir;

    fn fixtures_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .unwrap()
            .parent()
            .unwrap()
            .join("tests/fixtures")
    }

    #[test]
    fn compare_identical() {
        let fixtures = fixtures_dir().join("fusioncatcher_minimal");
        let temp = tempdir().unwrap();
        let out_tsv = temp.path().join("compare.tsv");

        // Compare fixture against itself
        let args = rinfuse_cli::args::CompareArgs {
            fc: fixtures.clone(),
            rs: fixtures.clone(),
            out: out_tsv.clone(),
            focus_gene: vec![],
        };

        rinfuse_cli::commands::compare::run(args).unwrap();

        assert!(out_tsv.exists());
        let content = fs::read_to_string(&out_tsv).unwrap();
        // Should have 2 both, 0 only
        assert!(content.contains("both"));
        assert!(!content.contains("only_fc"));
        assert!(!content.contains("only_rs"));

        let md_path = out_tsv.with_extension("md");
        assert!(md_path.exists());
        let md_content = fs::read_to_string(md_path).unwrap();
        assert!(md_content.contains("- **Shared**: 2"));
    }

    #[test]
    fn compare_diff() {
        let fixtures = fixtures_dir().join("fusioncatcher_minimal");
        let temp = tempdir().unwrap();

        // Create a custom TSV for RS
        let rs_tsv = temp.path().join("rs_candidates.tsv");
        fs::write(&rs_tsv, "Gene_5p\tGene_3p\tSource\tRaw\nONLY_RS\tGENE\tsource\tONLY_RS\tGENE\t1\t1\nBCR\tABL1\tsource\tBCR\tABL1\t10\t5\n").unwrap();

        let out_tsv = temp.path().join("compare.tsv");
        let args = rinfuse_cli::args::CompareArgs {
            fc: fixtures,
            rs: rs_tsv,
            out: out_tsv.clone(),
            focus_gene: vec!["ONLY_RS".to_string()],
        };

        rinfuse_cli::commands::compare::run(args).unwrap();

        let content = fs::read_to_string(&out_tsv).unwrap();
        assert!(content.contains("both")); // BCR-ABL1
        assert!(content.contains("only_fc")); // ETV6-RUNX1 (from fixture)
        assert!(content.contains("only_rs")); // ONLY_RS-GENE

        let md_content = fs::read_to_string(out_tsv.with_extension("md")).unwrap();
        assert!(md_content.contains("## Focus Genes (ONLY_RS)"));
        assert!(md_content.contains("| GENE | ONLY_RS | only_rs |"));
    }
}
