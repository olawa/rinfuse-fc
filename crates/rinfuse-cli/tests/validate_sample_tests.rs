#[cfg(test)]
mod validate_sample_tests {
    use rinfuse_cli::args::ValidateSampleArgs;
    use rinfuse_cli::commands::validate_sample::run;
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
    fn test_validate_sample() {
        let temp = tempdir().unwrap();
        let fc_out = fixtures_dir().join("fusioncatcher_minimal");
        let star_candidates = temp.path().join("star_candidates.jsonl");

        // Create mock star candidates
        // Shared: BCR-ABL1
        // STAR-only: DUX4-IGH
        let c1 = r#"{"gene_a":"ABL1","gene_b":"BCR","gene_id_a":"E1","gene_id_b":"E2","chrom_a":"chr9","chrom_b":"chr22","support_junction_count":10,"unique_read_count":8,"max_overhang":40,"junction_types":[1],"example_reads":["R1"],"source":"STAR"}"#;
        let c2 = r#"{"gene_a":"DUX4","gene_b":"IGH","gene_id_a":"E3","gene_id_b":"E4","chrom_a":"chr4","chrom_b":"chr14","support_junction_count":5,"unique_read_count":5,"max_overhang":30,"junction_types":[2],"example_reads":["R2"],"source":"STAR"}"#;
        fs::write(&star_candidates, format!("{}\n{}\n", c1, c2)).unwrap();

        let out = temp.path().join("validation_out");

        let args = ValidateSampleArgs {
            fc_out,
            star_candidates,
            reads: vec![],
            out: out.clone(),
            focus_gene: vec!["DUX4".to_string(), "ETV6".to_string()], // DUX4 is STAR only, ETV6 is FC only
        };

        run(args).unwrap();

        assert!(out.join("sample_validation_summary.md").exists());
        assert!(out.join("fc_candidates.tsv").exists());
        assert!(out.join("fc_vs_star.tsv").exists());
        assert!(out.join("focus_fc_vs_star.tsv").exists());
        assert!(out.join("missing_from_star.tsv").exists());
        assert!(out.join("recovered_by_star.tsv").exists());
        assert!(out.join("manifest.json").exists());

        let summary = fs::read_to_string(out.join("sample_validation_summary.md")).unwrap();
        assert!(summary.contains("- **FusionCatcher Total**: 2")); // ETV6-RUNX1, BCR-ABL1
        assert!(summary.contains("- **STAR Candidates Total**: 2")); // BCR-ABL1, DUX4-IGH
        assert!(summary.contains("- **Shared**: 1")); // BCR-ABL1
        assert!(summary.contains("- **Only FusionCatcher**: 1")); // ETV6-RUNX1
        assert!(summary.contains("- **Only STAR (rinfuse)**: 1")); // DUX4-IGH

        // ETV6 is missing from STAR
        assert!(summary.contains("- ETV6 -- RUNX1"));

        let missing = fs::read_to_string(out.join("missing_from_star.tsv")).unwrap();
        assert!(missing.contains("ETV6\tRUNX1\tonly_fc"));

        let recovered = fs::read_to_string(out.join("recovered_by_star.tsv")).unwrap();
        assert!(recovered.contains("DUX4\tIGH\tonly_star"));
    }
}
