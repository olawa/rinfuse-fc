#[cfg(test)]
mod validate_cohort_tests {
    use rinfuse_cli::args::ValidateCohortArgs;
    use rinfuse_cli::commands::validate_cohort::run;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_validate_cohort() {
        let temp = tempdir().unwrap();

        // Create sample1 validation output
        let s1_dir = temp.path().join("sample1");
        fs::create_dir_all(&s1_dir).unwrap();
        let s1_manifest = r#"{"counts":{"fc_candidates":2,"star_candidates":2,"shared":1,"only_fc":1,"only_star":1}}"#;
        fs::write(s1_dir.join("manifest.json"), s1_manifest).unwrap();
        let s1_tsv = "gene_a\tgene_b\tstatus\tfc_spanning\tstar_unique_reads\tfc_source\n\
        ABL1\tBCR\tshared\t5\t8\tfile1\n\
        ETV6\tRUNX1\tonly_fc\t3\t-\tfile2\n\
        DUX4\tIGH\tonly_star\t-\t5\t-\n";
        fs::write(s1_dir.join("fc_vs_star.tsv"), s1_tsv).unwrap();

        // Create sample2 validation output
        let s2_dir = temp.path().join("sample2");
        fs::create_dir_all(&s2_dir).unwrap();
        let s2_manifest = r#"{"counts":{"fc_candidates":1,"star_candidates":0,"shared":0,"only_fc":1,"only_star":0}}"#;
        fs::write(s2_dir.join("manifest.json"), s2_manifest).unwrap();
        let s2_tsv = "gene_a\tgene_b\tstatus\tfc_spanning\tstar_unique_reads\tfc_source\n\
        DUX4\tIGH\tonly_fc\t10\t-\tfile3\n";
        fs::write(s2_dir.join("fc_vs_star.tsv"), s2_tsv).unwrap();

        let out = temp.path().join("cohort_out");

        let args = ValidateCohortArgs {
            validation_dir: vec![s1_dir, s2_dir],
            out: out.clone(),
            focus_gene: vec!["DUX4".to_string()],
        };

        run(args).unwrap();

        assert!(out.join("cohort_summary.tsv").exists());
        assert!(out.join("cohort_summary.md").exists());
        assert!(out.join("all_missing_from_star.tsv").exists());
        assert!(out.join("all_recovered_by_star.tsv").exists());
        assert!(out.join("focus_missing.tsv").exists());

        // Check cohort summary tsv
        let summary_tsv = fs::read_to_string(out.join("cohort_summary.tsv")).unwrap();
        assert!(summary_tsv.contains("sample1\t2\t2\t1\t1\t1\t0\t0")); // sample1: DUX4 is only_star, not missing, shared=0
        assert!(summary_tsv.contains("sample2\t1\t0\t0\t1\t0\t1\t0")); // sample2: DUX4 is missing (only_fc), focus_missing=1

        // Check all_missing_from_star
        let missing = fs::read_to_string(out.join("all_missing_from_star.tsv")).unwrap();
        assert!(missing.contains("sample1\tETV6\tRUNX1\tonly_fc"));
        assert!(missing.contains("sample2\tDUX4\tIGH\tonly_fc"));

        // Check focus_missing
        let focus = fs::read_to_string(out.join("focus_missing.tsv")).unwrap();
        assert!(!focus.contains("ETV6")); // Not a focus gene
        assert!(focus.contains("sample2\tDUX4\tIGH\tonly_fc"));

        // Check cohort summary md
        let md = fs::read_to_string(out.join("cohort_summary.md")).unwrap();
        assert!(md.contains("Total Samples Processed: 2"));
        assert!(md.contains("- **Total Focus Candidates Missing (FC Only)**: 1"));
    }
}
