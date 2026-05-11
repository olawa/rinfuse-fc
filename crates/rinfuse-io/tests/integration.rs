#[cfg(test)]
mod extract_reads_integration {
    use rinfuse_io::{ReadRegistry, ReadRegistryBuildOptions};
    use std::path::PathBuf;

    fn fixtures_dir() -> PathBuf {
        // Relative to workspace root
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent() // crates/
            .unwrap()
            .parent() // workspace root
            .unwrap()
            .join("tests/fixtures")
    }

    #[test]
    fn extract_reads_finds_both_mates() {
        let fixtures = fixtures_dir();
        let r1 = fixtures.join("reads_R1.fq");
        let r2 = fixtures.join("reads_R2.fq");
        let ids = fixtures.join("read_ids.txt");

        let mut registry = ReadRegistry::from_read_id_file(&ids).unwrap();
        let opts = ReadRegistryBuildOptions::default();
        registry.collect_from_fastq_paths(&[r1, r2], &opts).unwrap();

        // READ_A and READ_C match; READ_MISSING does not.
        assert_eq!(registry.requested_base_count(), 3);
        assert_eq!(
            registry.found_record_count(),
            4,
            "expected R1+R2 for READ_A and READ_C"
        );
        assert_eq!(
            registry.missing_base_count(),
            1,
            "READ_MISSING should be missing"
        );
    }

    #[test]
    fn inspect_fc_parsing() {
        use rinfuse_io::fc_intermediates::FcOutputDir;
        let fixtures = fixtures_dir().join("fusioncatcher_minimal");

        let fc_dir = FcOutputDir::detect(&fixtures, false, 3).unwrap();
        assert_eq!(fc_dir.candidate_reports.len(), 1);
        assert_eq!(fc_dir.supporting_read_files.len(), 1);

        let candidates = fc_dir.parse_all_candidates().unwrap();
        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[0].gene_5p, "BCR");
        assert_eq!(candidates[0].gene_3p, "ABL1");

        let tokens = fc_dir.collect_read_id_tokens().unwrap();
        // Fixture has READ_A/1, READ_B/2, READ_MISSING/1
        assert_eq!(tokens.len(), 3);
        assert!(tokens.contains(&"READ_A".to_string()));
        assert!(tokens.contains(&"READ_B".to_string()));
        assert!(tokens.contains(&"READ_MISSING".to_string()));
    }

    #[test]
    fn inspect_fc_recursive_discovery() {
        use rinfuse_io::fc_intermediates::FcOutputDir;
        let fixtures = fixtures_dir().join("fusioncatcher_minimal");
        // Create a nested structure for testing
        let temp = tempfile::tempdir().unwrap();
        let nested = temp.path().join("subdir");
        std::fs::create_dir_all(&nested).unwrap();
        std::fs::copy(
            fixtures.join("final-list_candidate-fusion-genes.txt"),
            nested.join("final-list_candidate-fusion-genes.txt"),
        )
        .unwrap();

        let fc_dir = FcOutputDir::detect(temp.path(), true, 3).unwrap();
        assert_eq!(fc_dir.candidate_reports.len(), 1);
        assert!(fc_dir.candidate_reports[0]
            .to_string_lossy()
            .contains("subdir/final-list_candidate-fusion-genes.txt"));
    }
}
