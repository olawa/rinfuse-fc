#[cfg(test)]
mod extract_reads_integration {
    use rinfuse_io::{ReadRegistry, ReadRegistryBuildOptions};
    use std::path::PathBuf;

    fn fixtures_dir() -> PathBuf {
        // Relative to workspace root
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .parent()  // crates/
            .unwrap()
            .parent()  // workspace root
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
        assert_eq!(registry.found_record_count(), 4, "expected R1+R2 for READ_A and READ_C");
        assert_eq!(registry.missing_base_count(), 1, "READ_MISSING should be missing");
    }
}
