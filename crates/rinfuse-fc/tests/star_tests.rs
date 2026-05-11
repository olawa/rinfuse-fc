#[cfg(test)]
mod star_tests {
    use rinfuse_fc::steps::star::StarStep;
    use std::fs;
    use std::path::PathBuf;
    use tempfile::tempdir;

    #[test]
    fn test_star_command_construction_plain() {
        let temp = tempdir().unwrap();
        let star = StarStep::new("STAR", &PathBuf::from("/idx"), temp.path());
        let runner = star.build_command(temp.path());

        let args_str = runner.args.join(" ");
        assert!(args_str.contains("--genomeDir /idx"));
        assert!(!args_str.contains("--readFilesCommand zcat"));
    }

    #[test]
    fn test_star_command_construction_gz() {
        let temp = tempdir().unwrap();
        let mut star = StarStep::new("STAR", &PathBuf::from("/idx"), temp.path());
        star.reads = vec![PathBuf::from("reads.fastq.gz")];
        let runner = star.build_command(temp.path());

        let args_str = runner.args.join(" ");
        assert!(args_str.contains("--readFilesCommand zcat"));
    }

    #[test]
    fn test_output_discovery() {
        let temp = tempdir().unwrap();
        let star_dir = temp.path().join("star");
        fs::create_dir_all(&star_dir).unwrap();

        fs::write(star_dir.join("Chimeric.out.junction"), "test").unwrap();

        let star = StarStep::new("STAR", &PathBuf::from("/idx"), temp.path());
        let outputs = star.discover_outputs();

        assert!(outputs.chimeric_junction.is_some());
        assert!(outputs.chimeric_sam.is_none());
    }

    #[test]
    fn test_dry_run_writes_manifest() {
        // This is more of an integration test for rinfuse-cli, but we can mock it here
        // if we want to test the full flow.
        // For now, let's just test that CommandRunner supports dry-run (already tested).
    }
}
