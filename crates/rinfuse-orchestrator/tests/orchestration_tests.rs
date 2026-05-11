#[cfg(test)]
mod orchestration_tests {
    use rinfuse_orchestrator::{timing, workdir::WorkDir, CommandRunner, RunManifest};
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_successful_command() {
        let temp = tempdir().unwrap();
        let wd = WorkDir::new(temp.path()).unwrap();
        wd.init().unwrap();

        let runner = CommandRunner::new("echo", &["hello"], temp.path())
            .with_label("test_echo")
            .with_log_dir(&wd.log_dir());

        let res = runner.run().unwrap();
        assert_eq!(res.exit_code, Some(0));
        assert!(wd.log_dir().join("test_echo.stdout.log").exists());
        let stdout = fs::read_to_string(wd.log_dir().join("test_echo.stdout.log")).unwrap();
        assert_eq!(stdout.trim(), "hello");
    }

    #[test]
    fn test_failing_command_default_policy() {
        let temp = tempdir().unwrap();
        let wd = WorkDir::new(temp.path()).unwrap();
        wd.init().unwrap();

        let runner = CommandRunner::new("false", &[], temp.path()).with_log_dir(&wd.log_dir());

        let res = runner.run();
        assert!(res.is_err());
        if let Err(e) = res {
            assert!(e.to_string().contains("failed with exit code Some(1)"));
        }
    }

    #[test]
    fn test_failing_command_allow_nonzero() {
        let temp = tempdir().unwrap();
        let wd = WorkDir::new(temp.path()).unwrap();
        wd.init().unwrap();

        let runner = CommandRunner::new("false", &[], temp.path()).with_allow_nonzero(true);

        let res = runner.run().unwrap();
        assert_eq!(res.exit_code, Some(1));
    }

    #[test]
    fn test_labels_prevent_overwrite() {
        let temp = tempdir().unwrap();
        let wd = WorkDir::new(temp.path()).unwrap();
        wd.init().unwrap();

        let r1 = CommandRunner::new("echo", &["A"], temp.path())
            .with_label("cmd1")
            .with_log_dir(&wd.log_dir());
        let r2 = CommandRunner::new("echo", &["B"], temp.path())
            .with_label("cmd2")
            .with_log_dir(&wd.log_dir());

        r1.run().unwrap();
        r2.run().unwrap();

        assert_eq!(
            fs::read_to_string(wd.log_dir().join("cmd1.stdout.log"))
                .unwrap()
                .trim(),
            "A"
        );
        assert_eq!(
            fs::read_to_string(wd.log_dir().join("cmd2.stdout.log"))
                .unwrap()
                .trim(),
            "B"
        );
    }

    #[test]
    fn test_env_vars() {
        let temp = tempdir().unwrap();
        let wd = WorkDir::new(temp.path()).unwrap();
        wd.init().unwrap();

        let runner = CommandRunner::new("sh", &["-c", "echo $MY_VAR"], temp.path())
            .with_label("test_env")
            .with_env("MY_VAR", "VAL")
            .with_log_dir(&wd.log_dir());

        runner.run().unwrap();
        let stdout = fs::read_to_string(wd.log_dir().join("test_env.stdout.log")).unwrap();
        assert_eq!(stdout.trim(), "VAL");
    }

    #[test]
    fn test_dry_run() {
        let temp = tempdir().unwrap();
        let wd = WorkDir::new(temp.path()).unwrap();
        wd.init().unwrap();

        let runner = CommandRunner::new("echo", &["hello"], temp.path()).with_dry_run(true);

        let res = runner.run().unwrap();
        assert_eq!(res.exit_code, None);
        assert!(res.dry_run);
        // Should not create log files if not executed
        assert!(!wd.log_dir().join("echo.stdout.log").exists());
    }

    #[test]
    fn test_manifest_roundtrip() {
        let manifest = RunManifest {
            workflow_name: "test".to_string(),
            ..Default::default()
        };
        let json = serde_json::to_string(&manifest).unwrap();
        let decoded: RunManifest = serde_json::from_str(&json).unwrap();
        assert_eq!(decoded.workflow_name, "test");
    }

    #[test]
    fn test_timings_creation() {
        let temp = tempdir().unwrap();
        let tsv_path = temp.path().join("timings.tsv");

        let runner = CommandRunner::new("echo", &["hi"], temp.path());
        let res = runner.run().unwrap();

        timing::write_timings_tsv(&tsv_path, &[res]).unwrap();
        assert!(tsv_path.exists());
        let content = fs::read_to_string(tsv_path).unwrap();
        assert!(content.contains("echo"));
        assert!(content.contains("walltime_ms"));
    }
}
