use rinfuse_cli::args::{RunCommandArgs, RunStarArgs};
use std::fs;
use tempfile::tempdir;

#[test]
fn run_command_false_returns_err_and_writes_failed_manifest() {
    let temp = tempdir().unwrap();
    let out = temp.path().join("run_command_false");

    let result = rinfuse_cli::commands::run_command::run(RunCommandArgs {
        program: "false".to_string(),
        arg: vec![],
        out: out.clone(),
        dry_run: false,
    });

    assert!(result.is_err(), "failing command should return Err");
    assert!(out.join("manifest.json").exists());
    assert!(out.join("timings.tsv").exists());
    assert!(out.join("logs/false.stdout.log").exists());
    assert!(out.join("logs/false.stderr.log").exists());

    let manifest: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["steps"][0]["status"], "Failed");
    assert_eq!(manifest["steps"][0]["commands"][0]["exit_code"], 1);
}

#[test]
fn run_star_false_returns_err_and_writes_failed_manifest() {
    let temp = tempdir().unwrap();
    let out = temp.path().join("run_star_false");
    let star_index = temp.path().join("star_index");
    fs::create_dir_all(&star_index).unwrap();

    let result = rinfuse_cli::commands::run_star::run(RunStarArgs {
        reads: vec![],
        star_index,
        out: out.clone(),
        threads: 1,
        star_bin: "false".to_string(),
        dry_run: false,
        extra_star_arg: vec![],
        parse: false,
        genes: None,
    });

    assert!(result.is_err(), "failing STAR command should return Err");
    assert!(out.join("manifest.json").exists());
    assert!(out.join("timings.tsv").exists());
    assert!(out.join("star_outputs.json").exists());
    assert!(out.join("logs/star_chimeric.stdout.log").exists());
    assert!(out.join("logs/star_chimeric.stderr.log").exists());

    let manifest: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(out.join("manifest.json")).unwrap()).unwrap();
    assert_eq!(manifest["steps"][0]["status"], "Failed");
    assert_eq!(manifest["steps"][0]["commands"][0]["exit_code"], 1);
}
