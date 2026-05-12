use crate::args::RunCommandArgs;
use anyhow::{anyhow, Result};
use rinfuse_orchestrator::{
    manifest::{StepManifest, StepStatus},
    timing,
    workdir::WorkDir,
    CommandRunner, OrchestratorError, RunManifest,
};
use std::fs;
use std::io::BufWriter;

pub fn run(args: RunCommandArgs) -> Result<()> {
    let wd = WorkDir::new(&args.out)?;
    wd.init()?;

    let runner = CommandRunner::new(
        &args.program,
        &args.arg.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
        &wd.root,
    )
    .with_log_dir(&wd.log_dir())
    .with_dry_run(args.dry_run);

    let start_time = chrono::Utc::now();
    let res_result = runner.run();
    let end_time = chrono::Utc::now();

    let mut command_failure = None;
    let (res, status) = match res_result {
        Ok(r) => (r, StepStatus::Completed),
        Err(OrchestratorError::CommandFailed {
            program,
            exit_code,
            result,
        }) => {
            command_failure = Some(anyhow!(
                "command '{}' failed with exit code {:?}",
                program,
                exit_code
            ));
            (*result, StepStatus::Failed)
        }
        Err(e) => return Err(e.into()),
    };

    // Create manifest
    let mut manifest = RunManifest {
        workflow_name: "run-command-test".to_string(),
        start_time: Some(start_time),
        end_time: Some(end_time),
        total_walltime_ms: res.walltime_ms,
        ..Default::default()
    };

    let step = StepManifest {
        name: "execute".to_string(),
        commands: vec![res.clone()],
        status,
    };
    manifest.steps.push(step);

    // Write manifest.json
    let mut mw = BufWriter::new(fs::File::create(args.out.join("manifest.json"))?);
    serde_json::to_writer_pretty(&mut mw, &manifest)?;

    // Write timings.tsv
    timing::write_timings_tsv(&args.out.join("timings.tsv"), &[res])?;

    if let Some(err) = command_failure {
        return Err(err);
    }

    Ok(())
}
