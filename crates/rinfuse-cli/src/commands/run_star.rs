use crate::args::RunStarArgs;
use anyhow::Result;
use rinfuse_fc::steps::star::StarStep;
use rinfuse_io::star::{
    parse_chimeric_junctions, write_junctions_jsonl, write_junctions_tsv, write_parse_summary,
};
use rinfuse_fc::steps::aggregate::aggregate_star_junctions;
use crate::commands::aggregate_star::write_outputs as write_candidates;
use rinfuse_orchestrator::{
    manifest::{StepManifest, StepStatus},
    timing, OrchestratorError, RunManifest,
    workdir::WorkDir,
};
use std::fs;
use std::io::BufWriter;

pub fn run(args: RunStarArgs) -> Result<()> {
    let wd = WorkDir::new(&args.out)?;
    wd.init()?;

    let star_dir = args.out.join("star");
    if !star_dir.exists() {
        fs::create_dir_all(&star_dir)?;
    }

    let mut star_step = StarStep::new(&args.star_bin, &args.star_index, &args.out);
    star_step.reads = args.reads;
    star_step.threads = args.threads;
    star_step.extra_args = args.extra_star_arg;

    let runner = star_step
        .build_command(&wd.root)
        .with_log_dir(&wd.log_dir())
        .with_dry_run(args.dry_run);

    let start_time = chrono::Utc::now();
    let res_result = runner.run();
    let end_time = chrono::Utc::now();

    let (res, status) = match res_result {
        Ok(r) => (r, StepStatus::Completed),
        Err(OrchestratorError::CommandFailed { result, .. }) => (*result, StepStatus::Failed),
        Err(e) => return Err(e.into()),
    };

    let outputs = star_step.discover_outputs();

    let mut manifest = RunManifest {
        workflow_name: "star-run".to_string(),
        start_time: Some(start_time),
        end_time: Some(end_time),
        total_walltime_ms: res.walltime_ms,
        ..Default::default()
    };

    let step = StepManifest {
        name: "star_alignment".to_string(),
        commands: vec![res.clone()],
        status: status.clone(),
    };
    manifest.steps.push(step);

    let mut mw = BufWriter::new(fs::File::create(args.out.join("manifest.json"))?);
    serde_json::to_writer_pretty(&mut mw, &manifest)?;

    timing::write_timings_tsv(&args.out.join("timings.tsv"), &[res])?;

    let mut ow = BufWriter::new(fs::File::create(args.out.join("star_outputs.json"))?);
    serde_json::to_writer_pretty(&mut ow, &outputs)?;

    // Optionally parse Chimeric.out.junction
    if args.parse {
        match &outputs.chimeric_junction {
            Some(junction_path) => {
                let evidence_dir = args.out.join("evidence");
                if !evidence_dir.exists() {
                    fs::create_dir_all(&evidence_dir)?;
                }
                let (junctions, report) = parse_chimeric_junctions(junction_path)?;
                write_junctions_jsonl(&evidence_dir.join("star_junctions.jsonl"), &junctions)?;
                write_junctions_tsv(&evidence_dir.join("star_junctions.tsv"), &junctions)?;
                write_parse_summary(&evidence_dir.join("star_parse_summary.json"), &report)?;
                eprintln!(
                    "Parsed {} junctions ({} warnings) → {}",
                    report.parsed_ok,
                    report.parse_warnings.len(),
                    evidence_dir.display()
                );

                if let Some(genes_path) = &args.genes {
                    let candidates_dir = args.out.join("candidates");
                    if !candidates_dir.exists() {
                        fs::create_dir_all(&candidates_dir)?;
                    }
                    let candidates = aggregate_star_junctions(junction_path, genes_path)?;
                    write_candidates(&candidates_dir, "star", &candidates)?;
                    eprintln!(
                        "Aggregated {} candidates → {}",
                        candidates.len(),
                        candidates_dir.display()
                    );
                }
            }
            None => {
                eprintln!(
                    "{}--parse requested but Chimeric.out.junction not found (dry-run={}).",
                    if args.dry_run { "[dry-run] " } else { "" },
                    args.dry_run
                );
            }
        }
    }

    eprintln!("STAR run complete. Status: {:?}", status);

    Ok(())
}
