use crate::command::CommandResult;
use anyhow::Result;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn write_timings_tsv(path: &Path, results: &[CommandResult]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "program\twalltime_ms\tstart_time\texit_code\tdry_run"
    )?;
    for res in results {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            res.program,
            res.walltime_ms,
            res.start_time,
            res.exit_code
                .map(|c| c.to_string())
                .unwrap_or_else(|| "-".to_string()),
            res.dry_run
        )?;
    }
    writer.flush()?;
    Ok(())
}
