use crate::args::ExtractReadsArgs;
use anyhow::Result;
use rinfuse_io::{ReadRegistry, ReadRegistryBuildOptions};

pub fn run(args: ExtractReadsArgs) -> Result<()> {
    let mut registry = ReadRegistry::from_read_id_file(&args.read_ids)?;
    let opts = ReadRegistryBuildOptions::default();

    registry.collect_from_fastq_paths(&args.reads, &opts)?;
    registry.write_fastq(&args.out)?;

    if let Some(missing_out) = args.missing_out.as_deref() {
        registry.write_missing(missing_out)?;
    }

    eprintln!(
        "requested_bases={} found_records={} missing_bases={}",
        registry.requested_base_count(),
        registry.found_record_count(),
        registry.missing_base_count()
    );

    Ok(())
}
