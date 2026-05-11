use crate::args::ParseStarArgs;
use anyhow::Result;
use rinfuse_io::star::{
    parse_chimeric_junctions, write_junctions_jsonl, write_junctions_tsv, write_parse_summary,
};
use std::fs;

pub fn run(args: ParseStarArgs) -> Result<()> {
    if !args.out.exists() {
        fs::create_dir_all(&args.out)?;
    }

    let (junctions, report) = parse_chimeric_junctions(&args.junction)?;

    write_junctions_jsonl(&args.out.join("star_junctions.jsonl"), &junctions)?;
    write_junctions_tsv(&args.out.join("star_junctions.tsv"), &junctions)?;
    write_parse_summary(&args.out.join("star_parse_summary.json"), &report)?;

    eprintln!(
        "Parsed {} junctions from {} ({} warnings)",
        report.parsed_ok,
        args.junction.display(),
        report.parse_warnings.len()
    );

    Ok(())
}
