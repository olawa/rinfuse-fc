use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "rinfuse-fc")]
#[command(about = "Rust FusionCatcher-compatible orchestrator and evidence tools")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Extract read records by read ID from FASTQ files.
    ExtractReads(ExtractReadsArgs),

    /// Inspect an existing FusionCatcher output directory.
    InspectFc(InspectFcArgs),

    /// Compare FusionCatcher and rinfuse-fc outputs.
    Compare(CompareArgs),
}

#[derive(Debug, Args)]
pub struct ExtractReadsArgs {
    /// Comma-separated FASTQ paths, for example R1.fq.gz,R2.fq.gz.
    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,

    /// Text file with one read ID per line.
    #[arg(long)]
    pub read_ids: PathBuf,

    /// Output FASTQ path. First MVP writes plain FASTQ.
    #[arg(long)]
    pub out: PathBuf,

    /// Output missing read IDs here.
    #[arg(long)]
    pub missing_out: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct InspectFcArgs {
    #[arg(long)]
    pub fc_out: PathBuf,

    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,

    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, Args)]
pub struct CompareArgs {
    #[arg(long)]
    pub fc: PathBuf,

    #[arg(long)]
    pub rs: PathBuf,

    #[arg(long)]
    pub out: PathBuf,
}
