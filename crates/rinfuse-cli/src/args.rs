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
    /// FusionCatcher output directory to inspect.
    #[arg(long)]
    pub fc_out: PathBuf,

    /// Comma-separated FASTQ paths for read extraction.
    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,

    /// Output directory for inspected results.
    #[arg(long)]
    pub out: PathBuf,

    /// Scan fc-out recursively for reports and read IDs.
    #[arg(long, default_value_t = false)]
    pub recursive: bool,

    /// Maximum depth for recursive scanning.
    #[arg(long, default_value_t = 3)]
    pub max_depth: usize,

    /// Focus on specific genes (repeatable or comma-separated).
    #[arg(long, value_delimiter = ',')]
    pub focus_gene: Vec<String>,
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
