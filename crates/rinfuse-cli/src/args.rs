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
    /// Extract reads from FASTQ based on read IDs.
    ExtractReads(ExtractReadsArgs),

    /// Inspect a FusionCatcher output directory.
    InspectFc(InspectFcArgs),

    /// Compare candidate fusions between two reports.
    Compare(CompareArgs),

    /// [DEV/TEST] Run an external command and record manifest/logs.
    RunCommand(RunCommandArgs),

    /// Run STAR aligner to collect chimeric reads.
    RunStar(RunStarArgs),

    /// Parse a STAR Chimeric.out.junction file into typed evidence.
    ParseStar(ParseStarArgs),

    /// Aggregate STAR junctions into fusion candidates using gene annotations.
    AggregateStar(AggregateStarArgs),

    /// Validate sample by comparing FusionCatcher outputs against rinfuse-fc STAR candidates.
    ValidateSample(ValidateSampleArgs),

    /// Summarize multiple sample validation reports into a cohort-level report.
    ValidateCohort(ValidateCohortArgs),
}

#[derive(Debug, Args)]
pub struct ValidateCohortArgs {
    /// Validation directories to aggregate (repeatable).
    #[arg(long)]
    pub validation_dir: Vec<PathBuf>,

    /// Output directory for cohort reports.
    #[arg(long)]
    pub out: PathBuf,

    /// Optional focus genes (repeatable or comma-separated).
    #[arg(long, value_delimiter = ',')]
    pub focus_gene: Vec<String>,
}

#[derive(Debug, Args)]
pub struct ValidateSampleArgs {
    /// FusionCatcher output directory to inspect.
    #[arg(long)]
    pub fc_out: PathBuf,

    /// Path to parsed STAR candidates TSV or JSONL file.
    #[arg(long)]
    pub star_candidates: PathBuf,

    /// Optional comma-separated FASTQ paths to extract supporting reads.
    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,

    /// Output directory for validation summary and reports.
    #[arg(long)]
    pub out: PathBuf,

    /// Focus on specific genes (repeatable or comma-separated).
    #[arg(long, value_delimiter = ',')]
    pub focus_gene: Vec<String>,
}

#[derive(Debug, Args)]
pub struct AggregateStarArgs {
    /// Path to parsed STAR junctions JSONL or raw Chimeric.out.junction.
    #[arg(long)]
    pub junctions: PathBuf,

    /// Path to gene intervals TSV.
    #[arg(long)]
    pub genes: PathBuf,

    /// Output directory for candidate TSV/JSONL.
    #[arg(long)]
    pub out: PathBuf,

    /// Optional specific genes to focus on (comma-separated or repeatable).
    #[arg(long, value_delimiter = ',')]
    pub focus_gene: Vec<String>,
}

#[derive(Debug, Args)]
pub struct ParseStarArgs {
    /// Path to STAR Chimeric.out.junction file.
    #[arg(long)]
    pub junction: PathBuf,

    /// Output directory for junction JSONL, TSV and parse summary.
    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, Args)]
pub struct RunStarArgs {
    /// Input FASTQ files (comma-separated or repeatable).
    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,

    /// Path to STAR genome index directory.
    #[arg(long)]
    pub star_index: PathBuf,

    /// Output directory.
    #[arg(long)]
    pub out: PathBuf,

    /// Number of threads for STAR.
    #[arg(long, default_value_t = 1)]
    pub threads: u32,

    /// Path to STAR binary.
    #[arg(long, default_value = "STAR")]
    pub star_bin: String,

    /// Dry-run mode.
    #[arg(long, default_value_t = false)]
    pub dry_run: bool,

    /// Extra arguments for STAR.
    #[arg(long)]
    pub extra_star_arg: Vec<String>,

    /// Parse Chimeric.out.junction after alignment and write evidence files.
    #[arg(long, default_value_t = false)]
    pub parse: bool,

    /// Optional path to gene intervals TSV. If provided with --parse, will automatically aggregate candidates.
    #[arg(long)]
    pub genes: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct RunCommandArgs {
    /// Program to execute.
    #[arg(long)]
    pub program: String,

    /// Arguments for the program (repeatable).
    #[arg(long)]
    pub arg: Vec<String>,

    /// Output directory for manifest and logs.
    #[arg(long)]
    pub out: PathBuf,

    /// Dry-run mode (do not execute).
    #[arg(long, default_value_t = false)]
    pub dry_run: bool,
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
    /// FusionCatcher output directory or report file.
    #[arg(long)]
    pub fc: PathBuf,

    /// rinfuse-fc output directory or candidates.tsv.
    #[arg(long)]
    pub rs: PathBuf,

    /// Output comparison TSV path.
    #[arg(long)]
    pub out: PathBuf,

    /// Focus on specific genes (repeatable or comma-separated).
    #[arg(long, value_delimiter = ',')]
    pub focus_gene: Vec<String>,
}
