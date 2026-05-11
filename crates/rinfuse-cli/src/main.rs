mod args;
mod commands;

use anyhow::Result;
use args::{Cli, Commands};
use clap::Parser;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::ExtractReads(args) => commands::extract_reads::run(args),
        Commands::InspectFc(args) => commands::inspect_fc::run(args),
        Commands::Compare(args) => commands::compare::run(args),
    }
}
