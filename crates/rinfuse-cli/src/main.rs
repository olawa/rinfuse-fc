use anyhow::Result;
use clap::Parser;
use rinfuse_cli::args::{Cli, Commands};
use rinfuse_cli::commands;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::ExtractReads(args) => commands::extract_reads::run(args),
        Commands::InspectFc(args) => commands::inspect_fc::run(args),
        Commands::Compare(args) => commands::compare::run(args),
        Commands::RunCommand(args) => commands::run_command::run(args),
        Commands::RunStar(args) => commands::run_star::run(args),
        Commands::ParseStar(args) => commands::parse_star::run(args),
    }
}
