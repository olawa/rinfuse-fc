use crate::args::InspectFcArgs;
use anyhow::{bail, Result};

pub fn run(_args: InspectFcArgs) -> Result<()> {
    bail!("inspect-fc is not implemented yet. Start with extract-reads and read-name normalization.")
}
