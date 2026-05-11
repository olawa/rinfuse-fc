use crate::args::CompareArgs;
use anyhow::{bail, Result};

pub fn run(_args: CompareArgs) -> Result<()> {
    bail!("compare is not implemented yet. It will compare FusionCatcher and rinfuse-fc candidate outputs.")
}
