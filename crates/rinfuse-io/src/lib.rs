pub mod fastq;
pub mod read_registry;

pub use fastq::{open_maybe_gz, FastqRecord, FastqReader};
pub use read_registry::{ReadRegistry, ReadRegistryBuildOptions};
