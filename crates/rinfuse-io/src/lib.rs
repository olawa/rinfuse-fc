pub mod fastq;
pub mod fc_intermediates;
pub mod read_registry;

pub use fastq::{open_maybe_gz, FastqReader, FastqRecord};
pub use read_registry::{ReadRegistry, ReadRegistryBuildOptions};
