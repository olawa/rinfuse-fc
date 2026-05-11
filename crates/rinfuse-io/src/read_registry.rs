use crate::fastq::{FastqReader, FastqRecord};
use anyhow::{Context, Result};
use rinfuse_core::{normalize_read_name, MateSide};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone)]
pub struct ReadRegistryBuildOptions {
    /// Also emit the mate if one mate is requested.
    pub include_mates: bool,
}

impl Default for ReadRegistryBuildOptions {
    fn default() -> Self {
        Self { include_mates: true }
    }
}

#[derive(Debug, Default)]
pub struct ReadRegistry {
    requested_bases: HashSet<String>,
    found: HashMap<(String, MateSide), FastqRecord>,
    missing_requested_bases: HashSet<String>,
}

impl ReadRegistry {
    pub fn from_read_id_file(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("failed to open read id file {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut requested_bases = HashSet::new();

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let normalized = normalize_read_name(trimmed);
            requested_bases.insert(normalized.base);
        }

        let missing_requested_bases = requested_bases.clone();
        Ok(Self {
            requested_bases,
            found: HashMap::new(),
            missing_requested_bases,
        })
    }

    pub fn collect_from_fastq_paths(&mut self, read_paths: &[PathBuf], _opts: &ReadRegistryBuildOptions) -> Result<()> {
        for path in read_paths {
            self.collect_from_fastq_path(path)?;
        }
        Ok(())
    }

    fn collect_from_fastq_path(&mut self, path: &Path) -> Result<()> {
        let mut reader = FastqReader::from_path(path)?;
        while let Some(record) = reader.next_record()? {
            let normalized = normalize_read_name(&record.header);
            if self.requested_bases.contains(&normalized.base) {
                self.missing_requested_bases.remove(&normalized.base);
                self.found.insert((normalized.base, normalized.mate), record);
            }
        }
        Ok(())
    }

    pub fn write_fastq(&self, out: &Path) -> Result<()> {
        let file = File::create(out).with_context(|| format!("failed to create {}", out.display()))?;
        let mut writer = BufWriter::new(file);

        let mut keys: Vec<_> = self.found.keys().cloned().collect();
        keys.sort_by(|a, b| a.0.cmp(&b.0).then(format!("{:?}", a.1).cmp(&format!("{:?}", b.1))));

        for key in keys {
            if let Some(record) = self.found.get(&key) {
                record.write_to(&mut writer)?;
            }
        }

        writer.flush()?;
        Ok(())
    }

    pub fn write_missing(&self, out: &Path) -> Result<()> {
        let file = File::create(out).with_context(|| format!("failed to create {}", out.display()))?;
        let mut writer = BufWriter::new(file);
        let mut missing: Vec<_> = self.missing_requested_bases.iter().collect();
        missing.sort();
        for id in missing {
            writeln!(writer, "{}", id)?;
        }
        Ok(())
    }

    pub fn found_record_count(&self) -> usize {
        self.found.len()
    }

    pub fn requested_base_count(&self) -> usize {
        self.requested_bases.len()
    }

    pub fn missing_base_count(&self) -> usize {
        self.missing_requested_bases.len()
    }
}
