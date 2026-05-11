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
        Self {
            include_mates: true,
        }
    }
}

#[derive(Debug, Default)]
pub struct ReadRegistry {
    /// Maps base name to a set of requested mates. None means "any/all".
    requested: HashMap<String, HashSet<Option<MateSide>>>,
    found: HashMap<(String, MateSide), FastqRecord>,
    missing_requested_bases: HashSet<String>,
}

impl ReadRegistry {
    pub fn from_read_id_file(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("failed to open read id file {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut requested: HashMap<String, HashSet<Option<MateSide>>> = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            let normalized = normalize_read_name(trimmed);
            let entry = requested.entry(normalized.base).or_default();
            if normalized.mate == MateSide::Unknown {
                entry.insert(None);
            } else {
                entry.insert(Some(normalized.mate));
            }
        }

        let missing_requested_bases = requested.keys().cloned().collect();
        Ok(Self {
            requested,
            found: HashMap::new(),
            missing_requested_bases,
        })
    }

    /// Build a registry from already-normalized base names (e.g. from fc_intermediates).
    pub fn from_tokens(bases: &[String]) -> Self {
        let mut requested: HashMap<String, HashSet<Option<MateSide>>> = HashMap::new();
        for base in bases {
            if !base.is_empty() {
                requested.entry(base.clone()).or_default().insert(None);
            }
        }
        let missing_requested_bases = requested.keys().cloned().collect();
        Self {
            requested,
            found: HashMap::new(),
            missing_requested_bases,
        }
    }

    pub fn collect_from_fastq_paths(
        &mut self,
        read_paths: &[PathBuf],
        opts: &ReadRegistryBuildOptions,
    ) -> Result<()> {
        for path in read_paths {
            self.collect_from_fastq_path(path, opts)?;
        }
        Ok(())
    }

    fn collect_from_fastq_path(
        &mut self,
        path: &Path,
        opts: &ReadRegistryBuildOptions,
    ) -> Result<()> {
        let mut reader = FastqReader::from_path(path)?;
        while let Some(record) = reader.next_record()? {
            let normalized = normalize_read_name(&record.header);
            if let Some(desired_mates) = self.requested.get(&normalized.base) {
                let should_include = if opts.include_mates {
                    true
                } else {
                    // If include_mates is false, only include if explicit mate matches or was generic
                    desired_mates.contains(&None) || desired_mates.contains(&Some(normalized.mate))
                };

                if should_include {
                    self.missing_requested_bases.remove(&normalized.base);
                    self.found
                        .insert((normalized.base, normalized.mate), record);
                }
            }
        }
        Ok(())
    }

    pub fn write_fastq(&self, out: &Path) -> Result<()> {
        let file =
            File::create(out).with_context(|| format!("failed to create {}", out.display()))?;
        let mut writer = BufWriter::new(file);

        let mut keys: Vec<_> = self.found.keys().cloned().collect();
        keys.sort_by(|a, b| {
            a.0.cmp(&b.0)
                .then(format!("{:?}", a.1).cmp(&format!("{:?}", b.1)))
        });

        for key in keys {
            if let Some(record) = self.found.get(&key) {
                record.write_to(&mut writer)?;
            }
        }

        writer.flush()?;
        Ok(())
    }

    pub fn write_missing(&self, out: &Path) -> Result<()> {
        let file =
            File::create(out).with_context(|| format!("failed to create {}", out.display()))?;
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
        self.requested.len()
    }

    pub fn missing_base_count(&self) -> usize {
        self.missing_requested_bases.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_include_mates_true() {
        let mut reg = ReadRegistry::from_tokens(&["READ_A".to_string()]);
        let temp = tempdir().unwrap();
        let fq = temp.path().join("test.fq");
        fs::write(&fq, "@READ_A/1\nACGT\n+\n####\n@READ_A/2\nTGCA\n+\n####\n").unwrap();

        let opts = ReadRegistryBuildOptions {
            include_mates: true,
        };
        reg.collect_from_fastq_paths(&[fq], &opts).unwrap();

        assert_eq!(reg.found_record_count(), 2);
    }

    #[test]
    fn test_include_mates_false_honors_explicit() {
        // Request specifically /1
        let temp = tempdir().unwrap();
        let id_file = temp.path().join("ids.txt");
        fs::write(&id_file, "READ_A/1\n").unwrap();

        let mut reg = ReadRegistry::from_read_id_file(&id_file).unwrap();

        let fq = temp.path().join("test.fq");
        fs::write(&fq, "@READ_A/1\nACGT\n+\n####\n@READ_A/2\nTGCA\n+\n####\n").unwrap();

        let opts = ReadRegistryBuildOptions {
            include_mates: false,
        };
        reg.collect_from_fastq_paths(&[fq], &opts).unwrap();

        // Should only have 1 (the /1 mate)
        assert_eq!(reg.found_record_count(), 1);
        assert!(reg
            .found
            .contains_key(&("READ_A".to_string(), MateSide::R1)));
        assert!(!reg
            .found
            .contains_key(&("READ_A".to_string(), MateSide::R2)));
    }

    #[test]
    fn test_include_mates_false_generic_still_finds_both() {
        let mut reg = ReadRegistry::from_tokens(&["READ_A".to_string()]);
        let temp = tempdir().unwrap();
        let fq = temp.path().join("test.fq");
        fs::write(&fq, "@READ_A/1\nACGT\n+\n####\n@READ_A/2\nTGCA\n+\n####\n").unwrap();

        let opts = ReadRegistryBuildOptions {
            include_mates: false,
        };
        reg.collect_from_fastq_paths(&[fq], &opts).unwrap();

        assert_eq!(reg.found_record_count(), 2);
    }
}
