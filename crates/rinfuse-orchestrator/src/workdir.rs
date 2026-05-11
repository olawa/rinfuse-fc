use anyhow::{Context, Result};
use std::fs;
use std::path::{Path, PathBuf};

pub struct WorkDir {
    pub root: PathBuf,
}

impl WorkDir {
    pub fn new(root: &Path) -> Result<Self> {
        if !root.exists() {
            fs::create_dir_all(root)
                .with_context(|| format!("failed to create workdir {}", root.display()))?;
        }
        Ok(Self {
            root: root.to_path_buf(),
        })
    }

    pub fn log_dir(&self) -> PathBuf {
        self.root.join("logs")
    }

    pub fn init(&self) -> Result<()> {
        let ld = self.log_dir();
        if !ld.exists() {
            fs::create_dir_all(&ld)?;
        }
        Ok(())
    }
}
