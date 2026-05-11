use anyhow::Result;
use chrono::Utc;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitStatus};
use std::time::Instant;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum OrchestratorError {
    #[error("command '{program}' failed with exit code {exit_code:?}")]
    CommandFailed {
        program: String,
        exit_code: Option<i32>,
        result: Box<CommandResult>,
    },
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Other(#[from] anyhow::Error),
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CommandResult {
    pub label: String,
    pub program: String,
    pub args: Vec<String>,
    pub exit_code: Option<i32>,
    pub start_time: chrono::DateTime<Utc>,
    pub walltime_ms: u128,
    pub stdout_log: Option<PathBuf>,
    pub stderr_log: Option<PathBuf>,
    pub dry_run: bool,
}

pub struct CommandRunner {
    pub label: String,
    pub program: String,
    pub args: Vec<String>,
    pub cwd: PathBuf,
    pub envs: HashMap<String, String>,
    pub dry_run: bool,
    pub allow_nonzero: bool,
    pub log_dir: Option<PathBuf>,
}

impl CommandRunner {
    pub fn new(program: &str, args: &[&str], cwd: &Path) -> Self {
        let base = Path::new(program)
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("cmd");
        Self {
            label: base.to_string(),
            program: program.to_string(),
            args: args.iter().map(|s| s.to_string()).collect(),
            cwd: cwd.to_path_buf(),
            envs: HashMap::new(),
            dry_run: false,
            allow_nonzero: false,
            log_dir: None,
        }
    }

    pub fn with_label(mut self, label: &str) -> Self {
        self.label = label.to_string();
        self
    }

    pub fn with_log_dir(mut self, dir: &Path) -> Self {
        self.log_dir = Some(dir.to_path_buf());
        self
    }

    pub fn with_dry_run(mut self, dry_run: bool) -> Self {
        self.dry_run = dry_run;
        self
    }

    pub fn with_allow_nonzero(mut self, allow: bool) -> Self {
        self.allow_nonzero = allow;
        self
    }

    pub fn with_env(mut self, key: &str, value: &str) -> Self {
        self.envs.insert(key.to_string(), value.to_string());
        self
    }

    pub fn with_envs<I, K, V>(mut self, envs: I) -> Self
    where
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: AsRef<str>,
    {
        for (k, v) in envs {
            self.envs
                .insert(k.as_ref().to_string(), v.as_ref().to_string());
        }
        self
    }

    pub fn run(&self) -> Result<CommandResult, OrchestratorError> {
        let start_utc = Utc::now();
        let start_inst = Instant::now();

        let mut stdout_path = None;
        let mut stderr_path = None;

        if let Some(ref dir) = self.log_dir {
            if !dir.exists() {
                fs::create_dir_all(dir)?;
            }
            stdout_path = Some(dir.join(format!("{}.stdout.log", self.label)));
            stderr_path = Some(dir.join(format!("{}.stderr.log", self.label)));
        }

        let exit_code = if self.dry_run {
            tracing::info!(
                "[DRY-RUN] {} {} (label={})",
                self.program,
                self.args.join(" "),
                self.label
            );
            None
        } else {
            let mut cmd = Command::new(&self.program);
            cmd.args(&self.args).current_dir(&self.cwd);
            for (k, v) in &self.envs {
                cmd.env(k, v);
            }

            if let Some(ref p) = stdout_path {
                let f = fs::File::create(p)?;
                cmd.stdout(f);
            }
            if let Some(ref p) = stderr_path {
                let f = fs::File::create(p)?;
                cmd.stderr(f);
            }

            let status: ExitStatus = cmd.status().map_err(|e| {
                anyhow::anyhow!("failed to execute command '{}': {}", self.program, e)
            })?;

            status.code()
        };

        let walltime_ms = start_inst.elapsed().as_millis();

        let result = CommandResult {
            label: self.label.clone(),
            program: self.program.clone(),
            args: self.args.clone(),
            exit_code,
            start_time: start_utc,
            walltime_ms,
            stdout_log: stdout_path,
            stderr_log: stderr_path,
            dry_run: self.dry_run,
        };

        if !self.allow_nonzero && !self.dry_run {
            if let Some(code) = exit_code {
                if code != 0 {
                    return Err(OrchestratorError::CommandFailed {
                        program: self.program.clone(),
                        exit_code,
                        result: Box::new(result),
                    });
                }
            } else {
                // Terminated by signal
                return Err(OrchestratorError::CommandFailed {
                    program: self.program.clone(),
                    exit_code: None,
                    result: Box::new(result),
                });
            }
        }

        Ok(result)
    }
}
