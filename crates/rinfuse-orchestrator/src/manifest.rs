use crate::command::CommandResult;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Debug, Serialize, Deserialize, Clone, Default)]
pub struct RunManifest {
    pub workflow_name: String,
    pub start_time: Option<chrono::DateTime<chrono::Utc>>,
    pub end_time: Option<chrono::DateTime<chrono::Utc>>,
    pub total_walltime_ms: u128,
    pub steps: Vec<StepManifest>,
    pub tool_versions: Vec<ToolVersion>,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct StepManifest {
    pub name: String,
    pub commands: Vec<CommandResult>,
    pub status: StepStatus,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub enum StepStatus {
    Pending,
    Running,
    Completed,
    Failed,
    Skipped,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ToolVersion {
    pub name: String,
    pub version: String,
}

#[derive(Debug, Serialize, Deserialize, Clone, Default)]
pub struct InputManifest {
    pub files: Vec<PathBuf>,
}
