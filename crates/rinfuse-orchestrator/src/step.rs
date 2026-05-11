use crate::command::{CommandResult, CommandRunner};
use crate::manifest::StepStatus;

pub struct Step {
    pub name: String,
    pub commands: Vec<CommandRunner>,
}

impl Step {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            commands: Vec::new(),
        }
    }

    pub fn add_command(&mut self, cmd: CommandRunner) {
        self.commands.push(cmd);
    }
}

pub struct StepResult {
    pub name: String,
    pub command_results: Vec<CommandResult>,
    pub status: StepStatus,
}
