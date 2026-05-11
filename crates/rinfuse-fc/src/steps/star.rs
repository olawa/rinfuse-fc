use rinfuse_orchestrator::CommandRunner;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Debug, Serialize, Deserialize, Clone, Default)]
pub struct StarOutputs {
    pub chimeric_junction: Option<PathBuf>,
    pub chimeric_sam: Option<PathBuf>,
    pub aligned_bam: Option<PathBuf>,
    pub log_final: Option<PathBuf>,
}

pub struct StarStep {
    pub bin: String,
    pub reads: Vec<PathBuf>,
    pub genome_dir: PathBuf,
    pub out_dir: PathBuf,
    pub threads: u32,
    pub extra_args: Vec<String>,
}

impl StarStep {
    pub fn new(bin: &str, genome_dir: &Path, out_dir: &Path) -> Self {
        Self {
            bin: bin.to_string(),
            reads: Vec::new(),
            genome_dir: genome_dir.to_path_buf(),
            out_dir: out_dir.to_path_buf(),
            threads: 1,
            extra_args: Vec::new(),
        }
    }

    pub fn build_command(&self, work_dir: &Path) -> CommandRunner {
        let mut args = vec![
            "--genomeDir".to_string(),
            self.genome_dir.display().to_string(),
            "--readFilesIn".to_string(),
        ];

        for r in &self.reads {
            args.push(r.display().to_string());
        }

        args.extend(vec![
            "--runThreadN".to_string(),
            self.threads.to_string(),
            "--outFileNamePrefix".to_string(),
            self.out_dir.join("star/").display().to_string(),
            "--chimSegmentMin".to_string(),
            "12".to_string(),
            "--chimJunctionOverhangMin".to_string(),
            "12".to_string(),
            "--chimOutType".to_string(),
            "Junctions".to_string(),
            "SeparateSAMold".to_string(),
            "--outSAMtype".to_string(),
            "BAM".to_string(),
            "Unsorted".to_string(),
        ]);

        if self.reads.iter().any(|r| r.extension().is_some_and(|e| e == "gz")) {
            args.push("--readFilesCommand".to_string());
            args.push("zcat".to_string());
        }

        for extra in &self.extra_args {
            args.push(extra.clone());
        }

        let arg_refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
        CommandRunner::new(&self.bin, &arg_refs, work_dir).with_label("star_chimeric")
    }

    pub fn discover_outputs(&self) -> StarOutputs {
        let star_out = self.out_dir.join("star");
        let mut outputs = StarOutputs::default();

        let candidates = [
            ("Chimeric.out.junction", &mut outputs.chimeric_junction),
            ("Chimeric.out.sam", &mut outputs.chimeric_sam),
            ("Aligned.out.bam", &mut outputs.aligned_bam),
            ("Log.final.out", &mut outputs.log_final),
        ];

        for (name, target) in candidates {
            let p = star_out.join(name);
            if p.exists() {
                *target = Some(p);
            }
        }

        outputs
    }
}
