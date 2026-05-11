# rinfuse-fc: Rust-orchestrator för FusionCatcher-liknande fusiondetektion

## Syfte

`rinfuse-fc` är tänkt som en modern Rust-baserad orchestrator och evidensmotor för RNA-seq-fusioner. Första målet är inte att ersätta alla algoritmer i FusionCatcher, utan att reproducera FusionCatchers praktiskt viktiga recall med en mer kontrollerad, snabbare och testbar implementation.

Första fasen använder fortfarande externa aligners:

* Bowtie 1
* STAR
* Bowtie2
* BLAT, optional

Men Rust-koden ska äga:

* steg-grafen
* körmanifest
* loggning
* read registry
* parsing av intermediärer
* normaliserad evidence model
* caching/restart
* rapportering
* jämförelse mot original-FusionCatcher

## Problem vi försöker lösa

FusionCatcher har hög sensitivitet, särskilt för svåra fusioner som DUX4/IGH-liknande fall, men implementationen är svår att underhålla och ineffektiv:

* Python 2-monolit
* många fulla läsningar av FASTQ/intermediärer
* mycket temporär disk-I/O
* många externa kommandon utan tydlig typed state
* svårt att förstå exakt vilken evidens som ledde till en call
* svårt att targeta dyra steg mot kliniskt relevanta gener

`rinfuse-fc` ska bli en bro mellan nuvarande FusionCatcher-användning och en framtida native fusion-aware RNA-seq-aligner i `stars2`/`rmap`.

## Strategisk princip

Vi ska inte porta FusionCatcher rad-för-rad.

Vi ska istället reproducera beteendet på evidensnivå:

```text
FASTQ/BAM
  -> read registry
  -> externa aligners
  -> parser till normaliserad FusionEvidence
  -> candidate aggregator
  -> filter/labels
  -> rapport
```

I början är målet compatibility och observability, inte ny metodik.

## Föreslagen repo-struktur

```text
rinfuse-fc/
  Cargo.toml
  README.md
  LICENSE
  docs/
    architecture.md
    fusioncatcher_compatibility.md
    evidence_model.md
    mvp_plan.md
    dux4_rescue_notes.md
    file_formats.md
    cli.md
  crates/
    rinfuse-core/
      Cargo.toml
      src/
        lib.rs
        read_id.rs
        evidence.rs
        candidate.rs
        annotation.rs
        scoring.rs
        filters.rs
    rinfuse-io/
      Cargo.toml
      src/
        lib.rs
        fastq.rs
        bam.rs
        sam.rs
        psl.rs
        fc_intermediates.rs
        jsonl.rs
        parquet.rs
    rinfuse-orchestrator/
      Cargo.toml
      src/
        lib.rs
        step.rs
        command.rs
        manifest.rs
        cache.rs
        workdir.rs
        timing.rs
    rinfuse-fc/
      Cargo.toml
      src/
        lib.rs
        config.rs
        pipeline.rs
        steps/
          mod.rs
          fastq_prep.rs
          bowtie.rs
          star.rs
          bowtie2.rs
          blat.rs
          inspect_fc.rs
          aggregate.rs
          report.rs
    rinfuse-cli/
      Cargo.toml
      src/
        main.rs
        args.rs
        commands/
          mod.rs
          run.rs
          inspect_fc.rs
          compare.rs
          extract_reads.rs
  tests/
    fixtures/
      fusioncatcher_minimal/
      reads_small_R1.fq
      reads_small_R2.fq
    inspect_fc_smoke.rs
    read_registry.rs
    evidence_roundtrip.rs
```

## Crate-ansvar

### `rinfuse-core`

Ren domänmodell. Inga externa kommandon. Så lite I/O som möjligt.

Innehåller:

* read identifiers
* gene/transcript identifiers
* evidence structs
* fusion candidates
* labels och filters
* scoring primitives

### `rinfuse-io`

All parsing och serialisering.

Innehåller:

* FASTQ/BAM/SAM parsing
* PSL parsing för BLAT-output
* FusionCatcher intermediate readers
* JSONL writer/reader
* optional Parquet/Arrow writer

### `rinfuse-orchestrator`

Generic pipeline engine.

Innehåller:

* step trait
* command runner
* manifest
* checksum/cache
* workdir layout
* timing/loggning
* restart

### `rinfuse-fc`

FusionCatcher-compatible pipeline logic.

Innehåller:

* config
* pipeline graph
* wrappers runt Bowtie/STAR/Bowtie2/BLAT
* aggregation från externa aligner outputs
* compatibility-layer för FusionCatcher-output

### `rinfuse-cli`

CLI-binär.

Kommandon:

```bash
rinfuse-fc run
rinfuse-fc inspect-fc
rinfuse-fc compare
rinfuse-fc extract-reads
```

## CLI-utkast

### Körning

```bash
rinfuse-fc run \
  --input R1.fastq.gz,R2.fastq.gz \
  --data /path/to/fusioncatcher-data \
  --out out/ \
  --threads 32 \
  --memory-mode evidence \
  --aligners bowtie,star,bowtie2 \
  --skip-blat
```

### Inspektera befintlig FusionCatcher-körning

```bash
rinfuse-fc inspect-fc \
  --fc-out fusioncatcher_out/ \
  --reads R1.fastq.gz,R2.fastq.gz \
  --out inspected/
```

Output:

```text
inspected/
  evidence.jsonl
  supporting_reads.fq.gz
  candidates.tsv
  manifest.json
```

### Jämför original och Rust-pipeline

```bash
rinfuse-fc compare \
  --fc fusioncatcher_out/ \
  --rs rinfuse_out/ \
  --out comparison.tsv
```

### Extrahera reads

```bash
rinfuse-fc extract-reads \
  --reads R1.fastq.gz,R2.fastq.gz \
  --read-ids read_ids.txt \
  --out supporting_reads.fq.gz
```

## Workdir-layout

```text
out/
  manifest.json
  timings.tsv
  logs/
    bowtie.stderr.log
    star.stderr.log
    bowtie2.stderr.log
    blat.stderr.log
  registry/
    read_index.bin
    selected_reads.bin
  align/
    bowtie/
    star/
    bowtie2/
    blat/
  evidence/
    discordant_pairs.jsonl
    split_reads.jsonl
    junction_hits.jsonl
    validation_hits.jsonl
    all_evidence.jsonl
  candidates/
    raw_candidates.tsv
    filtered_candidates.tsv
  report/
    final-list_candidate-fusion-genes.tsv
    evidence_summary.tsv
```

## Memory modes

### `minimal`

Spara bara read names, mate relation och offsets där det går.

Bra för:

* väldigt stora WTS-körningar
* begränsat RAM

### `evidence`

Spara full sekvens/qual för reads som är eller kan bli intressanta:

* unmapped
* soft-clipped
* discordant
* supporting read IDs från intermediate files
* mates till intressanta reads

Detta bör vara default.

### `full`

Spara alla reads i RAM eller mmap-backed store.

Bra för:

* mindre paneler
* utveckling
* mycket RAM
* maximal minimering av FASTQ-rescans

## Core data model: utkast

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ReadId(pub u64);

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ReadName(pub Box<str>);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MateSide {
    R1,
    R2,
    Unknown,
}

#[derive(Debug, Clone)]
pub struct ReadRecordLite {
    pub id: ReadId,
    pub name: ReadName,
    pub mate: MateSide,
    pub len: u16,
    pub seq: Option<Vec<u8>>,
    pub qual: Option<Vec<u8>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GeneId(pub Box<str>);

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TranscriptId(pub Box<str>);

#[derive(Debug, Clone)]
pub struct AlignmentSegment {
    pub read_id: ReadId,
    pub chrom: Box<str>,
    pub start_0based: u64,
    pub end_0based: u64,
    pub strand: Strand,
    pub mapq: Option<u8>,
    pub cigar: Option<Box<str>>,
    pub gene: Option<GeneId>,
    pub transcript: Option<TranscriptId>,
}

#[derive(Debug, Clone)]
pub enum FusionEvidence {
    DiscordantPair(DiscordantPairEvidence),
    SplitRead(SplitReadEvidence),
    JunctionHit(JunctionHitEvidence),
    ValidationHit(ValidationHitEvidence),
}

#[derive(Debug, Clone)]
pub struct DiscordantPairEvidence {
    pub read_id_1: ReadId,
    pub read_id_2: ReadId,
    pub left: AlignmentSegment,
    pub right: AlignmentSegment,
    pub source: EvidenceSource,
}

#[derive(Debug, Clone)]
pub struct SplitReadEvidence {
    pub read_id: ReadId,
    pub left: AlignmentSegment,
    pub right: AlignmentSegment,
    pub anchor_left: u16,
    pub anchor_right: u16,
    pub mismatches: Option<u16>,
    pub source: EvidenceSource,
}

#[derive(Debug, Clone)]
pub struct JunctionHitEvidence {
    pub read_id: ReadId,
    pub gene_5p: GeneId,
    pub gene_3p: GeneId,
    pub transcript_5p: Option<TranscriptId>,
    pub transcript_3p: Option<TranscriptId>,
    pub anchor_min: u16,
    pub junction_sequence: Option<Vec<u8>>,
    pub source: EvidenceSource,
}

#[derive(Debug, Clone)]
pub struct ValidationHitEvidence {
    pub read_id: ReadId,
    pub candidate_id: CandidateId,
    pub aligner: ExternalAligner,
    pub score: Option<i32>,
    pub mismatches: Option<u16>,
    pub psl_or_sam_ref: Option<Box<str>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceSource {
    Bowtie,
    Star,
    Bowtie2,
    Blat,
    FusionCatcherIntermediate,
    RustDerived,
}
```

## Candidate model: utkast

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CandidateId(pub u64);

#[derive(Debug, Clone)]
pub struct FusionCandidate {
    pub id: CandidateId,
    pub gene_5p: GeneId,
    pub gene_3p: GeneId,
    pub symbol_5p: Option<Box<str>>,
    pub symbol_3p: Option<Box<str>>,
    pub evidence: Vec<FusionEvidence>,
    pub labels: Vec<FusionLabel>,
    pub filter_decisions: Vec<FilterDecision>,
    pub score: Option<FusionScore>,
}

#[derive(Debug, Clone)]
pub enum FusionLabel {
    KnownCancerFusion,
    KnownHealthyFusion,
    ReadthroughLike,
    ParalogLike,
    PseudogeneLike,
    ShortDistance,
    RepetitiveRegion,
    Dux4Class,
    IghClass,
    LowComplexity,
    Ambiguous,
}

#[derive(Debug, Clone)]
pub struct FilterDecision {
    pub label: FusionLabel,
    pub action: FilterAction,
    pub reason: Box<str>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FilterAction {
    Keep,
    SoftLabel,
    HideByDefault,
    Drop,
}
```

Viktig princip: många FusionCatcher-liknande filter ska först vara `SoftLabel`, inte `Drop`. Det är särskilt viktigt för DUX4/IGH där repetitivitet/paralogi annars kan ge falskt negativa.

## Orchestrator trait: utkast

```rust
pub trait Step {
    fn name(&self) -> &'static str;
    fn inputs(&self, ctx: &RunContext) -> anyhow::Result<Vec<PathBuf>>;
    fn outputs(&self, ctx: &RunContext) -> anyhow::Result<Vec<PathBuf>>;
    fn run(&self, ctx: &RunContext) -> anyhow::Result<StepResult>;
}

pub struct StepResult {
    pub status: StepStatus,
    pub walltime_ms: u128,
    pub outputs: Vec<PathBuf>,
    pub metrics: serde_json::Value,
}

pub enum StepStatus {
    Completed,
    SkippedCached,
    Failed,
}
```

## Manifest-utkast

```json
{
  "tool": "rinfuse-fc",
  "version": "0.1.0",
  "mode": "fc-compatible",
  "sample": "SAMPLE001",
  "inputs": {
    "reads": ["R1.fastq.gz", "R2.fastq.gz"],
    "data_dir": "/path/to/fusioncatcher-data"
  },
  "aligners": {
    "bowtie": { "path": "bowtie", "version": "1.2.3" },
    "star": { "path": "STAR", "version": "2.7.2b" },
    "bowtie2": { "path": "bowtie2", "version": "2.3.5.1" },
    "blat": { "path": null, "enabled": false }
  },
  "steps": [
    {
      "name": "read_registry",
      "status": "completed",
      "walltime_ms": 12000,
      "outputs": ["registry/read_index.bin"]
    }
  ]
}
```

## MVP-plan

### MVP 0: Repository skeleton

Leverans:

* Cargo workspace
* crate-struktur
* CLI med tomma kommandon
* docs enligt ovan
* CI med `cargo fmt`, `cargo clippy`, `cargo test`

Acceptans:

```bash
cargo test
rinfuse-fc --help
rinfuse-fc inspect-fc --help
```

### MVP 1: Read registry + extract-reads

Leverans:

* parse paired FASTQ
* bygg `ReadRegistry`
* läs `read_ids.txt`
* skriv matchande reads till FASTQ

Acceptans:

```bash
rinfuse-fc extract-reads \
  --reads tests/fixtures/reads_small_R1.fq,tests/fixtures/reads_small_R2.fq \
  --read-ids tests/fixtures/read_ids.txt \
  --out /tmp/supporting.fq
```

Testa:

* hittar R1/R2 med `/1`, `/2`, mellanslagsformat och rena read names
* output deterministisk
* missade read IDs rapporteras

### MVP 2: inspect-fc

Leverans:

* läs ett befintligt FusionCatcher-output
* hitta final report
* hitta supporting read IDs om filer finns
* skapa `evidence.jsonl`
* extrahera supporting reads med registry

Acceptans:

```bash
rinfuse-fc inspect-fc \
  --fc-out tests/fixtures/fusioncatcher_minimal \
  --reads tests/fixtures/reads_small_R1.fq,tests/fixtures/reads_small_R2.fq \
  --out /tmp/inspect
```

Output:

```text
/tmp/inspect/evidence.jsonl
/tmp/inspect/candidates.tsv
/tmp/inspect/supporting_reads.fq
/tmp/inspect/manifest.json
```

### MVP 3: External command runner

Leverans:

* robust command execution
* stdout/stderr logs
* exit code handling
* versionsinsamling
* manifest update

Acceptans:

* kör `echo`/mock-command i test
* fångar stderr
* misslyckad command ger tydligt fel
* `--dry-run` visar kommandon utan att köra

### MVP 4: STAR wrapper + parser

Leverans:

* kör STAR med befintligt index
* samla chimeric/split outputs
* parse till `SplitReadEvidence`

Acceptans:

* wrapper bygger kommando deterministiskt
* parser-test på liten STAR-output-fixture
* evidence JSONL skrivs

### MVP 5: Bowtie/Bowtie2/BLAT wrappers

Leverans:

* wrappers och parsers för SAM/PSL
* normalisering till evidence records

Acceptans:

* parsers testade på fixtures
* externa steg kan mockas i CI

### MVP 6: Candidate aggregation

Leverans:

* gruppera evidence per gene-pair
* summera spanning pairs, split reads, anchors, sources
* skriv `raw_candidates.tsv`

Acceptans:

* deterministisk sortering
* enkla fixtures ger förväntade candidates

### MVP 7: Compatibility comparison

Leverans:

```bash
rinfuse-fc compare --fc old/ --rs new/ --out compare.tsv
```

Jämför:

* gene-pair presence
* support counts
* split-read counts
* longest anchor
* aligner sources

Acceptans:

* tydlig diff för saknade/nya calls
* DUX4-fall kan granskas snabbt

## Initiala filer att skapa

### `README.md`

```markdown
# rinfuse-fc

Rust-based FusionCatcher-compatible orchestrator and fusion evidence engine.

First goal: reproduce FusionCatcher-style evidence using the same external aligners, while reducing redundant FASTQ/intermediate file reads and making the pipeline observable, restartable and testable.

This is not a line-by-line port of FusionCatcher.
```

### `docs/architecture.md`

Innehåll:

* motivation
* pipeline diagram
* crate responsibilities
* memory modes
* step graph
* workdir layout

### `docs/fusioncatcher_compatibility.md`

Innehåll:

* vilka delar vi försöker reproducera
* vilka delar som först bara läses från original-output
* vilka externa aligners används
* vilka outputfält ska jämföras
* known differences

### `docs/evidence_model.md`

Innehåll:

* `FusionEvidence`
* `FusionCandidate`
* source types
* filter labels
* scoring placeholders
* JSONL examples

### `docs/mvp_plan.md`

Innehåll:

* MVP 0–7 enligt ovan
* acceptanskriterier
* testdata
* risker

### `docs/dux4_rescue_notes.md`

Innehåll:

* varför DUX4 är svårt
* vilka samples i benchmarkset
* vilka FusionCatcher intermediates först visar DUX4-signalen
* vilka filter får inte vara hard drops
* targeted-mode-hypoteser

## Rekommenderad första implementation

Börja inte med att köra aligners.

Börja med:

```text
MVP 1: Read registry + extract-reads
MVP 2: inspect-fc
```

Motivering:

* ger direkt nytta vid analys av befintliga FusionCatcher-körningar
* minskar osäkerheten kring read ID-normalisering
* bygger grund för att undvika onödiga FASTQ-rescans
* gör det möjligt att reverse-engineera DUX4-fall systematiskt

## Risker

### Read name-normalisering

Olika verktyg representerar read names olika:

```text
READ/1
READ/2
READ 1:N:0:INDEX
READ 2:N:0:INDEX
READ
```

Detta måste lösas tidigt.

### Gzip FASTQ random access

Vanlig gzip är inte bra för random access. Första versionen bör streama och cachea intressanta reads. Senare kan vi stödja bgzip/index.

### FusionCatcher intermediate variation

Filnamn kan skilja mellan versioner/körlägen. `inspect-fc` bör vara tolerant och rapportera vilka filer som hittades/saknades.

### Filterportning

Filter ska först vara labels. Hard drop först när vi har benchmarkat mot kliniskt viktiga fall.

### BLAT-licens

BLAT ska vara optional. Manifestet ska tydligt visa om BLAT använts.

## Definition of done för fas 1

Fas 1 är klar när vi kan:

1. Läsa befintliga FusionCatcher outputs.
2. Extrahera supporting reads utan extra ad hoc-skript.
3. Skriva normaliserad evidence JSONL.
4. Jämföra FusionCatcher original-output mot rinfuse-fc-output.
5. Köra minst ett DUX4/FusionCatcher-only sample genom inspect/compare-flödet.
6. Visa vilka reads/evidenstyper som stödjer callen.

## Nästa steg efter fas 1

När inspect/registry fungerar:

1. lägg till command runner
2. wrapper för STAR
3. wrapper för Bowtie
4. parser till evidence model
5. candidate aggregation
6. compatibility compare
7. targeted DUX4 rescue mode

---

# Startpatch: konkret filutkast

Detta är ett praktiskt första patchförslag. Det ska ge ett kompilerande Rust workspace med:

* CLI-binär `rinfuse-fc`
* `extract-reads`-kommando
* read-name-normalisering
* enkel FASTQ-parser
* read registry som kan extrahera supporting reads
* plats för `inspect-fc` och `compare`

Syftet är att få ett första fungerande verktyg som löser ett verkligt problem: hämta ut read IDs från befintliga FusionCatcher-intermediärer utan ad hoc-skript och utan flera onödiga omläsningar.

## Filer att skapa i första commit

```text
rinfuse-fc/
  Cargo.toml
  README.md
  crates/
    rinfuse-core/
      Cargo.toml
      src/
        lib.rs
        read_id.rs
    rinfuse-io/
      Cargo.toml
      src/
        lib.rs
        fastq.rs
        read_registry.rs
    rinfuse-cli/
      Cargo.toml
      src/
        main.rs
        args.rs
        commands/
          mod.rs
          extract_reads.rs
          inspect_fc.rs
          compare.rs
  tests/
    fixtures/
      reads_R1.fq
      reads_R2.fq
      read_ids.txt
```

## Root `Cargo.toml`

```toml
[workspace]
members = [
    "crates/rinfuse-core",
    "crates/rinfuse-io",
    "crates/rinfuse-cli",
]
resolver = "2"

[workspace.package]
edition = "2021"
license = "GPL-3.0-or-later"
rust-version = "1.76"

[workspace.dependencies]
anyhow = "1"
clap = { version = "4", features = ["derive"] }
flate2 = "1"
serde = { version = "1", features = ["derive"] }
serde_json = "1"
thiserror = "1"
tracing = "0.1"
tracing-subscriber = { version = "0.3", features = ["env-filter"] }
```

## `README.md`

````markdown
# rinfuse-fc

Rust-based FusionCatcher-compatible orchestrator and fusion evidence engine.

The first goal is not to replace FusionCatcher's aligner strategy. Instead, this project starts by making FusionCatcher-like workflows more observable, testable and efficient:

- read registry
- supporting-read extraction
- FusionCatcher-output inspection
- normalized fusion evidence
- external aligner orchestration in later phases

## First milestone

```bash
rinfuse-fc extract-reads \
  --reads R1.fastq.gz,R2.fastq.gz \
  --read-ids read_ids.txt \
  --out supporting_reads.fq.gz
````

This extracts requested read IDs and their mates from paired FASTQ files.

## Status

Experimental. Intended initially for internal validation against FusionCatcher outputs.

````

## `crates/rinfuse-core/Cargo.toml`

```toml
[package]
name = "rinfuse-core"
version = "0.1.0"
edition.workspace = true
license.workspace = true
rust-version.workspace = true

[dependencies]
serde.workspace = true
thiserror.workspace = true
````

## `crates/rinfuse-core/src/lib.rs`

```rust
pub mod read_id;

pub use read_id::{normalize_read_name, MateSide, NormalizedReadName, ReadId};
```

## `crates/rinfuse-core/src/read_id.rs`

```rust
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ReadId(pub u64);

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum MateSide {
    R1,
    R2,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NormalizedReadName {
    pub base: String,
    pub mate: MateSide,
}

/// Normalize common Illumina/FusionCatcher/read-aligner read-name variants.
///
/// Handles examples like:
/// - `READ/1`
/// - `READ/2`
/// - `READ 1:N:0:INDEX`
/// - `READ 2:N:0:INDEX`
/// - `@READ/1`
/// - `READ_supports_fusion_junction/1` should first be stripped by caller if needed
pub fn normalize_read_name(raw: &str) -> NormalizedReadName {
    let mut s = raw.trim();

    if let Some(rest) = s.strip_prefix('@') {
        s = rest;
    }

    // FASTQ header may contain metadata after first whitespace.
    let mut parts = s.split_whitespace();
    let first = parts.next().unwrap_or("");
    let second = parts.next();

    if let Some(meta) = second {
        if meta.starts_with("1:") {
            return NormalizedReadName {
                base: first.to_string(),
                mate: MateSide::R1,
            };
        }
        if meta.starts_with("2:") {
            return NormalizedReadName {
                base: first.to_string(),
                mate: MateSide::R2,
            };
        }
    }

    if let Some(base) = first.strip_suffix("/1") {
        return NormalizedReadName {
            base: base.to_string(),
            mate: MateSide::R1,
        };
    }

    if let Some(base) = first.strip_suffix("/2") {
        return NormalizedReadName {
            base: base.to_string(),
            mate: MateSide::R2,
        };
    }

    NormalizedReadName {
        base: first.to_string(),
        mate: MateSide::Unknown,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_slash_mates() {
        let r1 = normalize_read_name("@READ123/1");
        assert_eq!(r1.base, "READ123");
        assert_eq!(r1.mate, MateSide::R1);

        let r2 = normalize_read_name("READ123/2");
        assert_eq!(r2.base, "READ123");
        assert_eq!(r2.mate, MateSide::R2);
    }

    #[test]
    fn normalizes_illumina_metadata_mates() {
        let r1 = normalize_read_name("@READ123 1:N:0:ACGT");
        assert_eq!(r1.base, "READ123");
        assert_eq!(r1.mate, MateSide::R1);

        let r2 = normalize_read_name("@READ123 2:N:0:ACGT");
        assert_eq!(r2.base, "READ123");
        assert_eq!(r2.mate, MateSide::R2);
    }

    #[test]
    fn keeps_unknown_when_no_mate_marker() {
        let r = normalize_read_name("READ123");
        assert_eq!(r.base, "READ123");
        assert_eq!(r.mate, MateSide::Unknown);
    }
}
```

## `crates/rinfuse-io/Cargo.toml`

```toml
[package]
name = "rinfuse-io"
version = "0.1.0"
edition.workspace = true
license.workspace = true
rust-version.workspace = true

[dependencies]
anyhow.workspace = true
flate2.workspace = true
rinfuse-core = { path = "../rinfuse-core" }
serde.workspace = true
thiserror.workspace = true
```

## `crates/rinfuse-io/src/lib.rs`

```rust
pub mod fastq;
pub mod read_registry;

pub use fastq::{open_maybe_gz, FastqRecord, FastqReader};
pub use read_registry::{ReadRegistry, ReadRegistryBuildOptions};
```

## `crates/rinfuse-io/src/fastq.rs`

```rust
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub header: String,
    pub seq: Vec<u8>,
    pub plus: String,
    pub qual: Vec<u8>,
}

pub fn open_maybe_gz(path: &Path) -> Result<Box<dyn Read>> {
    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    if path.extension().is_some_and(|ext| ext == "gz") {
        Ok(Box::new(MultiGzDecoder::new(file)))
    } else {
        Ok(Box::new(file))
    }
}

pub struct FastqReader<R: BufRead> {
    reader: R,
    line_no: usize,
}

impl FastqReader<BufReader<Box<dyn Read>>> {
    pub fn from_path(path: &Path) -> Result<Self> {
        let reader = BufReader::new(open_maybe_gz(path)?);
        Ok(Self::new(reader))
    }
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader, line_no: 0 }
    }

    pub fn next_record(&mut self) -> Result<Option<FastqRecord>> {
        let mut header = String::new();
        if self.reader.read_line(&mut header)? == 0 {
            return Ok(None);
        }
        self.line_no += 1;

        let mut seq = String::new();
        let mut plus = String::new();
        let mut qual = String::new();

        if self.reader.read_line(&mut seq)? == 0
            || self.reader.read_line(&mut plus)? == 0
            || self.reader.read_line(&mut qual)? == 0
        {
            bail!("truncated FASTQ record starting at line {}", self.line_no);
        }
        self.line_no += 3;

        if !header.starts_with('@') {
            bail!("expected FASTQ header at line {}, got {:?}", self.line_no - 3, header.trim_end());
        }
        if !plus.starts_with('+') {
            bail!("expected FASTQ plus line at line {}, got {:?}", self.line_no - 1, plus.trim_end());
        }

        let seq = seq.trim_end_matches(['
', '
']).as_bytes().to_vec();
        let qual = qual.trim_end_matches(['
', '
']).as_bytes().to_vec();

        if seq.len() != qual.len() {
            bail!(
                "FASTQ seq/qual length mismatch at record starting line {}: {} != {}",
                self.line_no - 3,
                seq.len(),
                qual.len()
            );
        }

        Ok(Some(FastqRecord {
            header: header.trim_end_matches(['
', '
']).to_string(),
            seq,
            plus: plus.trim_end_matches(['
', '
']).to_string(),
            qual,
        }))
    }
}

impl FastqRecord {
    pub fn write_to<W: Write>(&self, mut writer: W) -> io::Result<()> {
        writer.write_all(self.header.as_bytes())?;
        writer.write_all(b"
")?;
        writer.write_all(&self.seq)?;
        writer.write_all(b"
")?;
        writer.write_all(self.plus.as_bytes())?;
        writer.write_all(b"
")?;
        writer.write_all(&self.qual)?;
        writer.write_all(b"
")?;
        Ok(())
    }
}
```

## `crates/rinfuse-io/src/read_registry.rs`

```rust
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
```

Obs: första versionen skriver output som okomprimerad FASTQ även om filnamnet slutar på `.gz`. Lägg gzip-writer i nästa patch om ni vill. Första målet är correctness.

## `crates/rinfuse-cli/Cargo.toml`

```toml
[package]
name = "rinfuse-cli"
version = "0.1.0"
edition.workspace = true
license.workspace = true
rust-version.workspace = true

[[bin]]
name = "rinfuse-fc"
path = "src/main.rs"

[dependencies]
anyhow.workspace = true
clap.workspace = true
rinfuse-core = { path = "../rinfuse-core" }
rinfuse-io = { path = "../rinfuse-io" }
tracing.workspace = true
tracing-subscriber.workspace = true
```

## `crates/rinfuse-cli/src/main.rs`

```rust
mod args;
mod commands;

use anyhow::Result;
use args::{Cli, Commands};
use clap::Parser;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    let cli = Cli::parse();

    match cli.command {
        Commands::ExtractReads(args) => commands::extract_reads::run(args),
        Commands::InspectFc(args) => commands::inspect_fc::run(args),
        Commands::Compare(args) => commands::compare::run(args),
    }
}
```

## `crates/rinfuse-cli/src/args.rs`

```rust
use clap::{Args, Parser, Subcommand};
use std::path::PathBuf;

#[derive(Debug, Parser)]
#[command(name = "rinfuse-fc")]
#[command(about = "Rust FusionCatcher-compatible orchestrator and evidence tools")]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Extract read records by read ID from FASTQ files.
    ExtractReads(ExtractReadsArgs),

    /// Inspect an existing FusionCatcher output directory.
    InspectFc(InspectFcArgs),

    /// Compare FusionCatcher and rinfuse-fc outputs.
    Compare(CompareArgs),
}

#[derive(Debug, Args)]
pub struct ExtractReadsArgs {
    /// Comma-separated FASTQ paths, for example R1.fq.gz,R2.fq.gz.
    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,

    /// Text file with one read ID per line.
    #[arg(long)]
    pub read_ids: PathBuf,

    /// Output FASTQ path. First MVP writes plain FASTQ.
    #[arg(long)]
    pub out: PathBuf,

    /// Output missing read IDs here.
    #[arg(long)]
    pub missing_out: Option<PathBuf>,
}

#[derive(Debug, Args)]
pub struct InspectFcArgs {
    #[arg(long)]
    pub fc_out: PathBuf,

    #[arg(long, value_delimiter = ',')]
    pub reads: Vec<PathBuf>,

    #[arg(long)]
    pub out: PathBuf,
}

#[derive(Debug, Args)]
pub struct CompareArgs {
    #[arg(long)]
    pub fc: PathBuf,

    #[arg(long)]
    pub rs: PathBuf,

    #[arg(long)]
    pub out: PathBuf,
}
```

## `crates/rinfuse-cli/src/commands/mod.rs`

```rust
pub mod compare;
pub mod extract_reads;
pub mod inspect_fc;
```

## `crates/rinfuse-cli/src/commands/extract_reads.rs`

```rust
use crate::args::ExtractReadsArgs;
use anyhow::Result;
use rinfuse_io::{ReadRegistry, ReadRegistryBuildOptions};

pub fn run(args: ExtractReadsArgs) -> Result<()> {
    let mut registry = ReadRegistry::from_read_id_file(&args.read_ids)?;
    let opts = ReadRegistryBuildOptions::default();

    registry.collect_from_fastq_paths(&args.reads, &opts)?;
    registry.write_fastq(&args.out)?;

    if let Some(missing_out) = args.missing_out.as_deref() {
        registry.write_missing(missing_out)?;
    }

    eprintln!(
        "requested_bases={} found_records={} missing_bases={}",
        registry.requested_base_count(),
        registry.found_record_count(),
        registry.missing_base_count()
    );

    Ok(())
}
```

## `crates/rinfuse-cli/src/commands/inspect_fc.rs`

```rust
use crate::args::InspectFcArgs;
use anyhow::{bail, Result};

pub fn run(_args: InspectFcArgs) -> Result<()> {
    bail!("inspect-fc is not implemented yet. Start with extract-reads and read-name normalization.")
}
```

## `crates/rinfuse-cli/src/commands/compare.rs`

```rust
use crate::args::CompareArgs;
use anyhow::{bail, Result};

pub fn run(_args: CompareArgs) -> Result<()> {
    bail!("compare is not implemented yet. It will compare FusionCatcher and rinfuse-fc candidate outputs.")
}
```

## Testfixtures

### `tests/fixtures/reads_R1.fq`

```text
@READ_A/1
ACGTACGT
+
FFFFFFFF
@READ_B/1
TTTTCCCC
+
FFFFFFFF
@READ_C/1
GGGGAAAA
+
FFFFFFFF
```

### `tests/fixtures/reads_R2.fq`

```text
@READ_A/2
TGCATGCA
+
FFFFFFFF
@READ_B/2
AAAAGGGG
+
FFFFFFFF
@READ_C/2
CCCCTTTT
+
FFFFFFFF
```

### `tests/fixtures/read_ids.txt`

```text
READ_A/1
READ_C/2
READ_MISSING/1
```

Första körning:

```bash
cargo run -p rinfuse-cli --bin rinfuse-fc -- extract-reads \
  --reads tests/fixtures/reads_R1.fq,tests/fixtures/reads_R2.fq \
  --read-ids tests/fixtures/read_ids.txt \
  --out /tmp/supporting.fq \
  --missing-out /tmp/missing.txt
```

Förväntat:

```text
requested_bases=3 found_records=4 missing_bases=1
```

Varför 4 records? Eftersom `READ_A` och `READ_C` matchar basnamn och båda mates finns i input.



## Andra patch efter detta

När startpatchen fungerar, nästa patch bör vara:

```text
inspect-fc MVP:
  - detect common FusionCatcher final report files
  - detect supporting read ID files
  - collect read IDs
  - call ReadRegistry extraction
  - write candidates.tsv + manifest.json
```
