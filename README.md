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
```

This extracts requested read IDs and their mates from paired FASTQ files.

## Status

Experimental. Intended initially for internal validation against FusionCatcher outputs.
