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
