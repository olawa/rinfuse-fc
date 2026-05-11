use serde::{Deserialize, Serialize};

/// Strand orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Strand {
    #[serde(rename = "+")]
    Plus,
    #[serde(rename = "-")]
    Minus,
    #[serde(rename = ".")]
    Unknown,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Self::Plus,
            '-' => Self::Minus,
            _ => Self::Unknown,
        }
    }
    pub fn as_char(self) -> char {
        match self {
            Self::Plus => '+',
            Self::Minus => '-',
            Self::Unknown => '.',
        }
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

/// Source tool that produced the chimeric evidence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum EvidenceSource {
    Star,
    Other(String),
}

/// One end (segment) of a chimeric alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChimericSegment {
    pub chrom: String,
    pub genomic_pos: u64,
    pub strand: Strand,
    /// Alignment start within the read.
    pub read_start: u64,
    /// CIGAR string for this segment.
    pub cigar: String,
}

/// A single row from STAR Chimeric.out.junction.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StarChimericJunction {
    /// Left donor segment.
    pub seg1: ChimericSegment,
    /// Right acceptor segment.
    pub seg2: ChimericSegment,
    /// STAR junction type: -1 (not canonical), 0 (not determined), 1 (GT/AG), 2 (CT/AC).
    pub junction_type: i32,
    /// Repeat length on the left side of the junction.
    pub repeat_left: u32,
    /// Repeat length on the right side of the junction.
    pub repeat_right: u32,
    /// Read name that supports this junction.
    pub read_name: String,
    /// Number of chimeric reads spanning this junction.
    pub num_chimeric_reads: u32,
    /// Maximum overhang of supporting reads.
    pub max_overhang: u32,
    /// Source tool.
    pub source: EvidenceSource,
    /// Raw TSV fields, preserved for debugging.
    pub raw_fields: Vec<String>,
}

impl StarChimericJunction {
    /// TSV header matching tsv_row() output.
    pub fn tsv_header() -> &'static str {
        "chrom1\tpos1\tstrand1\tstart1\tcigar1\tchrom2\tpos2\tstrand2\tstart2\tcigar2\t\
junction_type\trepeat_left\trepeat_right\tread_name\tnum_chimeric_reads\tmax_overhang"
    }

    /// Produce a TSV row for this record.
    pub fn tsv_row(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.seg1.chrom,
            self.seg1.genomic_pos,
            self.seg1.strand,
            self.seg1.read_start,
            self.seg1.cigar,
            self.seg2.chrom,
            self.seg2.genomic_pos,
            self.seg2.strand,
            self.seg2.read_start,
            self.seg2.cigar,
            self.junction_type,
            self.repeat_left,
            self.repeat_right,
            self.read_name,
            self.num_chimeric_reads,
            self.max_overhang,
        )
    }
}
