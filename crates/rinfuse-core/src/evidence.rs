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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StarChimericSourceFormat {
    StarChimericV14,
    StarChimericExtended,
    UnknownExtraColumns,
}

/// One end (segment) of a chimeric alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChimericSegment {
    pub chrom: String,
    pub pos_1based: u64,
    pub strand: Strand,
    /// Alignment start within the read.
    pub segment_start_1based: u64,
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
    /// Number of chimeric reads spanning this junction when known.
    pub num_chimeric_reads: Option<u32>,
    /// Maximum overhang of supporting reads when known.
    pub max_overhang: Option<u32>,
    /// Source tool.
    pub source: EvidenceSource,
    /// STAR row layout detected by the parser.
    pub source_format: StarChimericSourceFormat,
    /// Raw TSV fields, preserved for debugging.
    pub raw_fields: Vec<String>,
}

impl StarChimericJunction {
    pub fn pos1_1based(&self) -> u64 {
        self.seg1.pos_1based
    }

    pub fn pos2_1based(&self) -> u64 {
        self.seg2.pos_1based
    }

    pub fn segment_start1_1based(&self) -> u64 {
        self.seg1.segment_start_1based
    }

    pub fn segment_start2_1based(&self) -> u64 {
        self.seg2.segment_start_1based
    }

    pub fn pos1_lookup_0based(&self) -> Option<u64> {
        self.pos1_1based().checked_sub(1)
    }

    pub fn pos2_lookup_0based(&self) -> Option<u64> {
        self.pos2_1based().checked_sub(1)
    }

    /// TSV header matching tsv_row() output.
    pub fn tsv_header() -> &'static str {
        "chrom1\tpos1_1based\tstrand1\tsegment_start1_1based\tcigar1\tchrom2\tpos2_1based\tstrand2\tsegment_start2_1based\tcigar2\tjunction_type\trepeat_left\trepeat_right\tread_name\tsource_format\tnum_chimeric_reads\tmax_overhang"
    }

    /// Produce a TSV row for this record.
    pub fn tsv_row(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.seg1.chrom,
            self.seg1.pos_1based,
            self.seg1.strand,
            self.seg1.segment_start_1based,
            self.seg1.cigar,
            self.seg2.chrom,
            self.seg2.pos_1based,
            self.seg2.strand,
            self.seg2.segment_start_1based,
            self.seg2.cigar,
            self.junction_type,
            self.repeat_left,
            self.repeat_right,
            self.read_name,
            match self.source_format {
                StarChimericSourceFormat::StarChimericV14 => "StarChimericV14",
                StarChimericSourceFormat::StarChimericExtended => "StarChimericExtended",
                StarChimericSourceFormat::UnknownExtraColumns => "UnknownExtraColumns",
            },
            self.num_chimeric_reads
                .map(|v| v.to_string())
                .unwrap_or_default(),
            self.max_overhang.map(|v| v.to_string()).unwrap_or_default(),
        )
    }
}
