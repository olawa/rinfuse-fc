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

        let seq = seq.trim_end_matches(['\n', '\r']).as_bytes().to_vec();
        let qual = qual.trim_end_matches(['\n', '\r']).as_bytes().to_vec();

        if seq.len() != qual.len() {
            bail!(
                "FASTQ seq/qual length mismatch at record starting line {}: {} != {}",
                self.line_no - 3,
                seq.len(),
                qual.len()
            );
        }

        Ok(Some(FastqRecord {
            header: header.trim_end_matches(['\n', '\r']).to_string(),
            seq,
            plus: plus.trim_end_matches(['\n', '\r']).to_string(),
            qual,
        }))
    }
}

impl FastqRecord {
    pub fn write_to<W: Write>(&self, mut writer: W) -> io::Result<()> {
        writer.write_all(self.header.as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&self.seq)?;
        writer.write_all(b"\n")?;
        writer.write_all(self.plus.as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(&self.qual)?;
        writer.write_all(b"\n")?;
        Ok(())
    }
}
