#[cfg(test)]
mod aggregate_tests {
    use rinfuse_fc::steps::aggregate::aggregate_star_junctions;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn test_aggregate_star_junctions() {
        let temp = tempdir().unwrap();

        // 1. Create fake gene intervals TSV
        let genes_path = temp.path().join("genes.tsv");
        let genes_content = "chrom\tstart_0based\tend_0based\tstrand\tgene_id\tgene_symbol\n\
chr1\t1000\t2000\t+\tENSG000001\tGENE_A\n\
chr2\t5000\t6000\t-\tENSG000002\tGENE_B\n\
chr3\t100\t500\t+\tENSG000003\tGENE_C\n";
        fs::write(&genes_path, genes_content).unwrap();

        // 2. Create fake STAR junctions JSONL
        let junctions_path = temp.path().join("junctions.jsonl");
        let j1 = r#"{"seg1":{"chrom":"chr1","genomic_pos":1500,"strand":"+","read_start":10,"cigar":"10M"},"seg2":{"chrom":"chr2","genomic_pos":5500,"strand":"-","read_start":20,"cigar":"20M"},"junction_type":1,"repeat_left":0,"repeat_right":0,"read_name":"READ1","num_chimeric_reads":1,"max_overhang":10,"source":"Star","raw_fields":[]}"#;
        // Same gene pair, different read
        let j2 = r#"{"seg1":{"chrom":"chr1","genomic_pos":1600,"strand":"+","read_start":10,"cigar":"10M"},"seg2":{"chrom":"chr2","genomic_pos":5600,"strand":"-","read_start":20,"cigar":"20M"},"junction_type":2,"repeat_left":0,"repeat_right":0,"read_name":"READ2","num_chimeric_reads":1,"max_overhang":20,"source":"Star","raw_fields":[]}"#;
        // Different gene pair, one intergenic (unknown)
        let j3 = r#"{"seg1":{"chrom":"chr1","genomic_pos":1500,"strand":"+","read_start":10,"cigar":"10M"},"seg2":{"chrom":"chrX","genomic_pos":9999,"strand":"-","read_start":20,"cigar":"20M"},"junction_type":1,"repeat_left":0,"repeat_right":0,"read_name":"READ3","num_chimeric_reads":1,"max_overhang":5,"source":"Star","raw_fields":[]}"#;

        fs::write(&junctions_path, format!("{}\n{}\n{}\n", j1, j2, j3)).unwrap();

        // 3. Run aggregation
        let mut candidates = aggregate_star_junctions(&junctions_path, &genes_path).unwrap();
        // Sort candidates by gene_a to make test deterministic
        candidates.sort_by(|a, b| a.gene_a.cmp(&b.gene_a));

        assert_eq!(candidates.len(), 2);

        // Candidate 1: GENE_A -- GENE_B
        assert_eq!(candidates[0].gene_a, "GENE_A");
        assert_eq!(candidates[0].gene_b, "GENE_B");
        assert_eq!(candidates[0].support_junction_count, 2);
        assert_eq!(candidates[0].unique_read_count, 2);
        assert_eq!(candidates[0].max_overhang, 20);
        let mut jtypes = candidates[0].junction_types.clone();
        jtypes.sort();
        assert_eq!(jtypes, vec![1, 2]);

        // Candidate 2: GENE_A -- UNKNOWN_chrX
        assert_eq!(candidates[1].gene_a, "GENE_A");
        assert_eq!(candidates[1].gene_b, "UNKNOWN_chrX");
        assert_eq!(candidates[1].support_junction_count, 1);
        assert_eq!(candidates[1].unique_read_count, 1);
        assert_eq!(candidates[1].max_overhang, 5);
        assert_eq!(candidates[1].junction_types, vec![1]);
    }
}
