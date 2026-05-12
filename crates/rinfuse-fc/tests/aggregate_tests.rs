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
        let j1 = r#"{"seg1":{"chrom":"chr1","pos_1based":1501,"strand":"+","segment_start_1based":10,"cigar":"10M"},"seg2":{"chrom":"chr2","pos_1based":5501,"strand":"-","segment_start_1based":20,"cigar":"20M"},"junction_type":1,"repeat_left":0,"repeat_right":0,"read_name":"READ1","num_chimeric_reads":null,"max_overhang":10,"source":"Star","source_format":"StarChimericV14","raw_fields":[]}"#;
        // Same gene pair, different read
        let j2 = r#"{"seg1":{"chrom":"chr1","pos_1based":1601,"strand":"+","segment_start_1based":10,"cigar":"10M"},"seg2":{"chrom":"chr2","pos_1based":5601,"strand":"-","segment_start_1based":20,"cigar":"20M"},"junction_type":2,"repeat_left":0,"repeat_right":0,"read_name":"READ2","num_chimeric_reads":null,"max_overhang":20,"source":"Star","source_format":"StarChimericV14","raw_fields":[]}"#;
        // Different gene pair, one intergenic (unknown)
        let j3 = r#"{"seg1":{"chrom":"chr1","pos_1based":1501,"strand":"+","segment_start_1based":10,"cigar":"10M"},"seg2":{"chrom":"chrX","pos_1based":10000,"strand":"-","segment_start_1based":20,"cigar":"20M"},"junction_type":1,"repeat_left":0,"repeat_right":0,"read_name":"READ3","num_chimeric_reads":null,"max_overhang":5,"source":"Star","source_format":"StarChimericV14","raw_fields":[]}"#;

        fs::write(&junctions_path, format!("{}\n{}\n{}\n", j1, j2, j3)).unwrap();

        // 3. Run aggregation
        let mut candidates = aggregate_star_junctions(&junctions_path, &genes_path).unwrap();
        // Sort candidates by gene_a to make test deterministic
        candidates.sort_by(|a, b| a.gene_5p.cmp(&b.gene_5p));

        assert_eq!(candidates.len(), 2);

        // Candidate 1: GENE_A -- GENE_B
        assert_eq!(candidates[0].gene_5p, "GENE_A");
        assert_eq!(candidates[0].gene_3p, "GENE_B");
        assert_eq!(candidates[0].support_junction_count, 2);
        assert_eq!(candidates[0].unique_read_count, 2);
        assert_eq!(candidates[0].max_overhang, 20);
        let mut jtypes = candidates[0].junction_types.clone();
        jtypes.sort();
        assert_eq!(jtypes, vec![1, 2]);

        // Candidate 2: GENE_A -- UNKNOWN_chrX
        assert_eq!(candidates[1].gene_5p, "GENE_A");
        assert_eq!(candidates[1].gene_3p, "UNKNOWN_chrX");
        assert_eq!(candidates[1].support_junction_count, 1);
        assert_eq!(candidates[1].unique_read_count, 1);
        assert_eq!(candidates[1].max_overhang, 5);
        assert_eq!(candidates[1].junction_types, vec![1]);
    }

    #[test]
    fn aggregate_star_keeps_opposite_orientations_distinct() {
        let temp = tempdir().unwrap();

        let genes_path = temp.path().join("genes.tsv");
        let genes_content = "chrom\tstart_0based\tend_0based\tstrand\tgene_id\tgene_symbol\n\
chr9\t1000\t2000\t+\tENSG_BCR\tBCR\n\
chr22\t5000\t6000\t-\tENSG_ABL1\tABL1\n";
        fs::write(&genes_path, genes_content).unwrap();

        let junctions_path = temp.path().join("junctions.jsonl");
        let bcr_to_abl1 = r#"{"seg1":{"chrom":"chr9","pos_1based":1501,"strand":"+","segment_start_1based":10,"cigar":"10M"},"seg2":{"chrom":"chr22","pos_1based":5501,"strand":"-","segment_start_1based":20,"cigar":"20M"},"junction_type":1,"repeat_left":0,"repeat_right":0,"read_name":"READ_FWD","num_chimeric_reads":null,"max_overhang":25,"source":"Star","source_format":"StarChimericV14","raw_fields":[]}"#;
        let abl1_to_bcr = r#"{"seg1":{"chrom":"chr22","pos_1based":5501,"strand":"-","segment_start_1based":10,"cigar":"10M"},"seg2":{"chrom":"chr9","pos_1based":1501,"strand":"+","segment_start_1based":20,"cigar":"20M"},"junction_type":1,"repeat_left":0,"repeat_right":0,"read_name":"READ_REV","num_chimeric_reads":null,"max_overhang":20,"source":"Star","source_format":"StarChimericV14","raw_fields":[]}"#;
        fs::write(
            &junctions_path,
            format!("{}\n{}\n", bcr_to_abl1, abl1_to_bcr),
        )
        .unwrap();

        let candidates = aggregate_star_junctions(&junctions_path, &genes_path).unwrap();

        assert_eq!(candidates.len(), 2);
        assert!(candidates
            .iter()
            .any(|c| c.gene_5p == "BCR" && c.gene_3p == "ABL1"));
        assert!(candidates
            .iter()
            .any(|c| c.gene_5p == "ABL1" && c.gene_3p == "BCR"));
    }

    #[test]
    fn aggregate_uses_zero_based_lookup_positions() {
        let temp = tempdir().unwrap();

        let genes_path = temp.path().join("genes.tsv");
        let genes_content = "chrom\tstart_0based\tend_0based\tstrand\tgene_id\tgene_symbol\n\
chr1\t100\t200\t+\tENSG000001\tGENE_A\n\
chr2\t500\t600\t-\tENSG000002\tGENE_B\n";
        fs::write(&genes_path, genes_content).unwrap();

        let hit_path = temp.path().join("hit.jsonl");
        let hit = r#"{"seg1":{"chrom":"chr1","pos_1based":101,"strand":"+","segment_start_1based":10,"cigar":"10M"},"seg2":{"chrom":"chr2","pos_1based":501,"strand":"-","segment_start_1based":20,"cigar":"20M"},"junction_type":1,"repeat_left":0,"repeat_right":0,"read_name":"READ_HIT","num_chimeric_reads":null,"max_overhang":null,"source":"Star","source_format":"StarChimericV14","raw_fields":[]}"#;
        fs::write(&hit_path, format!("{}\n", hit)).unwrap();

        let miss_path = temp.path().join("miss.jsonl");
        let miss = r#"{"seg1":{"chrom":"chr1","pos_1based":100,"strand":"+","segment_start_1based":10,"cigar":"10M"},"seg2":{"chrom":"chr2","pos_1based":500,"strand":"-","segment_start_1based":20,"cigar":"20M"},"junction_type":1,"repeat_left":0,"repeat_right":0,"read_name":"READ_MISS","num_chimeric_reads":null,"max_overhang":null,"source":"Star","source_format":"StarChimericV14","raw_fields":[]}"#;
        fs::write(&miss_path, format!("{}\n", miss)).unwrap();

        let hit_candidates = aggregate_star_junctions(&hit_path, &genes_path).unwrap();
        assert_eq!(hit_candidates.len(), 1);
        assert_eq!(hit_candidates[0].gene_5p, "GENE_A");
        assert_eq!(hit_candidates[0].gene_3p, "GENE_B");

        let miss_candidates = aggregate_star_junctions(&miss_path, &genes_path).unwrap();
        assert_eq!(miss_candidates.len(), 1);
        assert_eq!(miss_candidates[0].gene_5p, "UNKNOWN_chr1");
        assert_eq!(miss_candidates[0].gene_3p, "UNKNOWN_chr2");
    }
}
