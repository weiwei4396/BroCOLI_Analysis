library(Isosceles)

bam_file <- "/data/workdir/panw/Data/SIRV_Set2/mini230_sorted.bam"
gtf_file <- "/data/workdir/panw/Data/SIRV_Set2/SIRV_100.gtf"
genome_fasta_file <- "/data/workdir/panw/Data/SIRV_Set2/SIRV.fasta"

# 需要跑bulk的流程;
bam_files <- c(Sample = bam_file)
bam_parsed <- bam_to_read_structures(
    bam_files = bam_files
)

transcript_data <- prepare_transcripts(
    gtf_file = gtf_file,
    genome_fasta_file = genome_fasta_file,
    bam_parsed = bam_parsed,
    min_bam_splice_read_count = 2,
    min_bam_splice_fraction = 0.01
)

# de_novo_strict/strict/de_novo_loose 
se_tcc <- bam_to_tcc(
    bam_files = bam_files,
    transcript_data = transcript_data,
    run_mode = "strict",
    min_read_count = 20,
    min_relative_expression = 0
)

se_transcript <- tcc_to_transcript(
    se_tcc = se_tcc,
    use_length_normalization = TRUE
)

se_gene <- tcc_to_gene(
    se_tcc = se_tcc
)


counts_data <- assays(se_transcript)$counts
counts_data <- as.matrix(counts_data)
counts_data <- as.data.frame(counts_data)
counts_data$raw_name <- rowData(se_transcript)$compatible_tx
 
write.table(counts_data, file = "/data/workdir/panw/Data/SIRV_Set2/Quantification/Isosceles/counts_transcript.txt",sep = "\t",row.names = TRUE, col.names = TRUE, quote = FALSE)



