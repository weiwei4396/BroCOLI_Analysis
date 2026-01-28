
library(bambu)

fasta_file_path <- "/data/workdir/panw/Data/SIRV_Set2/SIRV.fasta"

gtf_file_path <- "/data/workdir/panw/Data/SIRV_Set2/SIRV_100.gtf"

annotations = prepareAnnotations(gtf_file_path)

bam_file_path <- "/data/workdir/panw/Data/SIRV_Set2/mini230_sorted.bam" 

se_bambu <- bambu(reads = bam_file_path, annotations = annotations, genome = fasta_file_path, NDR = 0.039)

writeBambuOutput(se_bambu, path = "/data/workdir/panw/Data/SIRV_Set2/Quantification/Bambu")
