library(seqinr)
data_path <- "~/dev/IBA_data_April2023/"

# Read in nochimera cluster data
nochimera_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_nochimera_clusters_SE_2019.tsv"),sep="\t",header=TRUE)

# Read in sequence data
cat("Reading asv sequences. This may take a while...\n")
seqs <- read.fasta(paste0(data_path,"data/Sweden/CO1_asv_seqs_SE_2019.fasta"))

# Extract representative asvs for the Lepidoptera_x clusters
asvs <- nochimera_clusters$ASV[nochimera_clusters$Family=="Lepidoptera_X" & nochimera_clusters$representative==1]

# Write the fasta file
write.fasta(seqs[asvs],names=asvs,file="family_x_seqs_for_blast.fasta",nbchar=500)

