library("seqinr")
data_path <- "~/dev/IBA_data_April2023/"

# Read in nochimera cluster data
nochimera_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_nochimera_clusters_SE_2019.tsv"),sep="\t",header=TRUE)

# Read in cleaned cluster data
cleaned_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_cleaned_nochimera_cluster_taxonomy_SE_2019.tsv"),sep="\t",header=TRUE)
colnames(cleaned_clusters)[1] <- "cluster"

# Read in sequence data
cat("Reading asv sequences. This may take a while...\n")
seqs <- read.fasta(paste0(data_path,"data/Sweden/CO1_asv_seqs_SE_2019.fasta"))

# Extract data for lep cluster representative ASVs
lep_cluster_reps <- nochimera_clusters[nochimera_clusters$Order=="Lepidoptera" & nochimera_clusters$representative==1,]

# Extract the unidentified lep clusters
lep_unidentified_clusters <- cleaned_clusters$cluster[cleaned_clusters$Order=="Lepidoptera" & !grepl(" ",cleaned_clusters$Species) & !grepl("BOL",cleaned_clusters$BOLD_bin)]

# Extract the asvs for the clusters we are interested in
asvs <- lep_cluster_reps$ASV[match(lep_unidentified_clusters,lep_cluster_reps$cluster)]

# Write the fasta file
write.fasta(seqs[asvs],names=asvs,file="unidentified_seqs_for_blast.fasta",nbchar=500)


