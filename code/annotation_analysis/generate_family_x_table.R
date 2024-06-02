# Enable the reading in of the asv counts table
library(data.table)

# Set path to data
data_path <- "~/dev/IBA_data_April2023/"

# Set path to blast responnse files
blast_response_path <- "blast_searches/"

# Read in useful functions
source("functions.R")

# Read in base info in table
nochimera_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_nochimera_clusters_SE_2019.tsv"),sep="\t",header=TRUE)
family_x <- nochimera_clusters[nochimera_clusters$Family=="Lepidoptera_X" & nochimera_clusters$representative==1,]

# Read in the asv counts (requires a lot of memory)
cat("Reading the asv counts. This may take a while...\n")
asv_counts <- fread(paste0(data_path,"data/Sweden/CO1_asv_counts_SE_2019.tsv"),sep="\t",header=TRUE)

# Read in the cluster counts (slow)
cat("Reading the cluster counts. This may take a while...\n")
cluster_counts <- fread(paste0(data_path,"data/Sweden/CO1_cleaned_nochimera_cluster_counts_SE_2019.tsv"),sep="\t",header=TRUE)

# Read in the lepidoptera cluster info
cleaned_cluster_taxonomy <- read.table(paste0(data_path,"data/Sweden/CO1_cleaned_nochimera_cluster_taxonomy_SE_2019.tsv"),sep="\t",header=TRUE)
colnames(cleaned_cluster_taxonomy)[1]<-"cluster"
lep_cleaned_cluster_taxonomy <- cleaned_cluster_taxonomy[cleaned_cluster_taxonomy$Order=="Lepidoptera",]

# Read in the blast response file table
blast_response_files <- read.table("blast_response_files.tsv")

# Initialize the columns for the info we are interested in
len <- nrow(family_x)
colnames(family_x)[1] <- "rep_asv"
family_x$blast_taxon <- character(len)
family_x$blast_idty <- numeric(len)
family_x$blast_e_value <- numeric(len)
family_x$mitochondrial <- logical(len)
family_x$n_samples <- numeric(len)
family_x$n_reads <- numeric(len)
family_x$iba_match <- logical(len)
family_x$n_samples_match <- numeric(len)
family_x$n_reads_match <- numeric(len)
family_x$n_samples_overlap <- numeric(len)
family_x$prop_overlap <- numeric(len)

# ... and cycle through the blast hits for the ufo asvs
# filling in the relevant info
for (i in 1:nrow(family_x)) {

    cluster <- family_x$cluster[i]
    asv <- family_x$rep_asv[i]

    # Get the file name for corresponding blast hit description table
    file_name <- blast_response_files$file_name[match(asv,blast_response_files$asv)]
    file_name <- paste0(blast_response_path,file_name)

    temp <- read.csv(file_name)

    blast_taxon <- temp$Scientific.Name[1] 
    cat("Getting blast hit info for: ", blast_taxon,"\n")

    family_x$blast_taxon[i] <- blast_taxon
    family_x$blast_idty[i] <- temp$Per..ident[1]
    family_x$blast_e_value[i] <- temp$E.value[1]

    family_x$mitochondrial[i] <- grepl("mitochondrial",temp$Description[1])

    x <- sample_counts(nochimera_clusters$ASV[nochimera_clusters$cluster==cluster],asv_counts)
    family_x$n_samples[i] <- sum(x!=0)
    family_x$n_reads[i] <- sum(x)

    if (blast_taxon %in% lep_cleaned_cluster_taxonomy$Species) {

        family_x$iba_match[i] <- TRUE
        
        y <- cluster_sample_counts(lep_cleaned_cluster_taxonomy$cluster[lep_cleaned_cluster_taxonomy$Species==blast_taxon],cluster_counts)

        family_x$n_samples_match[i] <- sum(y!=0)
        family_x$n_reads_match[i] <- sum(y)
        family_x$n_samples_overlap[i] <- sum(y&x)
        family_x$prop_overlap[i] <- sum(y&x) / sum(x!=0)
    }

    else
        family_x$iba_match[i] <- FALSE

}

# Write resulting table
write.table(family_x, "../output_tables/family_x.tsv")

