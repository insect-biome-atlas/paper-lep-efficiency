# Set path to data files
data_path <- "~/dev/IBA_data_April2023/"

# Set path to blast response files
blast_response_path <- "blast_searches/"

# Read in base info in table
cleaned_cluster_taxonomy <- read.table(paste0(data_path,"data/Sweden/CO1_cleaned_nochimera_cluster_taxonomy_SE_2019.tsv"),sep="\t",header=TRUE)
colnames(cleaned_cluster_taxonomy)[1]<-"cluster"
lep_cleaned_cluster_taxonomy <- cleaned_cluster_taxonomy[cleaned_cluster_taxonomy$Order=="Lepidoptera",]
unidentified <- lep_cleaned_cluster_taxonomy[!grepl(" ",lep_cleaned_cluster_taxonomy$Species) & !grepl("BOL",lep_cleaned_cluster_taxonomy$BOLD_bin),]

# Read in base info in table
# Get the rep asv info
nochimera_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_nochimera_clusters_SE_2019.tsv"),sep="\t",header=TRUE)
nochimera_cluster_reps <- nochimera_clusters[nochimera_clusters$representative==1,]
                                       
# Read in base info in table
# Match to the ufo clusters
unidentified$rep_asv <- nochimera_cluster_reps$ASV[match(unidentified$cluster,nochimera_cluster_reps$cluster)]

cat("Reding the cluster counts. This may take a while...\n")

# Read in base info in table
# Read in the cluster counts
cluster_counts <- read.table(paste0(data_path,"data/Sweden/CO1_cleaned_nochimera_cluster_counts_SE_2019.tsv"),sep="\t",header=TRUE)

# Read in table of blast response file names
blast_response_files <- read.table("blast_response_files.tsv")

# Initialize the columns for the info we are interested in
len <- nrow(unidentified)
unidentified$blast_taxon <- character(len)
unidentified$blast_idty <- numeric(len)
unidentified$blast_e_value <- numeric(len)
unidentified$mitochondrial <- logical(len)
unidentified$n_samples <- numeric(len)
unidentified$n_reads <- numeric(len)
unidentified$iba_match <- logical(len)
unidentified$n_samples_match <- numeric(len)
unidentified$n_reads_match <- numeric(len)
unidentified$n_samples_overlap <- numeric(len)
unidentified$prop_overlap <- numeric(len)

# ... and cycle through the blast hits for the ufo asvs
# filling in the relevant info
for (i in 1:nrow(unidentified)) {

    cluster <- unidentified$cluster[i]
    asv <- unidentified$rep_asv[i]

    # Extract file name for corresponding blast hit description table
    file_name <- blast_response_files$file_name[match(asv,blast_response_files$asv)]
    file_name <- paste0(blast_response_path,file_name)

    temp <- read.csv(file_name)

    blast_taxon <- temp$Scientific.Name[1] 
    cat("Getting blast hit info for: ", blast_taxon,"\n")

    unidentified$blast_taxon[i] <- blast_taxon
    unidentified$blast_idty[i] <- temp$Per..ident[1]
    unidentified$blast_e_value[i] <- temp$E.value[1]

    unidentified$mitochondrial[i] <- grepl("mitochondrial",temp$Description[1])

    x <- as.numeric(cluster_counts[match(cluster,cluster_counts$cluster),2:ncol(cluster_counts)])
    unidentified$n_samples[i] <- sum(x!=0)
    unidentified$n_reads[i] <- sum(x)

    if (blast_taxon %in% lep_cleaned_cluster_taxonomy$Species) {

        unidentified$iba_match[i] <- TRUE
        
        y <- rep(0,times=length(x))
        for (c in lep_cleaned_cluster_taxonomy$cluster[lep_cleaned_cluster_taxonomy$Species==blast_taxon]) {
            y <- y + as.numeric(cluster_counts[match(c,cluster_counts$cluster),2:ncol(cluster_counts)])
        }

        unidentified$n_samples_match[i] <- sum(y!=0)
        unidentified$n_reads_match[i] <- sum(y)
        unidentified$n_samples_overlap[i] <- sum(y&x)
        unidentified$prop_overlap[i] <- sum(y&x) / sum(x!=0)
    }

    else
        unidentified$iba_match[i] <- FALSE

}

# Write resulting table
write.table(unidentified, "../output_tables/unidentified.tsv")

