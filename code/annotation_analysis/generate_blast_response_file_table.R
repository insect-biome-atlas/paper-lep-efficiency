# This script generates a table of the file names for the blast responses
# to the unidentified, family_x and unclassified.family sets of lep clusters

# Set path to data
data_path <- "~/dev/IBA_data_April2023/"

# Read in base info in ufos table
cleaned_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_cleaned_nochimera_cluster_taxonomy_SE_2019.tsv"),sep="\t",header=TRUE)
colnames(cleaned_clusters)[1]<-"cluster"
ufos <- cleaned_clusters[cleaned_clusters$Order=="Lepidoptera" & !grepl(" ",cleaned_clusters$Species) & !grepl("BOL",cleaned_clusters$BOLD_bin),]

# Read in and extract the rep asv info
nochimera_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_nochimera_clusters_SE_2019.tsv"),sep="\t",header=TRUE)
nochimera_cluster_reps <- nochimera_clusters[nochimera_clusters$representative==1,]
 
# Match to the ufo clusters
ufos$rep_asv <- nochimera_cluster_reps$ASV[match(unidentified$cluster,nochimera_cluster_reps$cluster)]

# Specify blast search id
search_id <- "8E8XF8GY016"

# Initialize data frame variables
asv <- character()
cluster <- character()
file_name <- character()

# Generate response file table variables
for (i in 1:nrow(ufos)) {

    asv[i] <- ufos$rep_asv[i]
    cluster[i] <- ufos$cluster[i]
    file_name[i] <- paste0(search_id,"-Alignment-Descriptions(",i,").csv")
}

# Extract the "Lepidoptera_x" clusters
family_x <- nochimera_clusters[nochimera_clusters$Family=="Lepidoptera_X" & nochimera_clusters$representative==1,]
colnames(family_x)[1] <- "rep_asv"

# Specify blast search id
search_id <- "8GJEE4XV016"

# Extend response file table variables
j <- length(asv)
for (i in 1:nrow(family_x)) {

    asv[j+i] <- family_x$rep_asv[i]
    cluster[j+i] <- family_x$cluster[i]
    file_name[j+i] <- paste0(search_id,"-Alignment-Descriptions(",i,").csv")
}

# Extract the "unclassified.Lepidoptera" clusters
unclassified.family <- nochimera_clusters[nochimera_clusters$Family=="unclassified.Lepidoptera" & nochimera_clusters$representative==1,]
colnames(unclassified.family)[1] <- "rep_asv"

# Specify blast search ids (search divided into 10 batches)
search_id <-c("8GMN8PEH013","8GNGZJNE013","8GNTKHYS013","8GP03A6U016","8GP51B9F016","8GPBDEPP013")

# Extend response file table variables
j <- length(asv)
for (i in 1:nrow(unclassified.family)) {

    asv[j+i] <- unclassified.family$rep_asv[i]
    cluster[j+i] <- unclassified.family$cluster[i]

    # Construct file name for corresponding blast hit description table
    # taking batch id into account
    batch <- floor((i-1)/10) + 1
    file_name[j+i] <- paste0(search_id[batch],"-Alignment-Descriptions(",i-10*(batch-1),").csv")
}

blast_response_files <- data.frame(asv=asv,cluster=cluster,file_name=file_name)
write.table(blast_response_files,"blast_response_files.tsv")

