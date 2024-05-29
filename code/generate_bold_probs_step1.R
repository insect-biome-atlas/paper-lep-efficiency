# Set path to data
data_path <- "~/dev/IBA_data_April2023/"

# Define a useful bold search function
library("bold")

# Read in useful functions
source("functions.R")

# Read in the nochimera cluster info
nochimera_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_nochimera_clusters_SE_2019.tsv"),sep="\t",header=TRUE)

# Read in the resolved taxonomy info for cleaned clusters
cleaned_cluster_taxonomy <- read.table(paste0(data_path,"data/Sweden/CO1_cleaned_nochimera_cluster_taxonomy_SE_2019.tsv"),sep="\t",header=TRUE)
colnames(cleaned_cluster_taxonomy)[1]<-"cluster"

# Create the basic table info
attach(cleaned_cluster_taxonomy)
bold_probs <- cleaned_cluster_taxonomy[Order=="Lepidoptera" & grepl("BOL",BOLD_bin) & !grepl(" ",Species),]
detach(cleaned_cluster_taxonomy)

# Add updated bold annotations
len <- nrow(bold_probs)
updated_BIN_annotation <- character(len)
updated_BIN_annotation_n_records <- numeric(len)
updated_BIN_annotation_prop <- numeric(len)

for (i in 1:len) {

    bold_id <- paste0("BOLD:",substring(bold_probs$BOLD_bin[i],first=4))

    cat ("Updating annotation for: ", bold_id, "\n")

    res <- bold_annotation(bold_id)

    if (!is.null(res$taxon)) {
        updated_BIN_annotation[i] <- res$taxon
        updated_BIN_annotation_n_records[i] <- res$n_records
        updated_BIN_annotation_prop[i] <- res$prop
    }
}

bold_probs$updated_BIN_annotation <- updated_BIN_annotation
bold_probs$updated_BIN_annotation_n_records <- updated_BIN_annotation_n_records
bold_probs$updated_BIN_annotation_prop <- updated_BIN_annotation_prop

# Write table and file for dyntaxa matching
write.table(bold_probs,file="bold_probs_step1.tsv")
cat(unique(bold_probs$updated_BIN_annotation[!is.na(bold_probs$updated_BIN_annotation)]),sep=";",file="bold_probs_for_dyntaxa.txt")

