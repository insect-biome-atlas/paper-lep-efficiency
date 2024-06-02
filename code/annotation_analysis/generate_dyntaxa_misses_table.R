# Set path to data
data_path <- "~/dev/IBA_data_April2023/"

# Set path to dyntaxa response file
dyntaxa_response_path <- "dyntaxa_searches/"

# Define a useful bold search function
library("bold")

# Get useful function for bold annotation
source("functions.R")

# Read in the nochimera cluster info
nochimera_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_nochimera_clusters_SE_2019.tsv"),sep="\t",header=TRUE)

# Read in the resolved taxonomy info for cleaned clusters
cleaned_cluster_taxonomy <- read.table(paste0(data_path,"data/Sweden/CO1_cleaned_nochimera_cluster_taxonomy_SE_2019.tsv"),sep="\t",header=TRUE)
colnames(cleaned_cluster_taxonomy)[1]<-"cluster"

# Read in the dyntaxa response and rename the columns
dyntaxa_response <- read.table(paste0(dyntaxa_response_path,"lep_species_dyntaxa_response.csv"),header=TRUE,sep=";")
colnames(dyntaxa_response) <- c("IBA_species","match_status","taxon_id","sci_name","author","synonyms","SE_occurrence")

# Extract the entries not matching a species reproducing in Sweden ("Bofast och reproducerande")
dyntaxa_misses <- dyntaxa_response[dyntaxa_response$SE_occurrence!="Bofast och reproducerande",]

# Add useful information about the clusters
len <- nrow(dyntaxa_misses)
dyntaxa_misses$IBA_cluster <- character(len)
dyntaxa_misses$IBA_resolved_species <- character(len)
dyntaxa_misses$IBA_resolved_BOLD_bin <- character(len)
dyntaxa_misses$updated_BIN_annotation <- character(len)
dyntaxa_misses$updated_BIN_annotation_n_records <- numeric(len)
dyntaxa_misses$updated_BIN_annotation_prop <- numeric(len)

# Add cluster info (print warning if more than one)
for (i in 1:len) {

    species <- dyntaxa_misses$IBA_species[i]
    cat ("Getting info for: ",species,"\n")

    x <- length(unique(nochimera_clusters$cluster[nochimera_clusters$Species==species]))
    if (x > 1)
        cat("Warning: ", x, " matches for species: ", species, "\n")

    dyntaxa_misses$IBA_cluster[i] <- nochimera_clusters$cluster[match(species,nochimera_clusters$Species)]
    if (dyntaxa_misses$IBA_cluster[i] %in% cleaned_cluster_taxonomy$cluster) {

        m <- match(dyntaxa_misses$IBA_cluster[i],cleaned_cluster_taxonomy$cluster)
        dyntaxa_misses$IBA_resolved_species[i] <- cleaned_cluster_taxonomy$Species[m]
        dyntaxa_misses$IBA_resolved_BOLD_bin[i] <- cleaned_cluster_taxonomy$BOLD_bin[m]
        if (grepl("BOL",dyntaxa_misses$IBA_resolved_BOLD_bin[i])) {
            bold_full <- paste0("BOLD:",substring(dyntaxa_misses$IBA_resolved_BOLD_bin[i],first=4))
            res <- bold_annotation(bold_full)
            if (!is.null(res$taxon)) {
                dyntaxa_misses$updated_BIN_annotation[i] <- res$taxon
                dyntaxa_misses$updated_BIN_annotation_n_records[i] <- res$n_records
                dyntaxa_misses$updated_BIN_annotation_prop[i] <- res$prop
            }
        }
    }
}

# Write table
write.table (dyntaxa_misses, file="../../output_tables/dyntaxa_misses.tsv")

