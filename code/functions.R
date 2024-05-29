# Some useful functions for processing and
# extracting IBA data and for searching
# BOLD annotations


# Function for getting sample counts from asv table
sample_counts <- function(asvs, D) {

    x <- rep(0, times=ncol(D)-1)

    for (i in 1:length(asvs))
        x <- x + as.numeric(D[match(asvs[i],D$ASV_ID),2:ncol(D)])
    x
}

# Function for getting sample counts from cleaned cluster table
cluster_sample_counts <- function(clusters, D) {

    x <- rep(0, times=ncol(D)-1)

    for (i in 1:length(clusters))
        x <- x + as.numeric(D[match(clusters[i],D$cluster),2:ncol(D)])
    x
}

# Function for getting BOLD annotations
# On failure to find the BOLD id, taxon will be null
# so this can be tested with 'is.null(return_value$taxon)'
library("bold")
bold_annotation <- function(x) {

    res <- bold_stats(bin=x)

    tot_records <- res$total_records
    taxon <- res$species$drill_down$entity$name[1]
    n_records <- res$species$drill_down$entity$records[1]

    list(taxon=taxon, n_records=n_records, prop=n_records/tot_records)
}

