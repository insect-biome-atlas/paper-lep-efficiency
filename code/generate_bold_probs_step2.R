# Set path to dyntaxa response file
dyntaxa_response_path <- "dyntaxa_searches/"

# Read in data frame from step 1
bold_probs <- read.table("bold_probs_step1.tsv")

# Read in the dyntaxa response and rename the columns
dyntaxa_response <- read.table(paste0(dyntaxa_response_path,"bold_probs_dyntaxa_response.csv"),header=TRUE,sep=";")
colnames(dyntaxa_response) <- c("IBA_species","match_status","taxon_id","sci_name","author","synonyms","SE_occurrence")

# Complement the table with dyntaxa info
len <- nrow(bold_probs)
bold_probs$dyntaxa_match <- character(len)
bold_probs$dyntaxa_name <- character(len)
bold_probs$SE_occurrence <- character(len)

for (i in 1:nrow(bold_probs)) {

    m <- match(bold_probs$updated_BIN_annotation[i],dyntaxa_response$IBA_species)

    if (!is.na(m)) {
        bold_probs$dyntaxa_match[i] <- dyntaxa_response$match_status[m]
        bold_probs$dyntaxa_name[i] <- dyntaxa_response$sci_name[m]
        bold_probs$SE_occurrence[i] <- dyntaxa_response$SE_occurrence[m]
    }
}
 
# Write table
write.table(bold_probs,file="../output_tables/bold_probs.tsv")

