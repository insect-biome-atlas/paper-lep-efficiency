# Set path to data
data_path <- "~/dev/IBA_data_April2023/"

# Read in the lep species
nochimera_clusters <- read.table(paste0(data_path,"data/Sweden/CO1_nochimera_clusters_SE_2019.tsv"),sep="\t",header=TRUE)
lep_species <- unique(nochimera_clusters$Species[nochimera_clusters$Order=="Lepidoptera" & grepl(" ",nochimera_clusters$Species)])

# Write the species names in a format that can be used for online search in dyntaxa
cat (lep_species,sep=";",file="lep_species_for_dyntaxa.txt")

