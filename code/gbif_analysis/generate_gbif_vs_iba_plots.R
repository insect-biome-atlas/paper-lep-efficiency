# Script for generating figures 3, S1 and S2
# Also computes some numbers cited in text

# Define groups of Lepidoptera families
leps <- c('Adelidae','Alucitidae','Argyresthiidae','Autostichidae','Batrachedridae','Bedelliidae','Blastobasidae','Brahmaeidae','Bucculatricidae','Castniidae','Chimabachidae','Choreutidae','Coleophoridae','Cosmopterigidae','Cossidae','Crambidae','Depressariidae','Douglasiidae','Drepanidae','Dryadaulidae','Elachistidae','Endromidae','Epermeniidae','Erebidae','Eriocraniidae','Ethmiidae','Gelechiidae','Geometridae','Glyphipterigidae','Gracillariidae','Heliodinidae','Heliozelidae','Hepialidae','Hesperiidae','Incurvariidae','Lasiocampidae','Limacodidae','Lycaenidae','Lyonetiidae','Lypusidae','Meessiidae','Micropterigidae','Momphidae','Nepticulidae','Noctuidae','Nolidae','Notodontidae','Nymphalidae','Oecophoridae','Opostegidae','Papilionidae','Parametriotidae','Peleopodidae','Pieridae','Plutellidae','Praydidae','Prodoxidae','Psychidae','Pterophoridae','Pyralidae','Riodinidae','Roeslerstammiidae','Saturniidae','Schreckensteiniidae','Scythrididae','Scythropiidae','Sesiidae','Sphingidae','Stathmopodidae','Tineidae','Tischeriidae','Tortricidae','Urodidae','Yponomeutidae','Ypsolophidae','Zygaenidae')

butterflies <- c("Pieridae","Nymphalidae","Lycaenidae","Hesperiidae","Papilionidae","Riodinidae")

macroleps <- c('Brahmaeidae','Cossidae','Drepanidae','Endromidae','Erebidae','Geometridae','Hepialidae','Hesperiidae','Lasiocampidae','Limacodidae','Lycaenidae','Noctuidae','Nolidae','Notodontidae','Nymphalidae','Papilionidae','Pieridae','Psychidae','Riodinidae','Saturniidae','Sesiidae','Sphingidae','Zygaenidae')

macros_no_butterflies <- macroleps[!(macroleps %in% butterflies)]

microleps <- leps[!(leps %in% macroleps)]


# Set path to local copy of figshare repo for the IBA lep efficiency paper
path_gbif <- "~/dev/figshare-repos/iba/paper-lep-efficiency/"

# Set path to local copy of figshare repo 2 (cleaned data) for the IBA data paper, Swedish data
# path_iba <- "~/dev/figshare-repos/iba/paper-cleaned-data/data/Sweden/"
path_iba <- "~/Documents/Manuscripts/2023_IBA_data_paper/git-uppmax_IBA_data_April2023/data/Sweden/"

# Set path to ms repo for the IBA lep efficiency paper
path_lep <- "~/dev/ms-repos-iba/paper-lep-efficiency/output_tables/"

# Read in IBA cleaned cluster counts
cat ("Reading in cluster counts (this will take a while)\n")
cluster_counts <- read.delim(paste0(path_iba, "CO1_cleaned_nochimera_cluster_counts_SE_2019.tsv"),sep="\t",header=TRUE)

# Read in IBA cleaned cluster resolved taxonomy info
cluster_taxonomy <- read.delim(paste0(path_iba, "CO1_cleaned_nochimera_cluster_taxonomy_SE_2019.tsv"),sep="\t",header=TRUE)

# Fix problem with first column name
colnames(cluster_taxonomy)[1] <- "cluster"

# Get manual taxonomy annotations
man_taxonomy <- read.table(paste0(path_lep,"summary_manual_taxonomic_annotations.tsv"))

# Get all GBIF occurrences for Lepidoptera species
cat ("Reading in GBIF lep data\n")
gbif_lep_occurrences <- read.delim(paste0(path_gbif,"0055581-240506114902167.csv"),header=TRUE,sep="\t")
gbif_lep_occurrences <- gbif_lep_occurrences[gbif_lep_occurrences$taxonRank=="SPECIES" & gbif_lep_occurrences$taxonomicStatus=="ACCEPTED",]

# Get GBIF insect data for 2019 (only Artportalen in this version)
cat ("Reading in GBIF 2019 data (this will take a while)\n")
gbif_2019_data <- read.delim(paste0(path_gbif,"0055548-240506114902167.csv"),header=TRUE,sep="\t")




# Plot Fig. S1
# ============
#
# We want to compare GBIF occurrences 
# for IBA and non-IBA species of leps

# First get updated lep_species
lep_clusters <- cluster_taxonomy[cluster_taxonomy$Order=="Lepidoptera",]
lep_clusters$Species_updated <- lep_clusters$Species
for (i in 1:nrow(lep_clusters)) {
    if (lep_clusters$cluster[i] %in% man_taxonomy$cluster)
        lep_clusters$Species_updated[i] <- man_taxonomy$gbif_name[match(lep_clusters$cluster[i],man_taxonomy$cluster)]
}

gbif_lep_occurrences$in_iba <- gbif_lep_occurrences$species %in% lep_clusters$Species_updated
gbif_lep_occurrences$label <- "In IBA"
for (i in 1:nrow(gbif_lep_occurrences)) {
    if (!gbif_lep_occurrences$in_iba[i])
        gbif_lep_occurrences$label[i] <- "Not in IBA"
}

cat("Plotting Fig. S1\n")
pdf(file="../../figs/Fig_S1.pdf")
boxplot(gbif_lep_occurrences$numberOfOccurrences~gbif_lep_occurrences$label,ylab="Occurrence records",log="y",xlab="")
dev.off()

print(t.test(log(numberOfOccurrences)~label, data = gbif_lep_occurrences))



# Plot Fig. S2
# ============
#
# We want to compare GBIF occurrences
# and IBA occurrences in 2019 for the 
# major insect groups

# Get GBIF 2019 data for major groups
x <- table(gbif_2019_data$order)
lep <- x[["Lepidoptera"]]
col <- x[["Coleoptera"]]
hym <- x[["Hymenoptera"]]
dip <- x[["Diptera"]]
hem <- x[["Hemiptera"]]
other <- nrow(gbif_2019_data)-lep-col-hym-dip-hem

# Compute number of occurrence records for each iba cluster
cluster_occurrences <- numeric(nrow(cluster_counts))
cat("Computing occurrence records\n")
A <- as.matrix(cluster_counts[,2:ncol(cluster_counts)])
B <- A!=0
cluster_occurrences <- rowSums(B)

# Add the lep data to the lep_clusters dataframe generated above for future use
lep_clusters$occurrences <- cluster_occurrences[match(lep_clusters$cluster,cluster_counts$cluster)]

# Summarize the IBA occurrences by major group 
cat("Summarizing IBA occurrence records\n")
lep[2] <- sum(cluster_occurrences[cluster_taxonomy$Order[match(cluster_counts$cluster,cluster_taxonomy$cluster)]=="Lepidoptera"])
col[2] <- sum(cluster_occurrences[cluster_taxonomy$Order[match(cluster_counts$cluster,cluster_taxonomy$cluster)]=="Coleoptera"])
hym[2] <- sum(cluster_occurrences[cluster_taxonomy$Order[match(cluster_counts$cluster,cluster_taxonomy$cluster)]=="Hymenoptera"])
dip[2] <- sum(cluster_occurrences[cluster_taxonomy$Order[match(cluster_counts$cluster,cluster_taxonomy$cluster)]=="Diptera"])
hem[2] <- sum(cluster_occurrences[cluster_taxonomy$Order[match(cluster_counts$cluster,cluster_taxonomy$cluster)]=="Hemiptera"])
other[2] <- sum(cluster_occurrences)-lep[2]-col[2]-hym[2]-dip[2]-hem[2]

order_stats <- data.frame(list(Lepidoptera=lep,Coleoptera=col,Hymenoptera=hym,Diptera=dip,Hemiptera=hem,Other=other))
rownames(order_stats) <- c("GBIF","IBA")

cat("Plotting Fig. S2\n")
pdf(file="../../figs/Fig_S2.pdf")
barplot(as.matrix(order_stats), beside=TRUE,legend.text=TRUE,cex.names=0.8,ylab="Occurrence records")
dev.off()



# Plot Fig. 3
# ============
#
# We want to plot GBIF occurrences
# against IBA occurrences in 2019 for 
# families and species

# Get GBIF family data
gbif_fams <- data.frame(table(gbif_2019_data$family[gbif_2019_data$order=="Lepidoptera"]))
colnames(gbif_fams) <- c("family","occurrences")

# Get GBIF species data
gbif_species <- data.frame(table(gbif_2019_data$species[gbif_2019_data$order=="Lepidoptera"]))
colnames(gbif_species) <- c("species","occurrences")
gbif_species$family <- gbif_2019_data$family[match(gbif_species$species,gbif_2019_data$species)]

# Get IBA and GBIF data for families into same dataframe
fam_stat <- data.frame(xtabs(lep_clusters$occurrences~lep_clusters$Family))
colnames(fam_stat) <- c("family","iba_occurrences")
fam_stat$gbif_occurrences <- gbif_fams$occurrences[match(fam_stat$family,gbif_fams$family)]

# Get number of species per family
dyntaxa <- read.delim(sep=";",file="../../data/DynTaxa_Lepidoptera_2024-03-19_19.04.csv")
fam_diversity <- data.frame(table(dyntaxa$Familj[dyntaxa$Taxonkategori=="Art" & dyntaxa$Taxonstatus=="Accepterat"]))
colnames(fam_diversity) <- c("family","species_count")
fam_stat$species_count <- fam_diversity$species_count[match(fam_stat$family,fam_diversity$family)]
fam_stat$iba_rel_occ <- fam_stat$iba_occurrences / fam_stat$species_count
fam_stat$gbif_rel_occ <- fam_stat$gbif_occurrences / fam_stat$species_count

# Add column summarizing if this is a microlep, butterfly of other macrolep.
fam_stat$group <- "micros"
fam_stat[fam_stat$family %in% macros_no_butterflies, "group"] <- "other_macros"
fam_stat[fam_stat$family %in% butterflies, "group"] <- "butterflies"

f3A <-
  ggplot(fam_stat, aes(x=iba_rel_occ, y=gbif_rel_occ, color=group)) + 
  geom_point(size=3) + 
  geom_abline() +
  scale_x_continuous(trans='log2', breaks = c(0.2,0.5,1,2,5,10, 20,50, 100)) +
  scale_y_continuous(trans='log2', breaks = c(2,5,10,20,50,100,200,500, 1000)) +
  scale_color_manual(values=c( "tomato2", "#999999","dodgerblue3")) +
  labs(y = "GBIF occurrence records per species [log]", x="IBA occurrence records per species [log]") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
f3A
ggsave(filename = "./../../figs/Fig_3A_EIE.pdf", device="pdf", width = 6, height = 5.5, f3A)

# Get IBA and GBIF data for species into same dataframe
# Note that we rely on species name matching to remove IBA clusters that are not identified to species
# We ignore the duplicated species clusters here; they are represented more than one time in the sp_stat file
# To remove this problem entirely, clusters have to be aggregated in the computation of the occurrence records above
# Note that there are similar problems with lumped species, which cannot be resolved easily
# Overall, these data points represent a small fraction of the total number of data points
lep_clusters$gbif_occurrences <- gbif_species$occurrences[match(lep_clusters$Species_updated,gbif_species$species)]
sp_stat <- lep_clusters[!is.na(lep_clusters$gbif_occurrences),]

# Add column summarizing if this is a microlep, butterfly of other macrolep.
sp_stat$group <- "micros"
sp_stat[sp_stat$Family %in% macros_no_butterflies, "group"] <- "other_macros"
sp_stat[sp_stat$Family %in% butterflies, "group"] <- "butterflies"

f3B <-
  ggplot(sp_stat, aes(x=occurrences, y=gbif_occurrences, color=group)) + 
  geom_point(position = position_jitterdodge(0.3, dodge.width = .1), size=2) + 
  geom_abline() +
  scale_x_continuous(trans='log2',
                     breaks = c(1,5,10,50,100,500, 1000)) +
  scale_y_continuous(trans='log2',
                     breaks = c(1,10,100,1000,10000 )) +
  scale_color_manual(values=c( "tomato2", "#999999","dodgerblue3")) +
  labs(y = "GBIF occurrence records [log]", x="IBA occurrence records [log]") +
  theme_bw() +
  theme(legend.position="none", axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))

ggsave(filename = "./../../figs/Fig_3B_EIE.pdf", device="pdf", width = 6, height = 5.5, f3B)

cat("Pearson correlation coefficient (log(sp_stat$gbif_occurrences),log(sp_stat$occurrences))\n")
print(cor(log(sp_stat$gbif_occurrences),log(sp_stat$occurrences)))

# Compute number of recorded species and unique species (this is ignorant to whether the species represent new records)
cat("Total number of shared and unique species\n")
cat("Number of lep species recorded in GBIF in 2019:", length(unique(gbif_species$species[grepl(" ",gbif_species$species) & gbif_species$family!="Gerridae"])),"\n")
cat("Number of lep species recorded in IBA in 2019:", length(unique(lep_clusters$Species_updated[grepl(" ",lep_clusters$Species_updated)])), "\n")
cat("Number of shared lep species recorded in 2019:", length(unique(lep_clusters$Species_updated[lep_clusters$Species_updated %in% gbif_species$species])),"\n")

