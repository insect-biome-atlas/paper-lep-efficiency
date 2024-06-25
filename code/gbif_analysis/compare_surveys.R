
# Make lepidoptera maps from GBIF and IBA data -----------------------------------------------

# setup ---------------------------------------------------------------------------------------

# Load necessary libraries
library(tidyverse)       # Collection of R packages for data science
library(rnaturalearth)   # Provides map data for various countries
library(rnaturalearthdata) # Additional map data
library(sf)              # Simple features for R
library(patchwork)       # Combine multiple ggplot plots
library(stargazer)       # Create well-formatted regression tables

# Set a random seed for reproducibility
set.seed(10)

# data ----------------------------------------------------------------------------------------

# Read in IBA lepidoptera occurrence data
IBA_lep_occ <- readRDS("data/filtered_IBA_observations.rds") |> 
  select(Family, cluster, species = Species, lat = latitude_WGS84, lon = longitude_WGS84, trapID, sampleID_LAB, reads)

# Read in GBIF lepidoptera data
# GBIF.org (24 June 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.wyerdb
lep_data_GBIF <- readRDS("data/filtered_gbif_observations.rds") |> 
  select(Family = family, genus, species, lat = decimalLatitude, lon = decimalLongitude, year)

# Read in micro/macro classification of Lepidoptera families
lep_micro_macro <- read_delim("data/Lep_families_micro_macro.csv", delim = ";")

# Get map data for Sweden
sweden <- ne_countries(country = "sweden", scale = "large", returnclass = "sf")

# plots ---------------------------------------------------------------------------------------

# Prepare data for pooled richness for IBA
pr_IBA <- IBA_lep_occ |> drop_na() |> 
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(sweden)) |> 
  left_join(lep_micro_macro, by = "Family") |> select(species, Category) |> drop_na() |> distinct()

# Create a grid for the map
grid <- st_make_grid(pr_IBA, 1, crs = st_crs(pr_IBA), what = "polygons", square = FALSE, flat_topped = TRUE) |> 
  st_intersection(sweden) |> st_as_sf()

# useful indices ------------------------------------------------------------------------------

# Get unique species for GBIF data from 2019
lep_GBIF_2019_sp <- lep_data_GBIF |> 
  filter(year == 2019) |> 
  select(Family, species) |> 
  distinct()

# Get unique species for IBA data
lep_IBA_sp <- IBA_lep_occ |> 
  select(Family, species) |> 
  distinct()

# Find shared and unique species between GBIF and IBA data
IBA_sp    <- lep_IBA_sp$species
GBIF_sp   <- lep_GBIF_2019_sp$species
shared    <- intersect(IBA_sp, GBIF_sp)
IBA_only  <- setdiff(IBA_sp, GBIF_sp)
GBIF_only <- setdiff(GBIF_sp, IBA_sp)
total_sp  <- length(union(IBA_sp, GBIF_sp))

# Gbif ----------------------------------------------------------------------------------------

# Prepare GBIF data for mapping
lep_data_GBIF_2019 <- lep_data_GBIF |> 
                      filter(year == 2019) |> 
                      st_as_sf(coords = c("lon", "lat"), crs = st_crs(sweden)) |> 
                      left_join(lep_micro_macro, by = "Family") |> 
                      select(species, Category)

# Spatial join GBIF data with grid
cent_merge_gbif <- st_join(grid, lep_data_GBIF_2019, join = st_intersects, left = TRUE)


# Combine GBIF and IBA data
IBA_GBIF <- bind_rows(list(GBIF = lep_data_GBIF_2019, IBA = pr_IBA), .id = "id")

# Count species in GBIF data by category and grid cell
gbif_n_species <- cent_merge_gbif |> 
                  group_by(x, Category) |> 
                  count(species) |> 
                  select(-n)

# Set geometry column
st_geometry(gbif_n_species) <- "geometry"

# Combine and summarize species data
n_records_both <- gbif_n_species |>
                  group_by(geometry) |>
                  st_centroid() |> # Take only the centroid of each polygon in the gbif data
                  bind_rows(pr_IBA)

# Join grid with number of records in both
cent_merge_both <- st_join(grid, n_records_both, left = TRUE)

# Summarize macro and micro species data
macro_species_both <- cent_merge_both |> group_by(x, Category) |> 
  summarise(IBA = sum(species %in% IBA_only), GBIF = sum(species %in% GBIF_only)) |> 
  mutate(IBA = replace(IBA, IBA == 0, NA_real_)) |> 
  drop_na(Category) |> 
  mutate(total = IBA + GBIF, diff = (IBA - GBIF) / total) |> 
  select(diff, Category) |> 
  pivot_wider(values_from = diff, names_from = Category)

# Summarize total species data
n_species_both <- cent_merge_both |> group_by(x) |> 
  summarise(IBA = sum(species %in% IBA_only), GBIF = sum(species %in% GBIF_only)) |> 
  mutate(IBA = replace(IBA, IBA == 0, NA_real_)) |> 
  mutate(total = IBA + GBIF, n_species = (IBA - GBIF) / total) |> 
  select(n_species)

# Create a summary table for species data
n_species_tab <- cent_merge_both |> group_by(x) |> 
  summarise(IBA = sum(species %in% IBA_only), GBIF = sum(species %in% GBIF_only)) |> 
  rowwise() |> 
  mutate(total = sum(IBA, GBIF, na.rm = TRUE), IBA_proportion = IBA / total, GBIF_proportion = GBIF / total) |> 
  ungroup() |> 
  mutate(region_id = 1:n()) |> 
  st_drop_geometry()

# Save the summary table
saveRDS(n_species_tab, "figs/region_table.txt")

# difference ----------------------------------------------------------------------------------

# Total species

# Combine species data and summarize detections
all <- st_join(n_species_both, macro_species_both, join = st_equals) |> 
  pivot_longer(c("n_species", "MACRO", "micro"), names_to = "metric", values_to = "n_detections") |> 
  mutate(metric = factor(metric, labels = c("Total", "Macro", "Micro"), levels = c("n_species", "MACRO", "micro")))


# compare species detected bar charts ---------------------------------------------------------

# Get unique species for GBIF data from 2019
lep_GBIF_2019_sp <- lep_data_GBIF |> 
  filter(year == 2019) |> 
  select(Family, species) |> 
  distinct()

# Get unique species for IBA data
lep_IBA_sp <- IBA_lep_occ |> 
  select(Family, species) |> 
  distinct()

# Find shared and unique species between GBIF and IBA data
IBA_sp <- lep_IBA_sp$species
GBIF_sp <- lep_GBIF_2019_sp$species
shared <- intersect(IBA_sp, GBIF_sp)
IBA_only <- setdiff(IBA_sp, GBIF_sp)
GBIF_only <- setdiff(GBIF_sp, IBA_sp)

# Combine all species data and summarize by record type
all_sp <- bind_rows(lep_GBIF_2019_sp, lep_IBA_sp) |> distinct() |> 
  mutate(IBA_record = case_when(species %in% shared ~ "IBA_in_GBIF", species %in% IBA_only ~ "IBA_only", species %in% GBIF_only ~ "GBIF_only")) |> 
  group_by(IBA_record) |> 
  count(name = "Count")


# Combine all species data and summarize by family and record type
all_sp_family <- bind_rows(lep_GBIF_2019_sp, lep_IBA_sp) |> distinct() |> 
  left_join(lep_micro_macro, by = "Family") |> 
  mutate(IBA_record = case_when(species %in% shared ~ "IBA_in_GBIF", species %in% IBA_only ~ "IBA_only", species %in% GBIF_only ~ "GBIF_only"))

# figs -------------------------------------------------------------------------------------

# Create plots for species detected
p1 <- ggplot() +
  geom_sf(data = all, aes(fill = n_detections)) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "grey", breaks = c(-.99, -.5, 0, .5, .99)) +  
  labs(fill = "Species detected") +
  facet_wrap(~metric) +
  theme_linedraw(base_size = 15) +
  theme(legend.position = 'bottom', legend.title.position = "top", legend.text = element_text(size = 12), legend.title.align = 0.5, legend.key.width = unit(1, "cm"), strip.text = element_text(size = 20)) +
  scale_x_continuous(breaks = c(14, 18, 22))


# bar charts ----------------------------------------------------------------------------------

# Create bar chart for all species detected
p3 <- ggplot(all_sp, aes(x = "", y = Count, fill = IBA_record)) +
  geom_bar(stat = "identity") +
  labs(title = "a) All GBIF records", x = "", y = "Number of species") +
  scale_fill_manual(values = c("IBA_in_GBIF" = "ivory4", "GBIF_only" = "ivory2", "IBA_only" = "#7294D4"), labels = c("GBIF record not detected in metabarcoding", "GBIF record detected in metabarcoding", "Metabarcoding-specific OTUs")) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5)) +
  theme_minimal() +
  theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "vertical")

# Create bar chart for Macrolepidoptera
p4 <- all_sp_family |> filter(Category == "MACRO") |> 
  ggplot(aes(x = Family, fill = IBA_record)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "b) Macrolepidoptera", x = "Family", y = "Number of species") +
  scale_fill_manual(values = c("IBA_in_GBIF" = "ivory4", "GBIF_only" = "ivory2", "IBA_only" = "#7294D4"), labels = c("GBIF record not detected in metabarcoding", "GBIF record detected in metabarcoding", "Metabarcoding-specific OTUs")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  theme(legend.title = element_blank(), legend.position = "none")

# Create bar chart for Microlepidoptera
p5 <- all_sp_family |> filter(Category == "micro") |> 
  ggplot(aes(x = Family, fill = IBA_record)) +
  geom_bar() +
  theme_minimal() +
  labs(title = "c) Microlepidoptera", x = "Family", y = "Number of species") +
  scale_fill_manual(values = c("IBA_in_GBIF" = "ivory4", "GBIF_only" = "ivory2", "IBA_only" = "#7294D4"), labels = c("GBIF record not detected in metabarcoding", "GBIF record detected in metabarcoding", "Metabarcoding-specific OTUs")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  theme(legend.title = element_blank(), legend.position = "none")

# Save the bar charts
ggsave(filename = "figs/total_Counts_per_family.jpeg", device = "jpeg", width = 6, height = 8, p3)
ggsave(filename = "figs/MACROleps_Counts_per_family.jpeg", device = "jpeg", width = 5, height = 3, p4)
ggsave(filename = "figs/MICROleps_Counts_per_family.jpeg", device = "jpeg", width = 7, height = 4, p5)

# plots ---------------------------------------------------------------------------------------
library(patchwork)

# Save combined plot of range extensions
ggsave(filename = "figs/surv_comp.jpeg", device = "jpeg", width = 6, height = 8, p1)

# Open the saved plot
browseURL("figs/surv_comp.jpeg")
browseURL("figs/total_Counts_per_family.jpeg")
