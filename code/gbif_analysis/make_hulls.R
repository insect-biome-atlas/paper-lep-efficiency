#Construct alpha hulls of GBIF lepidoptera ------------------
# Packages

library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

set.seed(125)

# functions ---------------------------------------------------------------
hull_function <- function(species_locs){
  # species_locs = lat lon data for a single species
  sp_loc <- species_locs  # get all unique locations for a species
  hull <- chull(sp_loc$lat , sp_loc$lon)
  hull <- c(hull , hull[1])
  return(sp_loc[hull,]) # return alpha hull
}  # function to create a convex hull for all years of gbif occurrence records. 
hull_to_sf    <- function(hull , crsObj){ 
  #hull = convex hull for a single species
  #crs = crs object from st_crs
  st_as_sf(hull, coords = c("lon","lat") , crs = crsObj) %>%
    group_by(species , n) %>% 
    summarise(geometry=st_combine(geometry)) %>%
    st_cast("POLYGON")
} # function to turn hulls into sf spatial polygons

# data --------------------------------------------------------------------

# IBA lep occurrences with locations
IBA_lep_occ    <- readRDS("data/filtered_IBA_observations.rds") %>%
                  select(cluster , species = Species , 
                         lat = latitude_WGS84 , lon = longitude_WGS84 , trapID, sampleID_LAB , reads)

# Read in GBIF lepidoptera data
# GBIF.org (24 June 2024) GBIF Occurrence Download https://doi.org/10.15468/dl.wyerdb
swe_gbif_lep <- readRDS("data/filtered_gbif_observations.rds") %>% 
                select(genus , species , lat = decimalLatitude , lon = decimalLongitude , year)


# map of sweden + crs
sweden <- ne_countries(country = "sweden" , scale = "large" , returnclass = "sf")

# convert IBA lep occurrences to spatial data -----------------------------

# filter IBA occurrences to match gbif species & convert to spatial object
IBA_lep_sf <- IBA_lep_occ %>%
                distinct() %>% 
                drop_na() %>% 
                filter(species %in% unique(swe_gbif_lep$species)) %>% 
                st_as_sf(coords = c("lon" , "lat") , crs = st_crs(sweden)) %>% 
                arrange(species)

IBA_split <- IBA_lep_sf %>%  group_by(species , cluster) %>% group_split()

# alpha hull of gbif occurrences ------------------------------------------

# get a list of observations locations by species
swe_gbif_sub <- swe_gbif_lep %>% select(-year) %>% distinct() %>%               # remove duplicate observations at same coordinates (i.e. repeats in the same year) 
                  group_by(species) %>% mutate(n = n()) %>% filter(n >=3 ) %>%  # remove species where we can't calculate a hull
                  filter(species %in% IBA_lep_sf$species)  %>%                  # filter to IBA species
                  arrange(species)

# Split into lists by species
IBA_split <- IBA_lep_sf %>%
             filter(species %in% swe_gbif_sub$species) %>% # filter down to list of species we can calculate extents for
             group_by(species) %>% group_split()

gbif_split <- swe_gbif_sub %>% group_by(species) %>% group_split()


# convex hull of gbif species ---------------------------------------------

# calculate a hull for every species
all_hulls <- lapply(gbif_split , hull_function)

# convert hulls to spatial data 
all_hulls_sf <- lapply(all_hulls , hull_to_sf , crs = st_crs(sweden)) 

# calculate distance to known extents -------------------------------------

#  Loop over species & calculate distance between hull polygon and IBA observations
# -------------------------
distList <- vector(mode = "list" , length = length(IBA_split))
for(i in seq_along(IBA_split)){
  dists         <- st_distance(IBA_split[[i]] , all_hulls_sf[[i]]) 
  units(dists) <- NULL
  distList[[i]] <- data.frame(distance = dists , 
                              species  = IBA_split[[i]]$species[1], 
                              cluster  = IBA_split[[i]]$cluster[1],
                              reads = IBA_split[[i]]$reads,
                              sampleID_LAB = IBA_split[[i]]$sampleID_LAB,
                              trapID = IBA_split[[i]]$trapID,
                              listID = i) %>% 
                              st_set_geometry(IBA_split[[i]]$geometry)
 }
# -------------------------

allDist <- distList %>% bind_rows() 
rangeEx <- allDist %>% filter(distance > 0) %>% 
            mutate(lon = st_coordinates(.)[,1] , 
                   lat = st_coordinates(.)[,2]) %>% 
            st_drop_geometry() %>% 
            arrange(-distance)


# save preliminary data
write.csv(rangeEx , "data/prelim_IBA_range_extensions.csv")
saveRDS(all_hulls_sf , "data/all_hulls_sf.rds")
