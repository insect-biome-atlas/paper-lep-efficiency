library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)


# data --------------------------------------------------------------------

all_dists    <- read.csv("data/prelim_IBA_range_extensions.csv")
all_hulls_sf <- readRDS("data/all_hulls_sf.rds")
swe_gbif_lep <- readRDS("data/filtered_gbif_observations.rds") %>% 
                select(genus , species , lat = decimalLatitude , lon = decimalLongitude , year)

sweden       <- ne_countries(country = "sweden" , scale = "large" , returnclass = "sf")

leps_to_plot <- read_delim("data/leps_to_map.csv" , delim = ";")
IBA_lep_occ  <- readRDS("data/filtered_IBA_observations.rds") %>%
                select(cluster , species = Species , lat = latitude_WGS84 , lon = longitude_WGS84 , trapID) |> 
                distinct() |> 
                right_join(leps_to_plot)


# functions -----------------------------------------------------------------------------------
plot_fun <- function(category){

catObs   <- filter(IBA_lep_occ ,CATEGORY == category)  
spec_obs <- filter(all_dists ,species %in% catObs$species)
hull_df  <- all_hulls_sf %>% bind_rows() %>% filter(species %in% catObs$species)
obs_lep  <- swe_gbif_lep %>% filter(species %in% catObs$species)

p<-  ggplot()+
  geom_sf(data=sweden)+
  geom_sf(data=hull_df , colour = "black" , fill = "blue" , alpha = .1)+
  geom_point(data=obs_lep  , aes(lon,lat) , colour = "red" , size = .5)+
  geom_point(data=catObs, aes(lon,lat) , size = 2)+
  facet_wrap(~species , nrow = 4 , ncol = 5 , labeller = label_wrap_gen(width=10))+
  theme_linedraw()+
  scale_x_continuous(breaks = c(14, 18, 22))+
  theme(strip.text = element_text(size = 10) , 
        axis.text = element_text(size = 10))+
  labs(y = "Latitude" , x = "Longitude")

return(p)
}


# Filter and plot species in each category ----------------------------------------------------
# New species
p1 <- plot_fun(category = "NEW TO SWEDEN")
# Range expansions
p2 <- plot_fun(category = "RANGE EXPANSION")

# Putative new Species clusters
catObs   <- filter(IBA_lep_occ ,CATEGORY == "PUTATIVE NEW SPECIES")  |> drop_na(lat)

p3 <- ggplot()+
  geom_sf(data=sweden)+
  geom_point(data=catObs, aes(lon,lat) , size = 2)+
  facet_wrap(~cluster , nrow = 4 , ncol = 10 , labeller = label_wrap_gen(width=10))+
  theme_linedraw()+
  scale_x_continuous(breaks = c(14, 18, 22))+
  theme(strip.text = element_text(size = 10) , 
        axis.text = element_text(size = 10))+
  labs(y = "Latitude" , x = "Longitude")


# save ----------------------------------------------------------------------------------------

ggsave(filename = "figs/new_to_sweden.jpeg", device="jpeg", width = 8, height = 8, p1)
ggsave(filename = "figs/range_expansions.jpeg", device="jpeg", width = 8, height = 8, p2)
ggsave(filename = "figs/putative_new_species.jpeg", device="jpeg", width = 12, height = 8, p3)


browseURL("figs/new_to_sweden.jpeg")
browseURL("figs/range_expansions.jpeg")
browseURL("figs/putative_new_species.jpeg")


