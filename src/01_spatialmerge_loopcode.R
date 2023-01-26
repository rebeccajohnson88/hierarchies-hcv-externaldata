# 01_spatialmerge_loopcode.R
# Within each state, merge the PHA service area tracts 
# with the tract shapefiles


library(here)
source(here("src/utils.R"))
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)


###########################################################################
## Load census tract shapefiles and PHA shapefiles
###########################################################################
tract_df <- readRDS(here("data/intermediate/tracts_foroverlap.RDS"))
pha_df <- readRDS(here("data/intermediate/phas_foroverlap.RDS"))

print("read in files")

all_states = as.character(unique(pha_df$STUSAB))

tic("overlap loop")


###########################################################################
## Iterate by state and find intersection b/t PHA and tracts in that
## state (separate by state due to file size)
###########################################################################
for(i in 1:length(all_states)){
  
  ## filter tract and pha to same state
  tracts_onestate = tract_df %>%
    filter(state_abbrev == all_states[i])
  phas_onestate = pha_df %>% filter(STUSAB == all_states[i])
  print(sprintf("Estimating for %s, which has %s tracts and %s PHAS",
                all_states[i],
                nrow(tracts_onestate),
                nrow(phas_onestate)))
  
  ## get spatial intersection
  tracts_thatintersect = st_intersection(tracts_onestate,
                                         phas_onestate) # returns info on both but only pha geometry
  print("Found intersection")
  
  # specify directory to save result 
  save_directory = here("data/intermediate/PHA_tract_bystate/")
  filename = sprintf("%s_intersects.RDS",
                     all_states[i])
  saveRDS(tracts_thatintersect, 
          sprintf("%s%s", save_directory, filename))
  print("saved")
  
}

toc()
