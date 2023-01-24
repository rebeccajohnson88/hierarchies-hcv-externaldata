


## Packages
library(tigris)
library(ggplot2)
library(rgdal)
library(sf)
library(tigris)
library(tidyverse)
library(tictoc)


## 
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)

## function
## then, for each pha, get tracts that are fully contained within
contained_within = function(one_pha, tracts){
  
  ## logical flag for one pha
  tracts_inPHA_logical = map_lgl(st_within(tracts, one_pha),
                                 function(x) {
                                   if(length(x) == 1){
                                     return(TRUE)
                                   } else {
                                     return(FALSE)
                                   }
                                 })
  print(tracts_inPHA_logical)
  geoids_within = as.character(paste(tracts$GEOID[tracts_inPHA_logical],
                                     collapse = "; "))
  return(cbind.data.frame(PARTICIPAN = one_pha$PARTICIPAN,
                          geoids_within))
  
}


## read in data
tract_df = readRDS("../Data/Raw/spatial_derivedobjs/tracts_foroverlap.RDS")
pha_df = readRDS("../Data/Raw/spatial_derivedobjs/phas_foroverlap.RDS")
print("read in files")

all_states = as.character(unique(pha_df$STUSAB))

tic("overlap loop")

for(i in 1:length(all_states)){
  
  ## filter tract and pha to same stata
  tracts_onestate = tract_df %>%
    filter(state_abbrev == all_states[i])
  phas_onestate = pha_df %>% filter(STUSAB == all_states[i])
  print(sprintf("Estimating for %s, which has %s tracts and %s PHAS",
                all_states[i],
                nrow(tracts_onestate),
                nrow(phas_onestate)))
  
  ## first, get intersection
  tracts_thatintersect = st_intersection(tracts_onestate,
                                         phas_onestate) # returns info on both but only pha geometry
  print("Found intersection")
  
  ## try merging with one
  ## need to get indicator that it's within a tract
  #all_contain = list()
  #for(i in 1:length(phas_onestate)){
    
   # one_contain = contained_within(phas_onestate[i, ], tracts_onestate)
  #  all_contain[[i]] = one_contain
    
  #}
  
  ## bindi nto one df
  #all_contain_df  = do.call(rbind.data.frame, all_contain)
  #print("Found tracts contained within")
  
  ## merge 
  #tract_intersect_contain = merge(tracts_thatintersect,
                           #       all_contain_df,
                            #      by = "PARTICIPAN",
                             #     all.x = TRUE)
  #print("Merged")
  
  #tract_intersect_contain = tract_intersect_contain  %>%
   # mutate(geoid_char = as.character(geoids_within),
       #    tract_contained_within = unlist(lapply(GEOID, 
        #                                          function(x) ifelse(any(grepl(x, geoids_within)), 1, 0))))
  
 # print("created indicator var")
  
  save_directory = "/Users/raj2/Dropbox/EligibilityPaper/PreferenceInference/Data/Intermediate/PHA_tract_bystate/"
  filename = sprintf("%s_intersects.RDS",
                     all_states[i])
  saveRDS(tracts_thatintersect, 
          sprintf("%s%s", save_directory, filename))
  print("saved")
  
}

toc()
