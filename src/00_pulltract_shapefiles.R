

# 00_pulltract_shapefiles.R
# Pull shapefiles for the census tracts 

library(here)
source(here("src/utils.R"))
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)

READ_RAW_SHAPES <- TRUE
WRITE_NEW_INTERMEDIATE_FILES <- FALSE 

###########################################################################
## Load census tract shapefiles - created by 0helper script
###########################################################################

# Load census tracts to get the projection to 
# convert the PHA shapefile to 
# these are pulled by the 01_
tracts_shp <- readRDS(here("data/intermediate/tracts_sf_format.RDS"))

# Make sure converted to sf format
tracts_sf_format <- st_as_sf(tracts_shp)


###########################################################################
## Read in estimated housing service authority service areas
## Dataset from HUD described in supplementary materials
###########################################################################

if(READ_RAW_SHAPES){
  pha_shpformat <- readOGR(here("data/raw/Estimated_Housing_Authority_Service_Areas/"))
  pha_sf_format <- st_as_sf(pha_shpformat)
  ## align the projections to census tract projections
  st_crs(pha_sf_format) <- st_crs(tracts_sf_format)
  pha_sf_format <- st_transform(pha_sf_format,
                               crs = st_crs(tracts_sf_format))
} else{
  pha_sf_format <- readRDS(here("data/intermediate/PHA_sf_format.RDS"))
} 


###########################################################################
## Merge tract and PHA service area by state and fips codes
###########################################################################

# Convert pha state abbrevs to FIPS codes
fips_ids <- data.frame(state_abbrev = unique(pha_sf_format$STUSAB),
                      fips = fips(unique(pha_sf_format$STUSAB),
                                  to= "FIPS")) %>%
                    mutate(fips = ifelse(nchar(fips) < 2,
                                sprintf("0%s", fips),
                                fips))

# Left join to the ACS tracts data
tracts_sf_format_wstates = merge(tracts_sf_format,
                                 fips_ids,
                                 by.x = "STATEFP",
                                 by.y = "fips",
                                 all.x = TRUE)

# Add state info back in 
pha_sf_format_wstates = merge(pha_sf_format,
                              fips_ids,
                              by.x = "STUSAB",
                              by.y = "state_abbrev",
                              all.x = TRUE)


###########################################################################
## Save records of tracts and PHAs that might overlap because in same state
###########################################################################

if(WRITE_NEW_INTERMEDIATE_FILES){
  saveRDS(tracts_sf_format_wstates, here("data/intermediate/tracts_foroverlap.RDS"))
  saveRDS(pha_sf_format_wstates, here("data/intermediate/phas_foroverlap.RDS"))
}

