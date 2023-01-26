
# 02_spatialmerge_CoC.R
# Within each state, find intersection between CoC shapefiles
# and PHA shapefile 
library(here)
source(here("src/utils.R"))
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)

#######################################
# Read in data
#######################################

## read in PHA estimated service areas file
pha_shapes = readRDS(here("data/intermediate/phas_foroverlap.RDS"))
sprintf("There are %s unique PHAs in the estimated service area file",
        length(unique(pha_shapes$PARTICIPAN)))

## read in PHA point locations 
phas_gdp <- sf::st_read(dsn = here("data/raw/PHA_gdb/PHAs.gdb"))
sprintf("There are %s PHAs in the universe",
        length(unique(phas_gdp$PARTICIPANT_CODE)))

## read in continuum of care shapefiles
coc_shapefiles <- sf::st_read(dsn = here("data/raw/Continuum_of_Care_(CoC)_Grantee_Areas/"))
sprintf("There are %s coc", length(unique(coc_shapefiles$OBJECTID)))

#######################################
# Find intersection of PHA points and CoC shapes
# (PHA shapes and CoC shapes produces similar results)
#######################################

## convert coc coordinate projection to match phas
coc_coordchange <- st_buffer(st_transform(coc_shapefiles, st_crs(phas_gdp)), dist = 0)
                
## get centroids of cocs for labels
coc_centroids <- cbind.data.frame(coc_coordchange, st_coordinates(st_centroid(coc_coordchange)))

## uncomment to test with Illinois
# IL_phas <- phas_gdp %>% filter(grepl("^IL", PARTICIPANT_CODE))
# IL_coc <- coc_coordchange %>% filter(STUSAB == "IL")
# 
# ## get intersect
# IL_pointintersect <- st_intersection(IL_coc, IL_phas)
# 
# focal_pha <- "IL107"
# intersect_coc <- IL_pointintersect %>%
#   filter(PARTICIPANT_CODE == "IL107") %>%
#   pull(OBJECTID)
# 
# ## plot example
# ggplot(IL_coc) +
#   geom_sf() +
#   geom_sf(data = IL_phas %>%
#         filter(PARTICIPANT_CODE == "IL107"), size = 2, shape = 23, fill = "darkred") +
#   geom_text(data = coc_centroids %>%
#               filter(OBJECTID == intersect_coc), 
#   aes(x = X, y = Y, label = COCNAME), size = 2) 

## repeat for all; not constraining by state but should be state matched
all_pointintersect <- st_intersection(coc_coordchange, phas_gdp)
phas_dropped_intersect <- setdiff(phas_gdp$PARTICIPANT_CODE,
                                  all_pointintersect$PARTICIPANT_CODE)

#######################################
# Write results
#######################################

saveRDS(all_pointintersect, 
        here("data/intermediate/phapoint_coc_intersect.RDS"),
        compress = FALSE)
