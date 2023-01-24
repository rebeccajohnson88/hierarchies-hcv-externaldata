


#### packages
rm(list = ls())
library(tigris)
library(ggplot2)
library(sf)
library(tigris)
library(tidyverse)
library(gridExtra)
library(cdlTools)
library(stringr)
library(tictoc)
library(tidycensus)
library(reshape2)
library(readxl)
library(here)
library(data.table)
library(here)
library(yaml)


#### constants
PATH_DATA = "EligibilityPaper/PreferenceInference/Data"


options(tigris_class = "sf")
options(tigris_use_cache = TRUE)


####  Read in tract PHA intersections to get set of tracts to pull for

get_geoid <- function(x, path = here("PreferenceInference/Data/Intermediate/PHA_tract_bystate/")){
  
  data_loaded = readRDS(sprintf("%s%s", path, x))
  geoid = as.character(data_loaded$GEOID)
  return(geoid)
}


all_tracts_df = list.files(path = here("PreferenceInference/Data/Intermediate/PHA_tract_bystate/"),
                           pattern = "intersects")
tracts_df_list = lapply(all_tracts_df, get_geoid)
tracts_df_list_df = lapply(tracts_df_list, function(x) data.frame(GEOID = x))
state_names = gsub("\\_intersects\\.RDS",
                   "",
                   all_tracts_df)
names(tracts_df_list_df) = state_names 
tracts_df = do.call(rbind.data.frame, tracts_df_list_df)
tracts_df = tracts_df %>%
  mutate(state = gsub("\\..*", 
                      "", 
                      rownames(tracts_df)))

####  Pull ACS vars
acs_vars = load_variables(2019,
                          "acs5", cache = TRUE)

total_pop = c("B00001_001")
race = c("B02001_001",
         "B02001_002",
         "B02001_003", 
         "B02001_004", 
         "B02001_005",
         "B02001_006",
         "B02001_007",
         "B02001_008",
         "B03001_001",
         "B03001_002",
         "B03001_003",
         sprintf("B03002_00%s", seq(1, 9)),
         sprintf("B03002_01%s", seq(0, 9)),
         sprintf("B03002_02%s", seq(0, 1)))

educ_vars = c("B16010_001", "B16010_002", "B16010_015", # this skips numbers due to labor force breakdowns within each category
              "B16010_028", "B16010_041", "B19013_001")
housing_status = c("B25064_001",
                   "B25123_002",
                   "B25123_008",
                   "B25123_001",
                   "B07013_003", # renter occupied
                   "B07013_001",
                   "B07013_002", #
                   "B25002_001", # vacancy 
                   "B25002_003", 
                   "B25071_001",
                   "B25003_001",
                   "B25003_003") # rental burden 

poverty = c("B06012_001",
            "B06012_002",
            "B06012_003",
            "B06012_004",
            "B09010_001", # snap
            "B09010_002",
            "C18120_001", # total measured
            "C18120_002", # total in labor force
            "C18120_006") # unemployment  



veteran = c("B21001_002",
            "B21001_001")

vars_topull = c(total_pop, 
                race,
                educ_vars, housing_status,
                poverty, veteran)

stopifnot(length(unique(vars_topull)) == length(vars_topull))

View(acs_vars %>% filter(name %in% vars_topull))


creds = read_yaml(here("PreferenceInference/my_cred_testingV2.yaml"))
census_api_key(creds$census_api$api_key)


us <- unique(fips_codes$state)[1:51]
pull_all_tracts = lapply(us,  function(x) get_acs(geography = "tract",
                                                  variables = vars_topull,
                                                  state = x,
                                                  year = 2016))


## save copy 
saveRDS(pull_all_tracts,
        here("PreferenceInference/Data/Raw/spatial_derivedobjs/all_tracts_longernames_20221216.RDS"))

## if reading in
READ_TRACTS = FALSE
if(READ_TRACTS){
  pull_all_tracts = readRDS(here("PreferenceInference/Data/Raw/spatial_derivedobjs/all_tracts_longernames_20221216.RDS"))
}

pull_all_tracts_df = do.call(rbind.data.frame,
                             pull_all_tracts)


####  RENAME ACS VARS
vars_pulled = rbind.data.frame(acs_vars %>%
  filter(name %in% vars_topull) %>%
  mutate(cl1 = gsub("Estimate!!Total:(!!)?",
                                  "", label),
         cl2 = ifelse(cl1 == "", name,
                      sprintf("%s_%s",
                      name, gsub("\\s+", "_", cl1))),
         cl3 = gsub("\\:|\\!|\\'", "", cl2)) %>%
  select(name, cl3),
  data.frame(name = "B00001_001",
             cl3 = "total_population_unweighted"))

## make sure no varnames are duplicated
stopifnot(nrow(vars_pulled %>% group_by(cl3) %>% filter(n() > 1)) == 0)
stopifnot(setdiff(unique(pull_all_tracts_df$variable), unique(vars_pulled$name)) == 0)

## merge tract-level chars w/ names
pull_tracts_wnames = merge(pull_all_tracts_df,
                           vars_pulled,
                           by.x = "variable",
                           by.y = "name",
                           all.x = TRUE)  %>%
  dplyr::select(GEOID, estimate, cl3) %>%
  reshape2::dcast(GEOID ~ cl3, value.var = "estimate") 

colnames(pull_tracts_wnames) = sprintf("%s_acscount",
                                       colnames(pull_tracts_wnames))

fwrite(pull_tracts_wnames, here("PreferenceInference/Data/Intermediate/acs_20122016_tractcounts_longernames_20221216.csv"))

####  Merge w/ pha code data
get_geoid_pcode <- function(x, path =  here("PreferenceInference/Data/Intermediate/PHA_tract_bystate/")){
  
  data_loaded = readRDS(sprintf("%s%s", path, x))
  data_return = as.data.frame(data_loaded %>% dplyr::select(PARTICIPAN, 
                                                            GEOID)) %>%
    dplyr::select(-geometry)
  return(data_return)
}


tracts_df_list_wphacode = lapply(all_tracts_df, get_geoid_pcode)
tracts_wphacode = do.call(rbind.data.frame, tracts_df_list_wphacode)

## left join census dems
tracts_overlap = unique(intersect(tracts_wphacode$GEOID,
                                  pull_tracts_wnames$GEOID_acscount))

print(sprintf("%s of the %s tracts pulled are in the PHA data",
              length(tracts_overlap),
              length(unique(pull_tracts_wnames$GEOID_acscount))))


## left join onto pha data
tracts_wdem = merge(tracts_wphacode,
                    pull_tracts_wnames,
                    by.x = "GEOID",
                    by.y = "GEOID_acscount",
                    all.x = TRUE)

## store in spatial obs
saveRDS(tracts_wdem,
        here("PreferenceInference/Data/Intermediate/phas_wrawACScounts_longernames_20221216.RDS"))



