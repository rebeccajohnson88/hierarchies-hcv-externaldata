# 05_add_all_attributes_tobaseHUD_df
# Merge and clean different PHA-level characteristics
# Outputs: pha_wlocalattributes_final.[csv|RDS]

# Constants
GRAPH_DIR <- here("output/")
RAW_DIR <- here("data/raw/")
INTERMEDIATE_DIR <- here("data/intermediate/") 
PULL_NONPROFIT <- FALSE
WRITE_NEW_INTERMEDIATE_FILES <- FALSE

# Source others
source(here("src/utils.R"))
source(here("src/helperfunc_nonprofitdata.R"))

###########################################################################
## Load different input data sources
###########################################################################

# PHA point location data that contains universe of PHAs
phas_gdp <- sf::st_read(dsn = 
        here(RAW_DIR, "PHA_gdb/PHAs.gdb"))
colnames(phas_gdp) = sprintf("HUD_%s",
                             colnames(phas_gdp))


# Raw ACS counts of people in diff categories
phas_rawACScounts = readRDS(here(INTERMEDIATE_DIR, "phas_wrawACScounts_longernames_20221216.RDS"))


###########################################################################
## Clean ACS demographic data
###########################################################################

## first, aggregate demographics upwards from tract to PHA and convert counts
## to percentages where relevant. Agg rules:
## 1. For categorical variables, sum the counts of category members in each tract
## that intersects with the PHA. Eg if PHA 1
## intersects with Tract A that has 50 veterans and Tract B that has 150 veterans,
## sum is 100
## 2. For continuous (eg median household income), take the mean across tracts
acs_countbyPHA = phas_rawACScounts %>%
  group_by(PARTICIPAN) %>%
  summarise_at(vars(-GEOID, -`B25064_001_EstimateMedian_gross_rent_acscount`,
                    -`B25071_001_EstimateMedian_gross_rent_as_a_percentage_of_household_income_acscount`,
                    -`B19013_001_EstimateMedian_household_income_in_the_past_12_months_(in_2019_inflation-adjusted_dollars)_acscount`),
               sum) %>%
  left_join(phas_rawACScounts %>%
              group_by(PARTICIPAN) %>%
              summarise(mean_median_hhincome = mean(`B19013_001_EstimateMedian_household_income_in_the_past_12_months_(in_2019_inflation-adjusted_dollars)_acscount`,
                                                    na.rm = TRUE),
                        mean_median_grossrentpercincome = mean(`B25071_001_EstimateMedian_gross_rent_as_a_percentage_of_household_income_acscount`,
                                                               na.rm = TRUE)),
            by = 'PARTICIPAN')

## now that found sum, calculate the percentages for categorical variables
## from the numerator and denominator counts
acs_countbyPHA <- acs_countbyPHA %>%
  mutate(derived_acs_not_hispanic_or_latinowhite_alone_percent = B03002_003_Not_Hispanic_or_LatinoWhite_alone_acscount/B03002_001_acscount,
         derived_acs_not_hispanic_or_latinoblack_or_african_american_alone_percent = B03002_004_Not_Hispanic_or_LatinoBlack_or_African_American_alone_acscount/B03002_001_acscount,
         derived_acs_not_hispanic_or_latinoasian_alone_percent = B03002_006_Not_Hispanic_or_LatinoAsian_alone_acscount/B03002_001_acscount,
         derived_acs_hisp_any_perc = B03002_012_Hispanic_or_Latino_acscount/B03002_001_acscount,
         derived_acs_veteran_percent = B21001_002_Veteran_acscount/B21001_001_acscount,
         derived_acs_living_in_household_with_supplemental_security_income_ssi_cash_public_assistance_income_or_food_stampssnap_in_the_past_12_months_percent = 
           `B09010_002_Living_in_household_with_Supplemental_Security_Income_(SSI),_cash_public_assistance_income,_or_Food_Stamps/SNAP_in_the_past_12_months_acscount`/B09010_001_acscount,
         derived_acs_below_100_percent_of_the_poverty_level_percent = B06012_002_Below_100_percent_of_the_poverty_level_acscount/B06012_001_acscount,
         derived_acs_in_the_labor_forceunemployed_percent = C18120_006_In_the_labor_forceUnemployed_acscount/C18120_002_In_the_labor_force_acscount,
         derived_acs_total_unemployed_percent = C18120_006_In_the_labor_forceUnemployed_acscount/C18120_001_acscount,
         derived_acs_renter_occupied_percent = B25003_003_Renter_occupied_acscount/B25003_001_acscount,
         derived_median_hh_income = mean_median_hhincome,
         derived_median_rental_burden = mean_median_grossrentpercincome,
         derived_acs_vacant_percent = B25002_003_Vacant_acscount/B25002_001_acscount)


###########################################################################
## Clean and merge continuum of care data 
###########################################################################

coc_points <- readRDS(here(INTERMEDIATE_DIR, "phapoint_coc_intersect.RDS"))

coc_points_relcol <- coc_points %>%
                select(COCNAME, PARTICIPANT_CODE) %>%
                st_drop_geometry() %>%
                mutate(is_statewide_coc = case_when(grepl("Balance of State|Balance of Commonwealth",
                              COCNAME) ~ TRUE, 
                              TRUE ~ FALSE)) 

## get prefixes of statewide coc
has_statewide_coc <- unique(coc_points_relcol %>%
            filter(is_statewide_coc) %>%
            mutate(state_twodig = str_extract(PARTICIPANT_CODE, 
                                              "[A-Z][A-Z]")) %>%
            pull(state_twodig))

## left join onto main df and create binary indicator
acs_wcoc <- merge(acs_countbyPHA,
                                   coc_points_relcol,
                                   by.x = "PARTICIPAN",
                                   by.y = "PARTICIPANT_CODE",
                                   all.x = TRUE) %>%
                        mutate(COC_category_noimpute = 
                              case_when(is_statewide_coc ~ "Statewide COC",
                                     !is_statewide_coc & !is.na(COCNAME) ~ "Local COC",
                                     TRUE ~ "No COC"),
                              COC_category_impute = 
                              case_when(is_statewide_coc ~ "Statewide COC",
                                     !is_statewide_coc & !is.na(COCNAME) ~ "Local COC",
                                     is.na(COCNAME) & str_extract(as.character(PARTICIPAN),
                                    "[A-Z][A-Z]") %in% is_statewide_coc ~ "Statewide COC",
                                     TRUE ~ "No COC"),
                              COC_imputed = ifelse(COC_category_noimpute != COC_category_impute, TRUE, FALSE)) 


###########################################################################
## Clean Urban institute housing affordability data to prepare to merge
###########################################################################

# Load data
hai_vars <- read.csv(here(RAW_DIR,
                'HAI_map_201014.csv'))

# Add suffix to colnames
hai_vars_addsuffix <- sprintf("URBANINST_hai20102014_%s",
                  colnames(hai_vars))

colnames(hai_vars) <- hai_vars_addsuffix

# these data are merged later in script


###########################################################################
## Clean county-level election results to prepare to merge
###########################################################################

election_raw <- read.csv(here(RAW_DIR, "countypres_2000-2016.csv"))

## subset to 2016 and merge onto main data
election_raw_bycounty  <- election_raw %>%
          filter(year %in% c(2012, 2016)) %>%
          dplyr::select(state, FIPS, candidatevotes, 
                        totalvotes, 
                        party, year) %>%
          filter(party %in% c("democrat", 
                              "republican")) %>%
          mutate(party_year = sprintf("%s_%s",
                  party, year)) %>%
          dplyr::select(-party, -year)

election_bycounty_wide_init <- reshape(election_raw_bycounty,
                                 idvar = c("state", "FIPS"),
                                 timevar = c("party_year"), 
                                 direction = "wide", 
                                 sep = "_") 

## get percs
election_recent <- election_bycounty_wide_init %>%
            mutate(dem_2012 = candidatevotes_democrat_2012/totalvotes_democrat_2012,
                   dem_2016 = candidatevotes_democrat_2016/totalvotes_democrat_2016,
                   repub_2012 = 
candidatevotes_republican_2012/totalvotes_republican_2012,
                repub_2016 = 
candidatevotes_republican_2016/totalvotes_republican_2016) %>%
          dplyr::select(-contains("votes"))

## get state level aggregation
election_recent_statelevel <- election_bycounty_wide_init %>%
          group_by(state) %>%
          dplyr::select(-FIPS) %>% 
          summarise_if(is.numeric, mean, na.rm = TRUE) %>%
          mutate(dem_2012 = candidatevotes_democrat_2012/totalvotes_democrat_2012,
                   dem_2016 = candidatevotes_democrat_2016/totalvotes_democrat_2016,
                   repub_2012 = 
candidatevotes_republican_2012/totalvotes_republican_2012,
                repub_2016 = 
candidatevotes_republican_2016/totalvotes_republican_2016) %>%
          dplyr::select(-contains("votes")) %>%
          left_join(election_raw %>%
                dplyr::select(state, state_po) %>%
                filter(!duplicated(state)),
                by = "state")
          

###########################################################################
## Clean rurality data to merge
###########################################################################

# Read in data
ruca <- read_xlsx(here(RAW_DIR, "ruca2010revised.xlsx"))

# First row is colname and second row is data; subset
ruca_df <- ruca %>%
  slice(2:nrow(ruca))
ruca_colnames <- ruca[1, ]

# Clean up colnames
colnames(ruca_df) <- gsub("\\_+$", "", gsub("\\s+|Code.*|\\-", "_", ruca_colnames))

# Prep for merge
ruca_df_tomerge = ruca_df %>%
              mutate(urban = ifelse(Primary_RUCA %in% seq(from = 1, to = 3),
                              1, 0),
                     rural = ifelse(Primary_RUCA  %in% seq(from = 4, to = 10),
                              1, 0)) %>%
              dplyr::select(State_County_Tract_FIPS,
                            Primary_RUCA, urban, rural)  

## aggregate to percentage of tracts 
# a pha overlaps with that are classified as each type
#  merge with pha by tract info
pha_by_tract <- phas_rawACScounts %>% dplyr::select(GEOID, PARTICIPAN) 
sprintf("Of the %s tracts that overlap with PHAs, %s have ruca codes",
        length(unique(pha_by_tract$GEOID)),
        length(intersect(unique(pha_by_tract$GEOID), unique(ruca_df_tomerge$State_County_Tract_FIPS))))
pha_by_tract_wruca <- merge(pha_by_tract, 
                           ruca_df_tomerge,
                           all.x = TRUE, 
                           by.x = "GEOID",
                           by.y = "State_County_Tract_FIPS")

## calculate percent- just do perc_tract_urban
pha_ruca_perc <- pha_by_tract_wruca %>%
            group_by(PARTICIPAN) %>%
            summarise(perc_tracts_urban = mean(urban, na.rm = TRUE)) 
            


###########################################################################
## Clean nonprofit density data 
###########################################################################


## warning- getcorefile from Urban Institute function deprecated; see intermediate directory for the files
if(PULL_NONPROFIT){
  core2015pc = getcorefile(2015, "pc")
  core2015pf = getcorefile(2015, "pf")
  core2015co = getcorefile(2015, "co")
  saveRDS(core2015pc, here(INTERMEDIATE_DIR,"core2015pc.RDS"))
  saveRDS(core2015pf,  here(INTERMEDIATE_DIR,"core2015pf.RDS"))
  saveRDS(core2015pf,  here(INTERMEDIATE_DIR,"core2015co.RDS"))
} else{
  core2015pc <- readRDS(here(INTERMEDIATE_DIR,"core2015pc.RDS"))
  core2015pf <- readRDS(here(INTERMEDIATE_DIR,"core2015pf.RDS"))
  core2015co <- readRDS(here(INTERMEDIATE_DIR,"core2015co.RDS"))
}

# Store three datasets in list
all_nccs <- list(pc = core2015pc, pf = core2015pf, 
                other = core2015co)

# Iterate and clean the datasets
cleaned_nccs <- list()
for(i in 1:length(all_nccs)){
  
  which_type <- names(all_nccs)[1]
  cleaned_nccs[[i]] <- clean_select_nccs(all_nccs[[i]], which_type = which_type)

}


# Rbind back together
cleaned_nccs_all <- do.call(rbind.data.frame, cleaned_nccs) 


# See that there are more rows than unique eins, so filter out
# duplicates for the three columns we're using (tract; ein; nteecc code)
sprintf("There are %s rows in the NCCS data, and %s unique ids",
        nrow(cleaned_nccs_all),
        length(unique(cleaned_nccs_all$EIN)))
cleaned_nccs_dedup <- cleaned_nccs_all %>%
                distinct(EIN, CENSUSTRACT, NTEECC, .keep_all = TRUE)


# Create vectors for diff types of NTEE codes
sharkey_arts_notiv <- c("A23", "A25", "A27")
sharkey_crime <- c("I20", "I21", "F42", "I31", "I40",
                  "I43", "I44")
sharkey_neighbdev <- c("L25", "L30", "L80", 
                      "L81", "P28",
                      paste("S", 20:22, sep = ""),
                      paste("S", 30:31, sep = ""))

sharkey_substance <- paste("F", 20:22, sep = "")
sharkey_workforce <- c("J22", "J30")
sharkey_youthdev <- c("N60",  
                     paste("O", 20:23, sep = ""),
                     "O30", "O31", "O40",
                     paste("O", 50:55, sep = ""), 
                     "P27", "P30")
homelessness <- c("L41", "L99", "P85")

# Create dummy indicators
cleaned_nccs_dedup <- cleaned_nccs_dedup %>%
            mutate(ncee_first_letter = str_extract(NTEECC, "^[A-Z]"),
                   ncee_decile_code = str_extract(NTEECC, "^[A-Z][0-9][0-9]"),
                   sharkey_arts_iv = ifelse(ncee_first_letter %in% 
                  c("A", "C", "H") & 
                  !ncee_decile_code %in% sharkey_arts_notiv, 1, 0),
                  sharkey_crime = ifelse(ncee_decile_code %in% sharkey_crime, 1, 0),
                  sharkey_neighbdev = ifelse(ncee_decile_code %in% sharkey_neighbdev, 1, 0),
                  sharkey_substance = ifelse(ncee_decile_code %in% sharkey_substance, 1, 0),
                  sharkey_workforce = ifelse(ncee_decile_code %in% sharkey_workforce, 1, 0),
                  sharkey_youthdev = ifelse(ncee_decile_code %in% sharkey_youthdev, 1, 0),
                  sharkey_commorg = ifelse(sharkey_crime == 1 | 
                                    sharkey_neighbdev == 1 | 
                                    sharkey_substance == 1 | 
                                    sharkey_workforce == 1 | 
                                    sharkey_youthdev == 1, 1, 0),
                  human_services = ifelse(ncee_first_letter %in% 
                              c("I", "J", "K", "L", "M", "N", "O", "P"), 1, 0),
                  homelessness = ifelse(ncee_decile_code %in% homelessness, 1, 0))

# Aggregate from organization-tract dyad level to tract level to find
# organizations by tract
dummy_vars <- grep("sharkey|human|homelessness", colnames(cleaned_nccs_dedup), value = TRUE)
cleaned_nccs_forsummary <- cleaned_nccs_dedup[, c("EIN",
                                                 "CENSUSTRACT",
                                                 dummy_vars)] %>%
                          reshape2::melt(, id.vars = c("EIN", "CENSUSTRACT")) %>%
                          rename(type_nonprof = variable)
nonprofit_tract_summary <- cleaned_nccs_forsummary %>%
                  filter(CENSUSTRACT != "None" & !is.na(CENSUSTRACT)) %>%
                  group_by(CENSUSTRACT, type_nonprof) %>%
                  summarise(total_nonprof = sum(value)) %>%
                  ungroup() 
nonprofit_tract_summary_wide <- reshape2::dcast(nonprofit_tract_summary, CENSUSTRACT ~ type_nonprof, 
                                     value.var = "total_nonprof")

# Create dataframe with tracts in PHA data but that don't have any
# organizations in the nonprofit data so that these count as having 0 orgs
tracts_notin_nccs <- setdiff(unique(pha_by_tract$GEOID), 
                            nonprofit_tract_summary$CENSUSTRACT)


npsummary_cols <- setdiff(colnames(nonprofit_tract_summary_wide), "CENSUSTRACT")
notracts_cols <- data.frame(matrix(ncol = length(npsummary_cols),
                                  nrow = 1, rep(0, length(npsummary_cols))))
colnames(notracts_cols) <- npsummary_cols 
nonprofit_tract_summary_zero <- cbind.data.frame(data.frame(CENSUSTRACT = as.character(tracts_notin_nccs)),
                            notracts_cols)
nonprofit_tract_summary_tomerge <- rbind.data.frame(nonprofit_tract_summary_wide, 
                                                   nonprofit_tract_summary_zero)

# Merge with other information --- still at tract level
pha_by_tract_wpop <- phas_rawACScounts %>% dplyr::select(GEOID, PARTICIPAN) 
pha_wnonprofit_tractlevel <- merge(pha_by_tract_wpop,
                                  nonprofit_tract_summary_tomerge, 
                                  by.x = "GEOID",
                                  by.y = "CENSUSTRACT",
                                  all.x = TRUE)

# Aggregate from tract -> PHA level by finding average density per 1000 residents in poverty
nonprofit_vars <- grep("sharkey|human|homelessness", colnames(pha_wnonprofit_tractlevel), value = TRUE)
count_vars <- nonprofit_vars 

# For both measures, first find sum of organizations aggregating from tract -> PHA level
# this serves as the numerator
pha_wnonprofit_phalevel_sum <- pha_wnonprofit_tractlevel %>%
                      group_by(PARTICIPAN) %>%
                      summarise_at(all_of(count_vars),
                                   list(sum = sum)) %>%
                      ungroup() # statewide sums all tracts
colnames(pha_wnonprofit_phalevel_sum) = c("PARTICIPAN",
                                      sprintf("%s_sum", count_vars))

# Second- find denominator in terms of people <= 100% federal poverty line
# and merge it on
base_tomerge = as.data.frame(phas_gdp) %>%
  mutate(hud_char = as.character(HUD_COUNTY_LEVEL),
         hud_county_tomerge = ifelse(HUD_PARTICIPANT_CODE == "NJ114",
                                     "34023",
                                     hud_char)) # manually replace pha missing its county
pha_wnonprofit_denom <- merge(base_tomerge %>%
                                dplyr::select(HUD_PARTICIPANT_CODE, 
                                              HUD_PHA_TOTAL_UNITS),
                              pha_wnonprofit_phalevel_sum, 
                              by.x = "HUD_PARTICIPANT_CODE",
                              by.y  = "PARTICIPAN",
                              all.x = TRUE) %>%
                    left_join(acs_countbyPHA %>% dplyr::select(PARTICIPAN,
                                 B06012_002_Below_100_percent_of_the_poverty_level_acscount),
                          by = c('HUD_PARTICIPANT_CODE' = 'PARTICIPAN')) %>%
                    rename(count_100percfpl_below = B06012_002_Below_100_percent_of_the_poverty_level_acscount)

# Third, for each of the nonprofit variables, divide by that count of residents
# in poverty and multiply result by 1000 to get the rate per 1000 residents
nonprofit_vars_fordens <- grep("sum", colnames(pha_wnonprofit_phalevel_sum),
                              value = TRUE)
nonprofit_dens_bypov <- pha_wnonprofit_denom %>%
                    mutate_at(all_of(nonprofit_vars_fordens),
                              funs((./count_100percfpl_below)*1000))
colnames(nonprofit_dens_bypov) = gsub("\\_sum",
                                        "\\_densper_1000_inpov", colnames(nonprofit_dens_bypov))

# Fourth, merge back with pha-level data
pha_wnonprofit_wrates = pha_wnonprofit_phalevel_sum %>%
                        left_join(nonprofit_dens_bypov,
                          by = c("PARTICIPAN" = "HUD_PARTICIPANT_CODE")) %>%
                        mutate(state = str_extract(PARTICIPAN, 
                                                  "[A-Z][A-Z]")) 

# Fifth, restrict to variables needed to merge + used
# in analysis
pha_wnonprofit_tomerge <- pha_wnonprofit_wrates %>% dplyr::select(PARTICIPAN,
                          sharkey_commorg_densper_1000_inpov, 
                          human_services_densper_1000_inpov)

###########################################################################
## Merge these PHA context datasets
## Script 06 will merge those with the others 
###########################################################################

# Left join base of PHAs to ACS data that has
# coc information merged on 
pha_wdem <- merge(base_tomerge,
                 acs_wcoc,
                 by.x = "HUD_PARTICIPANT_CODE",
                 by.y = "PARTICIPAN",
                 all.x = TRUE)


# Left join RUCA
pha_wruca <- merge(pha_wdem, 
                  pha_ruca_perc,
                  by.x = "HUD_PARTICIPANT_CODE",
                  by.y = "PARTICIPAN",
                  all.x = TRUE)


# Left join nonprofit dens
pha_wnonprofit <- merge(pha_wruca, 
                       pha_wnonprofit_tomerge %>% rename(HUD_PARTICIPANT_CODE = PARTICIPAN),
                       by = "HUD_PARTICIPANT_CODE",
                       all.x = TRUE)

# Left join county-level housing affordability
pha_whai <- merge(pha_wnonprofit,
                 hai_vars,
                 by.x = "hud_county_tomerge",
                 by.y = "URBANINST_hai20102014_county",
                 all.x = TRUE)


# To prep for merging on elections data
# at either county or state level, merge on
# indicator for whether local or state PHA
pha_whai <- pha_whai %>%
              mutate(statelevel = grepl("^[A-Z][A-Z]9", HUD_PARTICIPANT_CODE))

pha_statewide <- pha_whai %>%
        filter(statelevel) %>%
        mutate(state_po = gsub("9[0-9][0-9]", "",
                               HUD_PARTICIPANT_CODE))

pha_county <- pha_whai %>%
      filter(!statelevel) 


# Merge on statewide results
pha_statewide_welect <- merge(pha_statewide,
                             election_recent_statelevel,
                             by = "state_po",
                             all.x = TRUE) %>%
                    dplyr::select(-state_po) #remove join var

# Check missingness, only guam, PR, and other territories
missing_2016 <- pha_statewide_welect %>%
            filter(is.na(dem_2016))
sprintf("Of the %s statewide PHAs, %s are missing 2016 election info, those are: %s",
        length(unique(pha_statewide_welect$HUD_PARTICIPANT_CODE)),
        length(unique(missing_2016$HUD_PARTICIPANT_CODE)),
        paste(missing_2016$HUD_FORMAL_PARTICIPANT_NAME, 
              collapse = ";"))

# Left join county-level results
pha_county_welect <- merge(pha_county,
          election_recent %>%
          mutate(hud_county_tomerge = str_pad(FIPS, 5,
                                side = "left", "0")) %>%
          dplyr::select(-FIPS),
                             by = "hud_county_tomerge",
                             all.x = TRUE)

missing_2016_county <- pha_county_welect %>%
            filter(is.na(dem_2016))

# See again that it's territories missing them
sprintf("Of the %s non-statewide PHAs, %s are missing 2016 election info, those are: %s",
        length(unique(pha_county_welect$HUD_PARTICIPANT_CODE)),
        length(unique(missing_2016_county$HUD_PARTICIPANT_CODE)),
        paste(missing_2016_county$HUD_FORMAL_PARTICIPANT_NAME, 
              collapse = ";"))


# Rowbind state and county
pha_welections <- rbind.data.frame(pha_statewide_welect,
                                  pha_county_welect)


###########################################################################
## Write the results to intermediate
###########################################################################

# Before writing, confirm no PHAs lost in joins
stopifnot(nrow(pha_welections) == nrow(phas_gdp))

# Write 
if(WRITE_NEW_INTERMEDIATE_FILES){
  saveRDS(pha_welections, 
          here(INTERMEDIATE_DIR, "pha_wlocalattributes_final.RDS"))
  fwrite(pha_welections,
         here(INTERMEDIATE_DIR, "pha_wlocalattributes_final.csv"))
}
