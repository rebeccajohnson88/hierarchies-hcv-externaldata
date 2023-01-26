# 06_join_covariates.R
# Joins all PHA level covariates, including
# 1) PHA characteristics from HUD
# 2) preference data from Abt Associates 
# 3) variables from external sources created in script 05

# Source others
library(here)
source(here("src/utils.R"))

###########################################################################
## Load different input data sources
###########################################################################

# Load coded preferences
codes <- readRDS(here("data/raw/reviewed_codes.RDS"))

# Load data on source of preference information 
codes_sources <- read_csv(here("data/raw/all_codes_source_info.csv"))

# Load HUD Voucher Management System reports 
vms_fnames <- list.files(here("data/raw/"), pattern = "^VMS_")

# Read in VMS data
read_vms <- function(fn) {
  print(fn)
  vms_report <- read_excel(paste0("data/raw/", fn), guess_max = 5000) %>%
    transmute(pha_code = `PHA Code`, pha_name = `PHA Name`, month = Month,
              total_vouchers = `Total Vouchers`) %>%
    filter(!is.na(pha_code))
  
  return(vms_report)
}

vms_reports <- bind_rows(lapply(vms_fnames, read_vms)) 

###########################################################################
## Clean the voucher management system reports - used for ##
## this variable hud_med_month_vouchers                   ##
###########################################################################

# Summarize max and min number of vouchers administered
# by each PHA monthly during this period 
vms_summary <- vms_reports %>%
  group_by(pha_code) %>%
  summarize(num_months = n(),
            min_monthly_vouchers = min(total_vouchers),
            med_month_vouchers = median(total_vouchers),
            max_monthly_vouchers = max (total_vouchers)) 

# Set universe to only PHAs with that consistently administered
# vouchers over period from April 2016-March 2018
pha_pop <- vms_summary %>%
  # Filter to only observations that administered vouchers each month
  filter(min_monthly_vouchers > 0 & num_months == 24) %>% 
  # Filter out Puerto Rico, Guam, Virgin Islands
  filter(!(str_sub(pha_code, 1, 2) %in% c("RQ", "GQ", "VQ", "TQ"))) 

# Filter out any preferences for housing authorities outside of universe
# (Results in removal of 5: 
# KS169, ND057, ND028, MN017, FL083)
codes_filtered <- codes %>%
  semi_join(pha_pop %>% dplyr::select(pha_code, 
                                      med_month_vouchers),
            by = "pha_code") 

############################################
## Load PHA characteristics from HUD data ##
############################################

# Load data from HUD about housing authority characteristics
pha_basedf <- st_read(here("data/raw/a00000009.gdbtable")) %>%
  as.data.frame() %>%
  select(-Shape) %>%
  mutate(pha_code = as.character(PARTICIPANT_CODE))

# Join HUD variables to universe
# note to remove: deleted a lot from this; if missing vars later, add back
pha_hud <- pha_pop %>%
  ungroup() %>%
  left_join(pha_basedf, by = "pha_code") %>%
  rename_all(.funs = tolower) %>%
  mutate_at(vars(starts_with("annl_expns_amnt")), funs(ifelse(. < 0, NA, .))) %>%
  transmute(pha_code,
            pha_name = formal_participant_name,
            hud_med_month_vouchers = med_month_vouchers, # Median monthly vouchers based on VMS reports
            hud_count_all_units = pha_total_units,
            hud_geo_fips_state = state2kx) 


#################################
## Load external PHA context data  ##
#################################

pha_spatial <- readRDS(here("data/intermediate/pha_wlocalattributes_final.RDS")) %>% 
  rename(pha_code = HUD_PARTICIPANT_CODE) %>%
  select(-starts_with("HUD_"))  

names(pha_spatial) <- make.names(names(pha_spatial))
names(pha_spatial) <- str_remove_all(names(pha_spatial), "\\.")

pha_spatial_tojoin <- pha_spatial %>%   
  mutate(urban_state_flag = ifelse(is.na(URBANINST_hai20102014_state_flag), 1, URBANINST_hai20102014_state_flag),
         # Affordable per 100 ELI renters, including subsidized units
         urban_affordable_per100 = ifelse(urban_state_flag == 1, 
                                          URBANINST_hai20102014_ST_per100, URBANINST_hai20102014_per100),
         # Naturally affordable per 100 ELI renters 
         urban_affordable_noassist_per100 = ifelse(urban_state_flag == 1, 
                                                   URBANINST_hai20102014_ST_per100_No_Assisted, URBANINST_hai20102014_per100_no_assisted)) %>%
  select(-starts_with("URBANINST")) %>%
  rename_all(tolower)

# Load Census region reference
census_regions <- read_excel(here("data/raw/state-geocodes-v2011.xls"), 
                             na = c("", "-99"),
                             skip = 5,
                             guess_max = 5000) %>%
  transmute(state_code = `State\n(FIPS)`, 
            census_region = case_when(Region == 1 ~ "Northeast",
                                      Region == 2 ~ "Midwest",
                                      Region == 3 ~ "South",
                                      Region == 4 ~ "West"))
  
# Join to main data frame
pha_all <- pha_hud %>%
  left_join(pha_spatial_tojoin, by = "pha_code") %>%
  left_join(census_regions, by = c("hud_geo_fips_state" = "state_code"))

# Filtering
pha_all_filtered <- pha_all %>% dplyr::select(pha_code, pha_name,
            hud_med_month_vouchers, 
            hud_count_all_units,
            derived_acs_not_hispanic_or_latinowhite_alone_percent, 
            derived_acs_not_hispanic_or_latinoblack_or_african_american_alone_percent, 
            derived_acs_not_hispanic_or_latinoasian_alone_percent,
            derived_acs_hisp_any_perc, derived_acs_veteran_percent, 
            derived_acs_below_100_percent_of_the_poverty_level_percent, 
            derived_acs_in_the_labor_forceunemployed_percent,
            derived_acs_living_in_household_with_supplemental_security_income_ssi_cash_public_assistance_income_or_food_stampssnap_in_the_past_12_months_percent,
            derived_median_hh_income, derived_median_rental_burden,
            derived_acs_renter_occupied_percent, derived_acs_vacant_percent,
            perc_tracts_urban, 
            urban_affordable_noassist_per100, 
            urban_affordable_per100, 
            sharkey_commorg_densper_1000_inpov, 
            human_services_densper_1000_inpov, 
            repub_2016, census_region, coc_category_noimpute)


## write to both repo-specific directory
## and to the replication archive 
### repo-specific cleaned folder
saveRDS(pha_all_filtered, here("data/cleaned/pha_all.RDS"))
### replication archive 
saveRDS(pha_all_filtered, here("../EligibilityPaper/PreferenceInference/hcv_replication/data/pha_all.RDS"))
fwrite(pha_all_filtered, here("../EligibilityPaper/PreferenceInference/hcv_replication/data/pha_all.csv"))

