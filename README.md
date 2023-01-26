# Constructing PHA-level characteristics 🏠

## Description

For Simone Zhang and Rebecca A. Johnson, 2023, `Hierarchies in the Decentralized Welfare State: Prioritization in the Housing Choice Voucher Program,` <a href="https://doi.org/10.1177/00031224221147899" target="_blank">doi.org/10.1177/00031224221147899</a>

Code for acquiring and cleaning the external data sources used to describe PHAs and their local contexts. Analytic replication code is here: https://doi.org/10.7910/DVN/O3YYII

## Authors

- [Simone Zhang](https://simonezhang.com/) 
  - [https://github.com/sxzh](https://github.com/sxzh)
- [Rebecca A. Johnson](https://www.rebeccajohnson.io/) 
  - [https://github.com/rebeccajohnson88](https://github.com/rebeccajohnson88)

## Input data

Files reference in `data/raw/` are available in [this Dropbox folder](https://www.dropbox.com/sh/kd0r5w9li1pmtl1/AAB0-zU9Hv8awOa0uhOViHSva?dl=0). Other files can be generated from those raw files using the scripts.

## Code 

### Scripts to be run in order

- [00a_pulltract_shapefiles.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/00a_pulltract_shapefiles.R)

  - Takes in:
    - `data/raw/Estimated_Housing_Authority_Service_Areas`: shapefiles from HUD on PHA estimated service areas
  - What it does: uses the Census API/[tigris wrapper](https://www.rdocumentation.org/packages/tigris/versions/1.6.1/topics/tracts) to pull the tract polygons in each state in which PHAs are located 
  - Outputs:
    - `data/intermediate/tracts_sf_format.RDS`: tract polygons stored in [sf](https://cran.r-project.org/web/packages/sf/index.html) format

- [00b_pha_tract_merge.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/00b_pha_tract_merge.R)

  - Takes in:
    - `data/intermediate/tracts_sf_format.RDS`
    - `data/raw/Estimated_Housing_Authority_Service_Areas`
  - What it does: restricts tract shapefiles and PHA shapefiles to polygons in same states
  - Outputs:
    - `data/intermediate/tracts_foroverlap.RDS`
    - `data/intermediate/phas_foroverlap.RDS`

- [01_spatialmerge_loopcode.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/01_spatialmerge_loopcode.R)
  - Takes in:
    - `data/intermediate/tracts_foroverlap.RDS`
    - `data/intermediate/phas_foroverlap.RDS`
  - What it does: iterates over state and finds the spatial intersection between the tract polygon and the PHA service area polygon; writes each state's output separately due to file sizes 
  - Outputs:
    - `data/intermediate/PHA_tract_bystate/[{State code}]_intersects.RDS`
    
- [02_spatialmerge_CoC.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/02_spatialmerge_CoC.R)
  - Takes in:
    - `data/raw/PHA_gdb/PHAs.gdb`
    - `data/raw/Continuum_of_Care_(CoC)_Grantee_Areas/`: shapefiles of CoC grantees
  - What it does: for each PHA, finds CoC that it intersects with
  - Outputs:
    - `data/intermediate/phapoint_coc_intersect.RDS`
    
- [03_pull_census_tractlevel.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/03_pull_census_tractlevel.R)
  - Takes in:
    - Spatial data from `data/intermediate/PHA_tract_bystate/[{State code}]_intersects.RDS`
    - `creds.yaml` file containing Census API key
  - What it does: specifies demographic variables to pull from the ACS 5-year estimates, uses the `get_ACS()` function in [tidycensus](https://github.com/walkerke/tidycensus/blob/master/man/get_acs.Rd) to pull tract-level counts, uses the ACS codebook to rename those estimated counts, and merges them back onto data at the PHA-tract dyad level (since one PHA can intersect with 1+ tracts)
  - Outputs:
    - `data/intermediate/phas_wrawACScounts_longernames_20221216.RDS`: a PHA-tract dyad level dataset with counts of people in different demographic categories/raw values for attributes like median household income

- [04_add_all_attributes_tobaseHUD_df.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/04_add_all_attributes_tobaseHUD_df.R)
  - Takes in:
    - `data/raw/HAI_map_201014.csv`
    - `data/raw/countypres_2000-2016.csv`
    - `data/raw/ruca2010revised.xlsx`
    - `data/intermediate/phas_wrawACScounts_longernames_20221216.RDS`
    - `data/intermediate/phapoint_coc_intersect.RDS`
  - What it does: reads in, cleans, and merges the contextual attributes of PHAs. See online supplement for more discussion.
  - Outputs: 
    - `data/intermediate/pha_wlocalattributes_final.RDS`
    
- [05_join_covariates.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/05_join_covariates.R)
  - Takes in:
    - `data/raw/reviewed_codes.RDS`: preference data - see paper + online supplement for description of coding
    - `data/raw/VMS_{}`: voucher management system data from HUD
    - `data/raw/a00000009.gdbtable`: housing authority characteristics from HUD
    - `data/intermediate/pha_wlocalattributes_final.RDS`
  - What it does: merges data sources together
  - Outputs:
    - `data/cleaned/pha_all.[RDS|csv]`

### Helper scripts sourced by above

- [utils.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/utils.R)
- [helperfunc_nonprofitdata.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/helperfunc_nonprofitdata.R)

