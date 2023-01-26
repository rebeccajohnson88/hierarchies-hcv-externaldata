# Constructing PHA-level characteristics 🏠

## Description

For `Hierarchies in the Decentralized Welfare State: Prioritization in the Housing Choice Voucher`

Code for acquiring and cleaning the external data sources used to describe PHAs and their local contexts. Analytic replication code is here: https://doi.org/10.7910/DVN/O3YYII

## Authors

- [Simone Zhang](https://simonezhang.com/) 
  - [https://github.com/sxzh](https://github.com/sxzh)
- [Rebecca A. Johnson](https://www.rebeccajohnson.io/) 
  - [https://github.com/rebeccajohnson88](https://github.com/rebeccajohnson88)

## Input data

Available in this Dropbox folder: [add link]

## Code 

`src/` run in order:
- [00a_pulltract_shapefiles.R](https://github.com/rebeccajohnson88/hierarchies-hcv-externaldata/blob/main/src/00a_pulltract_shapefiles.R)

  - Takes in:
    - `data/raw/Estimated_Housing_Authority_Service_Areas`: shapefiles from HUD on PHA estimated service areas
  - What it does: uses the Census API/[tigris wrapper](https://www.rdocumentation.org/packages/tigris/versions/1.6.1/topics/tracts) to pull the tract polygons in each state in which PHAs are located 
  - Outputs:
    - `data/intermediate/tracts_sf_format.RDS`: tract polygons stored in [sf](https://cran.r-project.org/web/packages/sf/index.html) format


