

## Packages
library(here)

### Data
library(tidyverse)
library(ggplot2)
library(stargazer)
library(readr)
library(data.table)
library(readxl)
library(httr)

### Spatial
library(rgdal)
library(cdlTools)
library(sf)
library(tigris)
library(sp)



## plotting theme
theme_new <- function(base_size = 16, base_family = "Helvetica"){
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.grid = element_blank(),   
      panel.border = element_rect(fill = NA, colour = "black", size=1),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = NA),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    )
}
