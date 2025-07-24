# Authors: Evan Funderburg, Emily Sitnik, Mark Ritchie, and Ray Bai

# Clear R global environment
rm(list = ls())

# Load libraries
library(rstan)
library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(sf)
library(viridis)
library(spdep)

# Read data
copd_dat<-read.csv("final_NC_2023_county_data.csv", header=T)
head(copd_dat)

# Add column for ED Visit rate per 100,000 residents
copd_dat$ED_visit_rate <- (copd_dat$ED_visit_count / copd_dat$total_pop) * 100000
# Add column for non-exponentiated median income in thousands
copd_dat$med_income <- exp(copd_dat$log_med_income)
# Check
head(copd_dat)

# Read shapefile
nc_map <- st_read("tl_2018_37_cousub.shp")
nc_map <- nc_map |>
  group_by(COUNTYFP) |>
  summarise(geometry = st_union(geometry))

# Merge county identifiers (countyfp) with corresponding data
nc_countyfp<-read.csv("nc_countyfp.csv", header=T)
copd_dat$COUNTYFP<-nc_countyfp$COUNTYFP
nc_map$COUNTYFP<-as.numeric(nc_map$COUNTYFP)

#Merge geometry data with COPD data
nc_map<-left_join(nc_map, copd_dat, by = "COUNTYFP")

# Make maps
#Plot for ED visit rate
ED_visit_map <- ggplot(nc_map) +
                geom_sf(aes(fill = ED_visit_rate), color = "black") + #Black describes border color
                scale_fill_viridis_c() +  #Color package
                theme_minimal() + #Optional GGplot setting
                labs(title = "COPD ED Visit Rates by NC County", fill = "ED Visit Rate")

ggsave("ED_visit_map.png", plot = ED_visit_map, width = 8, height = 3, dpi = 300, bg="white")

#Plot for Median Income
median_income_map <- ggplot(nc_map) +
                     geom_sf(aes(fill = med_income), color = "black") + #Black describes border color
                     scale_fill_viridis_c() +  #Color package
                     theme_minimal() + #Optional GGplot setting
                     labs(title = "Median Income by NC County", fill = "Median Income")

ggsave("median_income_map.png", plot = median_income_map, width = 8, height = 3, dpi = 300, bg="white")

#Plot for Vape shop density
vape_shop_map <- ggplot(nc_map) +
                 geom_sf(aes(fill = vape_shop_density), color = "black") + #Black describes border color
                 scale_fill_viridis_c() +  #Color package
                 theme_minimal() + #Optional GGplot setting
                 labs(title = "Vape Shop Density by NC County", fill = "Vape Shop Density")

ggsave("vape_shop_map.png", plot = vape_shop_map, width = 8, height = 3, dpi = 300, bg="white")


# Test for significant spatial autocorrelation
library(spdep)

nb <- poly2nb(nc_map, queen = TRUE) # a single shared boundary point meets the contiguity condition
nbw <- nb2listw(nb, style = "W")

# Test for median income
gmoran_log_med_income <- moran.test(nc_map$log_med_income, nbw, alternative = "greater")
gmoran_log_med_income
# Moran's I statistic: 0.37687393, p-value = 3.45e-10

# Test for vape shop density
# Only check for the counties with nonzero vape shop density
zero_shops <- which(nc_map$vape_shop_density==0)
reduced_nc_map <- nc_map[-zero_shops, ]

reduced_nb <- poly2nb(reduced_nc_map, queen = TRUE) # a single shared boundary point meets the contiguity condition
reduced_nbw <- nb2listw(reduced_nb, style = "W")

gmoran_log_vape_shop_density <- moran.test(log(reduced_nc_map$vape_shop_density), reduced_nbw, alternative = "greater")
gmoran_log_vape_shop_density
# Moran's I statistic: 0.131523611, p-value = 0.0327

# One-sided Global Moran's I tests for spatial autocorrelation
# Test for ED visit rate
gmoran_ED_visit_rate <- moran.test(nc_map$ED_visit_rate, nbw, alternative = "greater")
gmoran_ED_visit_rate
# Moran's I statistic: 0.270675467, p-value = 3.221e-06