# Authors: Evan Funderburg, Emily Sitnik, Mark Ritchie, and Ray Bai

# Clear R global environment
rm(list=ls())

# Load libraries
library(rstan)
library(parallel)

#Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(sf)
library(viridis)
library(spdep)

# Read shapefile
nc_map <- st_read("tl_2018_37_cousub.shp")
nc_map <- nc_map |>
  group_by(COUNTYFP) |>
  summarise(geometry = st_union(geometry))

# Create neighborhood matrix
neighbors <- poly2nb(nc_map)
W <- nb2mat(neighbors, style = "B")

# ✅ Force symmetry in W (important for Stan CAR model)
W <- (W + t(W)) / 2

# Remove rows and columns for counties with zero vape shops
nc_data <- read.csv("final_NC_2023_county_data.csv", header = TRUE)
zero_shops <- which(nc_data$vape_shop_density==0)
# Update W
W <- W[-zero_shops, -zero_shops]

# Load the data for only the counties with nonzero vape shop density
reduced_nc_data <- read.csv("reduced_NC_data_2023_with_gps.csv", header = TRUE)

# Save vape shop density
A <- reduced_nc_data$vape_shop_density
# Save generalized propensity scores
gps <- reduced_nc_data$gps

# Prepare covariate matrix X
X <- as.matrix(reduced_nc_data)
X <- X[, -c(1:3,13:14)]  # drop county name, ED visits, total population, vape_shop_density, gps
X <- as.numeric(X)
X <- matrix(unlist(X), nrow = 81, ncol = 9)

# Scale selected columns (excluding categorical variable in col 7)
X[, c(1:6, 8:9)] <- scale(X[, c(1:6, 8:9)])

n <- dim(X)[1]   # number of observations
p <- dim(X)[2]   # number of confounders

# Get response data
count_data <- read.csv("NC_SMR_AGE_COPD.csv", header=T)

dim(count_data)

# Look at first 6 columns of data
count_data[, c(1:6)]

# expected counts is in row 6, observed counts is in row 7

# Change column 'X' to be row.names
row.names(count_data) <- count_data[,1]
# Remove column 'X'
count_data <- count_data[,-1]
# Check
count_data[, c(1:6)]
dim(count_data)

colnames(count_data) # to make sure McDowell is the first county that begins with "M"

expected <- as.numeric(count_data[6, ])
observed <- as.numeric(count_data[7, ])

# We just need the data for the 81 counties that have nonzero vape shop density
expected <- expected[-zero_shops]
observed <- observed[-zero_shops]
expected
observed

min(A)
max(A)

# Let's get the curve for a from min(A) to max(A)
a_levels <- seq(from=min(A), to=max(A), length.out=15)

# Initialize vectors to store results
mu_a_mean <- rep(0, length(a_levels))
mu_a_lower <- rep(0, length(a_levels))
mu_a_upper <- rep(0, length(a_levels))


# Loop through and store results to mu_a, a_lower, and a_upper
for(i in 1:length(a_levels)) {

  print(paste("Exposure level a = ", a_levels[i]))
  
  # Save data to a list. Note that we set a_new = a_levels[i], 
  # so it estimates the average causal effect for the ith exposure level
  # in the vector a_levels
  dat <- list(n = n,
              p = p,
              y = observed,                 # observed number of cases
              W = W,
              X = X,
              A = A,
              gps = gps,
              a_new = a_levels[i],
              log_offset = log(expected)) 

  options(mc.cores = parallel::detectCores())
  fileName <- "./outcome_CAR_model.stan"
  stan_code <- readChar(fileName, file.info(fileName)$size)
  
  # Run Stan code
  # We'll run 3 chains for 3000 iterations and discard the first 1000 warm-up  
  # iterations as burn-in. Thus, we will have 6000 total MCMC samples
  resStan <- stan(model_code = stan_code, data = dat,
                  chains = 4, iter = 3000, warmup = 1000, refresh=100)

  ## Store results
  resStan_summaries <- summary(resStan)$summary
  
  # average causal effect
  mu_a_mean[i] <- resStan_summaries['mu_a', 'mean']
  mu_a_lower[i] <- resStan_summaries['mu_a', '2.5%']
  mu_a_upper[i] <- resStan_summaries['mu_a', '97.5%']
  
}

# Save results to a dataframe
results_df <- data.frame(cbind(a_levels, mu_a_mean, mu_a_lower, mu_a_upper))
colnames(results_df) <- c("a", "mu_a", "lowerCI", "upperCI")
head(results_df)

write.csv(results_df, "curve_results_vape_shops.csv", row.names=FALSE)

# Plotting
min_y <- min(mu_a_lower)
max_y <- max(mu_a_upper)

pdf(file="causal_ER_curve_vape_shops.pdf", width=6, height=6)
ggplot(results_df, aes(x = a, y = mu_a)) + theme_bw() +
  geom_line(color="red", size=1.25) +  # Plot the central line
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha = 0.2) + # Plot the confidence bands
  ggtitle("Causal Exposure-Response Curve with\n Vape Shop Density as the Exposure") + 
  theme(plot.title = element_text(hjust = 0.5, size=20)) +
  labs(y=expression(mu~(a)), x=expression(a)) + 
  coord_cartesian(ylim = c(0,5000)) 
dev.off()