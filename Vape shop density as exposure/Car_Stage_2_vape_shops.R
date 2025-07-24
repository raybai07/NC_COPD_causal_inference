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

# Make sure that these have length 81
length(A)
length(gps)

# Prepare confounder matrix X
X <- as.matrix(reduced_nc_data)
X <- X[, -c(1:3,13:14)]  # drop county name, ED visits, total population, vape_shop_density, gps
X <- as.numeric(X)
X <- matrix(unlist(X), nrow = 81, ncol = 9)

# Scale selected columns (excluding categorical variable in col 7)
head(X)
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

# We only need the data for the 81 counties that have nonzero vape shop density
expected <- expected[-zero_shops]
observed <- observed[-zero_shops]

# Ensure that these only have 81 entries
length(expected)
length(observed)

# Save data to a list
stan_data <- list(n = n,
                  p = p,
                  y = observed, 
                  W = W,
                  X = X,
                  A = A,
                  gps = gps,
                  a_new = 2,
                  log_offset = log(expected)) 

# Load and compile Stan model
options(mc.cores = parallel::detectCores())
fileName <- "./outcome_CAR_model.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)

# Run 4 chains of 5000 MCMC iterations per chain and discard first 2000 as warmup,
# leaving us with a total of 12,000 total MCMC samples.
resStan <- stan(model_code = stan_code, data = stan_data,
                chains = 4, iter = 5000, warmup = 2000)

# Check convergence
print(resStan, pars = c("beta0","beta_A","beta","beta_r","rho_phi","sigma_phi"))

png(filename="stage2_trace_plots_vape_shop_exposure.png")
traceplot(resStan, pars = c("beta0","beta_A","beta","beta_r","rho_phi","sigma_phi"), inc_warmup = FALSE)
dev.off()

# Results
resStan_summaries <- summary(resStan)$summary
head(resStan_summaries)

# Log-relative risk
beta_A <- resStan_summaries['beta_A', c('mean','2.5%','97.5%')]
beta_A
# posterior mean: -0.04073938
# 95% CI is (-0.18689537, 0.10633930)

# Exponentiate beta_A to get relative risk
exp(beta_A)
# posterior mean: 0.9600793
# 95% CI is (0.8295305, 1.1121992)

# stan_counterfactuals_2 <- resStan_summaries[grepl("y_a", rownames(resStan_summaries)),  ]
# head(stan_counterfactuals_2)

# stan_mean_counterfactuals_2 <- stan_counterfactuals_2[,1]
# These are the expected counterfactual mean number
# stan_mean_counterfactuals_2

## Average causal effect for all counties with log(vape shop density)=2
# resStan_summaries <- summary(resStan)$summary
# mu_2 <- resStan_summaries['mu_a', c('mean','2.5%','97.5%')]
# mu_2

# Get the MCMC samples
mcmc_samples <- as.matrix(resStan)

# Plot the posterior for beta_A (log-relative risk)
pdf("log_RR_vape_shop_density.pdf", height=5, width=5)
plot(density(mcmc_samples[,'beta_A']), main="Log-Relative Risk for Vape Shope Density", 
     xlab=expression(beta[A]),
     lwd=3,
     col='blue')
dev.off()

# Plot the posterior for the average causal effect
# plot(density(mcmc_samples[,'mu_a']), main=bquote(paste("Average Causal Effect for Counties with Log(Vape shop density) of 0")), 
#      xlab="Average Causal Effect", lwd=3, col="red")

# Save posterior samples. Might use them later.
posterior_samples <- rstan::extract(resStan)
write.csv(as.matrix(posterior_samples), file="outcome_regression_posterior_samples.csv", row.names=F)
