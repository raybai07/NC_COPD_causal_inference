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

# Read demographic data
nc_data <- read.csv("final_NC_2023_county_data.csv", header = TRUE)

# Prepare covariate matrix X
X <- as.matrix(nc_data)
X <- X[, -c(1:3,7)]  # drop county name, ED visits, total population, log_med_income
X <- as.numeric(X)
X <- matrix(unlist(X), nrow = 100, ncol = 9)

# Scale selected columns (excluding categorical variable in col 6)
X[, c(1:5, 7:9)] <- scale(X[, c(1:5, 7:9)])

# Sanity checks
mean(X[,1]); sd(X[,1])
mean(X[,9]); sd(X[,9])
X[,6]  # categorical variable, should not be scaled

# Read COPD count and expected data
count_data <- read.csv("NC_SMR_AGE_COPD.csv", header = TRUE)
row.names(count_data) <- count_data[, 1]
count_data <- count_data[, -1]
expected <- as.numeric(count_data[6, ])
observed <- as.numeric(count_data[7, ])

# Re-read data in case overwritten
nc_data <- read.csv("final_NC_2023_county_data.csv", header = TRUE)

# Extract treatment variable (median income)
A <- nc_data$log_med_income
head(A)

# Prepare data for Stan
n <- nrow(X)
p <- ncol(X)

stan_data <- list(n = n,
                  p = p,
                  X = X,
                  A = as.numeric(A),
                  W = W)

# Load and compile Stan model
options(mc.cores = parallel::detectCores())
fileName <- "./GPS_CAR_model.stan"
stan_code <- readChar(fileName, file.info(fileName)$size)

# Run 4 chains of 5000 MCMC iterations per chain and discard first 2000 as warmup,
# leaving us with a total of 12,000 total MCMC samples.
resStan <- stan(model_code = stan_code, data = stan_data,
                chains = 4, iter = 5000, warmup = 2000)

# Check convergence
print(resStan, pars = c("gamma0", "gamma", "sigma_e", "rho_e"))

png(filename="stage1_trace_plots_income_exposure.png")
traceplot(resStan, pars = c("gamma0", "gamma", "sigma_e", "rho_e"))
dev.off()

# Extract posterior samples
posterior_samples <- rstan::extract(resStan)

# Save posterior samples. Might use them later.
write.csv(as.matrix(posterior_samples), file="GPS_posterior_samples.csv", row.names=F)

# Check the normality assumption
gamma0_hat <- median(posterior_samples$gamma0)
gamma_hat <- apply(posterior_samples$gamma, 2, median)
# Residuals
resids <- A - (gamma0_hat + crossprod(t(X),gamma_hat))

pdf(file="income_histograms.pdf", width=10, height=4)
par(mfrow=c(1,3))
hist(exp(A), main="Median Income", xlab="Median income", prob=TRUE)
hist(A, main="Log-transformed Median Income", xlab="Median income (log-scale)", prob=TRUE)
hist(resids, main="Residuals for GPS Model", xlab="Residuals", prob=TRUE)
dev.off()

# Get the generalized propensity scores

# gps is an n x S matrix (S = number of posterior draws)
gps_samples <- posterior_samples$gps

# Get posterior means (or medians) across MCMC draws
gps_estimates <- apply(gps_samples, 2, mean)  # length n = 100

# Optional: look at first few values
head(gps_estimates)

# Created Updated File with GPS File 
# 1. Read original data (if not already loaded)
nc_data <- read.csv("final_NC_2023_county_data.csv", header = TRUE)

# 2. Add GPS estimates as new column
nc_data$gps <- gps_estimates  # gps_estimates is from your earlier code

# 3. Optional: check first few rows
head(nc_data[, c("County", "log_med_income", "gps")])

# 4. Write to new CSV
write.csv(nc_data, file = "final_NC_2023_with_gps.csv", row.names = FALSE)
