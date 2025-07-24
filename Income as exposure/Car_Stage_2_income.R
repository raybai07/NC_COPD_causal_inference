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

# Read demographic data
nc_data <- read.csv("final_NC_2023_with_gps.csv", header = TRUE)

# Save exposure
A <- nc_data$log_med_income
# Save generalized propensity scores
gps <- nc_data$gps

# Prepare covariate matrix X
X <- as.matrix(nc_data)
X <- X[, -c(1:3,7,14)]  # drop county name, ED visits, total population, log_med_income, gps
X <- as.numeric(X)
X <- matrix(unlist(X), nrow = 100, ncol = 9)

# Scale selected columns (excluding categorical variable in col 6)
X[, c(1:5, 7:9)] <- scale(X[, c(1:5, 7:9)])

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

# Check that we have all 100 counties and numbers look right
length(expected)
length(observed)
expected
observed

# Save data to a list
stan_data <- list(n = n,
                  p = p,
                  y = observed, 
                  W = W,
                  X = X,
                  A = A,
                  gps = gps,
                  a_new = 4,
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

png(filename="stage2_trace_plots_income_exposure.png")
traceplot(resStan, pars = c("beta0","beta_A","beta","beta_r","rho_phi","sigma_phi"), inc_warmup = FALSE)
dev.off()

# Results
resStan_summaries <- summary(resStan)$summary
head(resStan_summaries)

# Log-relative risk
beta_A <- resStan_summaries['beta_A', c('mean','2.5%','97.5%')]
beta_A
# Posterior mean: -1.0914279
# 95% CI is (-1.9244551, -0.2661778)
 
# Exponentiate beta_A to get relative risk
exp(beta_A)
# posterior mean: 0.3357368
# 95% CI is (0.1459553, 0.7663029)

# stan_counterfactuals_4 <- resStan_summaries[grepl("y_a", rownames(resStan_summaries)),  ]
# head(stan_counterfactuals_4)

# stan_mean_counterfactuals_4 <- stan_counterfactuals_4[,1]
# These are the expected counterfactual mean number
# stan_mean_counterfactuals_4

# Average causal effect for all counties with log(income)=4
# resStan_summaries <- summary(resStan)$summary
# mu_4 <- resStan_summaries['mu_a', c('mean','2.5%','97.5%')]
# mu_4

# Get the MCMC samples
mcmc_samples <- as.matrix(resStan)

# Plot the posterior for beta_A (log-relative risk)
pdf("log_RR_income.pdf", height=5, width=5)
plot(density(mcmc_samples[,'beta_A']), main="Log-Relative Risk for Median Income", 
     xlab=expression(beta[A]),
     lwd=3,
     col='blue')
dev.off()

# Plot the posterior for the average causal effect
# plot(density(mcmc_samples[,'mu_a']), main=bquote(paste("Average Causal Effect for Counties with Log(Median Income) of 4")), 
#      xlab="Average Causal Effect", lwd=3, col="red")

# Save posterior samples. Might use them later.
posterior_samples <- rstan::extract(resStan)
write.csv(as.matrix(posterior_samples), file="outcome_regression_posterior_samples.csv", row.names=F)
