##---load_packages----------------------------------------------------------
library(rstan)
library(dplyr)
library(lubridate)
library(ggplot2)
library(bayesplot)
library(loo)
library(shinystan)

theme_set(bayesplot::theme_default())

# seed for R's pseudo-RNGs, not Stan's
set.seed(1123) 

## ----load-data-----------------------------------------------------------
pest_data <- readRDS('data/pest_data.RDS')
str(pest_data)

## ----describe-data-------------------------------------------------------
N_buildings <- length(unique(pest_data$building_id))
N_buildings

## ----prep-data-----------------------------------------------------------
N_months <- length(unique(pest_data$date))

# Add some IDs for building and month
pest_data <- pest_data %>%
  mutate(
    building_fac = factor(building_id, levels = unique(building_id)),
    building_idx = as.integer(building_fac),
    ids = rep(1:N_months, N_buildings),
    mo_idx = lubridate::month(date)
  )

# Center and rescale the building specific data
building_data <- pest_data %>%
    select( 
      building_idx,
      live_in_super,
      age_of_building,
      total_sq_foot,
      average_tenant_age,
      monthly_average_rent
    ) %>%
    unique() %>% 
    arrange(building_idx) %>% 
    select(-building_idx) %>% 
    scale(scale=FALSE) %>%
    as.data.frame() %>%
    mutate( # scale by constants
      age_of_building = age_of_building / 10,
      total_sq_foot = total_sq_foot / 10000,
      average_tenant_age = average_tenant_age / 10,
      monthly_average_rent = monthly_average_rent / 1000
    ) %>%
    as.matrix()
str(building_data)
# Make data list for Stan
stan_dat_hier <-
  with(pest_data,
        list(complaints = complaints,
             traps = traps,
             N = length(traps),
             J = N_buildings,
             log_sq_foot = log(pest_data$total_sq_foot/1e4),
             building_data = building_data[,-3],
             K = 4,
             building_idx = building_idx
             )
        )

## ----comp-NB-hier-------------
comp_model_NB_hier <- stan_model('stan_programs/hier_NB_regression.stan')

## ----run-NB-hier---------------------------------------------------------
fitted_model_NB_hier <-
  sampling(
   comp_model_NB_hier,
   data = stan_dat_hier,
   chains = 4,
#   cores = 4,
   iter = 4000
  )

## We could also have used the function 'stan' 

## ----extract_summaries--------------------------------------------------------------------
samps_hier_NB <- rstan::extract(fitted_model_NB_hier)

## ----print-NB-hier-------------------------------------------------------
print(fitted_model_NB_hier, pars = c('sigma_mu','beta','alpha','phi','mu'))

## ------------------------------------------------------------------------
# use as.array to keep the markov chains separate for trace plots
mcmc_trace(
  as.array(fitted_model_NB_hier,pars = 'sigma_mu'),
  np = nuts_params(fitted_model_NB_hier),
  window = c(500,1000)
)

## ------------------------------------------------------------------------
# assign to object so we can compare to another plot later
scatter_with_divs <- mcmc_scatter(
  as.array(fitted_model_NB_hier),
  pars = c("mu[4]", 'sigma_mu'),
  transform = list('sigma_mu' = "log"),
  np = nuts_params(fitted_model_NB_hier)
)
scatter_with_divs

## ----parallel_plot--------------------------------------------------------------------
parcoord_with_divs <- mcmc_parcoord(
  as.array(fitted_model_NB_hier, pars = c("sigma_mu", "mu")),
  np = nuts_params(fitted_model_NB_hier)
)
parcoord_with_divs

## ----run-NB-hier-ncp, non_centered_parametrization-----------------------------------------------------
fitted_model_NB_hier_ncp <- stan('stan_programs/hier_NB_regression_ncp.stan',
                                     data = stan_dat_hier, 
                                     chains = 4, 
#                                     cores = 4,
                                     iter=4000)

## ----n-eff-NB-hier-ncp-check---------------------------------------------
print(fitted_model_NB_hier_ncp, pars = c('sigma_mu','beta','alpha','phi','mu'))

## ------------------------------------------------------------------------
scatter_no_divs <- mcmc_scatter(
  as.array(fitted_model_NB_hier_ncp),
  pars = c("mu[4]", 'sigma_mu'),
  transform = list('sigma_mu' = "log"),
  np = nuts_params(fitted_model_NB_hier_ncp)
)
bayesplot_grid(scatter_with_divs, scatter_no_divs,
               grid_args = list(ncol = 2), ylim = c(-11, 1))

## ----parallel_plot--------------------------------------------------------------------
parcoord_no_divs <- mcmc_parcoord(
  as.array(fitted_model_NB_hier_ncp, pars = c("sigma_mu", "mu")),
  np = nuts_params(fitted_model_NB_hier_ncp)
)
bayesplot_grid(parcoord_with_divs, parcoord_no_divs,
               ylim = c(-3, 3))

## ----samps-full-hier-----------------------------------------------------
samps_NB_hier_ncp <- rstan::extract(fitted_model_NB_hier_ncp, 
                                    pars = c('y_rep','inv_phi'))

## ----ppc-full-hier-------------------------------------------------------
y_rep <- as.matrix(fitted_model_NB_hier_ncp, pars = "y_rep")
ppc_dens_overlay(stan_dat_hier$complaints, y_rep[1:200,])

## ----ppc-group_means-hier------------------------------------------------
ppc_stat_grouped(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  group = pest_data$building_id,
  stat = 'mean',
  binwidth = 0.5
)


## ----proportion_of_zero--------------------------------------------------------------------
prop_zero <- function(x) mean(x == 0)
ppc_stat(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  stat = prop_zero,
  binwidth = 0.025
)

# plot separately for each building
ppc_stat_grouped(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  group = pest_data$building_id,
  stat = prop_zero,
  binwidth = 0.025
)

## ----ppc_intervals--------------------------------------------------------------------
ppc_intervals(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  x = stan_dat_hier$traps
) +
  labs(x = "Number of traps", y = "Number of complaints")

## ----standardized_residuals--------------------------------------------------------------------
mean_y_rep <- colMeans(y_rep)
mean_inv_phi <- mean(as.matrix(fitted_model_NB_hier_ncp, pars = "inv_phi"))
std_resid <- (stan_dat_hier$complaints - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)

## ------------------------------------------------------------------------
stan_dat_hier <- readRDS('data/pest_data_longer_stan_dat.RDS')

## ----run-NB-hier-slopes--------------------------------------------------
fitted_model_NB_hier_slopes <-
  stan(
    'stan_programs/hier_NB_regression_ncp_slopes_mod.stan',
    data = stan_dat_hier,
    chains = 4, #cores = 4,
    control = list(adapt_delta = 0.95)
  )

## ------------------------------------------------------------------------
mcmc_hist(
  as.matrix(fitted_model_NB_hier_slopes, pars = "sigma_kappa"),
  binwidth = 0.005
)

## ------------------------------------------------------------------------
print(fitted_model_NB_hier_slopes, pars = c('kappa','beta','alpha','phi','sigma_mu','sigma_kappa','mu'))

## ------------------------------------------------------------------------
mcmc_hist(
  as.matrix(fitted_model_NB_hier_slopes, pars = "beta"),
  binwidth = 0.005
)

## ----ppc-full-hier-slopes------------------------------------------------
y_rep <- as.matrix(fitted_model_NB_hier_slopes, pars = "y_rep")
ppc_dens_overlay(
  y = stan_dat_hier$complaints,
  yrep = y_rep[1:200,]
)


##----model_comparison----------------------------------------------------
## We are considering here the EXTENDED dataset

## Extract pointwise log-likelihood and
## compute loo and waic.

# data 
pest_data <- readRDS('data/pest_data_longer_stan_dat.RDS')
str(pest_data)

## ----describe-data-------------------------------------------------------
N_buildings <- length(unique(pest_data$building_idx))
N_buildings 

## Data acquisition

stan_dat_simple <- list(
  N = length(pest_data$complaints), 
  complaints = pest_data$complaints,
  traps = pest_data$traps
)


## Models fit

# simple poisson
comp_model_P <- stan_model('stan_programs/simple_poisson_regression.stan')
fit_P_real_data <- sampling(comp_model_P, data = stan_dat_simple)

# multiple poisson
stan_dat_simple$log_sq_foot <- pest_data$log_sq_foot
stan_dat_simple$live_in_super <- pest_data$pred_mat[,"live_in_super"]
comp_model_P_mult <- stan_model('stan_programs/multiple_poisson_regression.stan')
fit_model_P_mult_real <- sampling(comp_model_P_mult, data = stan_dat_simple)

# negative binomial
comp_model_NB <- stan_model('stan_programs/multiple_NB_regression.stan')
fitted_model_NB <- sampling(comp_model_NB, data = stan_dat_simple)

# -----------------------
stan_dat_hier <- pest_data
# hier NB regression
comp_model_NB_hier <- stan_model('stan_programs/hier_NB_regression.stan')
fitted_model_NB_hier <-
  sampling(
    comp_model_NB_hier,
    data = stan_dat_hier,
    chains = 4,
    cores = 4,
    iter = 4000
  )

# hier NCP NB regression
comp_model_NB_hier_ncp <- stan_model('stan_programs/hier_NB_regression_ncp.stan')
fitted_model_NB_hier_ncp <- sampling(comp_model_NB_hier_ncp, 
                                     data = stan_dat_hier, 
                                     chains = 4, 
                                     cores = 4,
                                     iter=4000)

# hier NB slopes
comp_model_NB_hier_slopes <- stan_model('stan_programs/hier_NB_regression_ncp_slopes_mod.stan')
fitted_model_NB_hier_slopes <-
  sampling(
    comp_model_NB_hier_slopes,
    data = stan_dat_hier,
    chains = 4, cores = 4,
    control = list(adapt_delta = 0.95)
  )

## Shinystan

y <- pest_data$complaints
launch_shinystan(fitted_model_NB_hier)

## Extract pointwise log-likelihood and
## compute loo and waic

# loo
log_lik_pois <- extract_log_lik(fit_P_real_data)
loo_pois <- loo(log_lik_pois)
waic_pois <- waic(log_lik_pois)

log_lik_mult_pois <- extract_log_lik(fit_model_P_mult_real)
loo_mult_pois <- loo(log_lik_mult_pois)
waic_mult_pois <- waic(log_lik_mult_pois)


log_lik_mult_NB <- extract_log_lik(fitted_model_NB)
loo_mult_NB <- loo(log_lik_mult_NB)
waic_mult_NB <- waic(log_lik_mult_NB)


log_lik_hier_NB <- extract_log_lik(fitted_model_NB_hier)
loo_hier_NB <- loo(log_lik_hier_NB)
waic_hier_NB <- waic(log_lik_hier_NB)

log_lik_hier_NB_ncp <- extract_log_lik(fitted_model_NB_hier_ncp)
loo_hier_NB_ncp <- loo(log_lik_hier_NB_ncp)
waic_hier_NB_ncp <- waic(log_lik_hier_NB_ncp)

log_lik_slopes <- extract_log_lik(fitted_model_NB_hier_slopes)
loo_slopes <- loo(log_lik_slopes)
waic_slopes <- waic(log_lik_slopes)


compare(loo_pois, loo_mult_pois,
        loo_mult_NB, loo_hier_NB,
        loo_hier_NB_ncp, loo_slopes)

looic1<-loo_pois$estimates[3,1]
looic2<-loo_mult_pois$estimates[3,1]
looic3<-loo_mult_NB$estimates[3,1]
looic4<-loo_hier_NB$estimates[3,1]
looic5<-loo_hier_NB_ncp$estimates[3,1]
looic6<-loo_slopes$estimates[3,1]
looics <-c(looic1, looic2, looic3, looic4, looic5,
           looic6)

waic1<-waic_pois$estimates[3,1]
waic2<-waic_mult_pois$estimates[3,1]
waic3<-waic_mult_NB$estimates[3,1]
waic4<-waic_hier_NB$estimates[3,1]
waic5<-waic_hier_NB_ncp$estimates[3,1]
waic6<-waic_slopes$estimates[3,1]
waics <-c(waic1, waic2, waic3, waic4, waic5,
          waic6)

par(xaxt="n", mfrow=c(1,2))
plot(looics, type="b", xlab="", ylab="LOOIC")
par(xaxt="s")
axis(1, c(1:6), c("Pois", "Pois mult", "NB",
                  "Hier NB", "Hier NB ncp", "NB slopes"), las=2)
par(xaxt="n")
plot(waics, type="b", xlab="", ylab="WAIC")
par(xaxt="s")
axis(1, c(1:6), c("Pois", "Pois mult", "NB",
                  "Hier NB", "Hier NB ncp", "NB slopes"), las=2)



