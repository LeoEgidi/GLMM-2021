## ----load-packages----------------------------------------------------
library(rstan)
library(dplyr)
library(lubridate)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default())

# seed for R's pseudo-RNGs, not Stan's
set.seed(1123) 

## ----load-data-----------------------------------------------------------

## Imagine that you are a statistician or data scientist working as an
## independent contractor. One of your clients is a company that owns many
## residential buildings throughout New York City. The property manager
## explains that they are concerned about the number of cockroach complaints
## that they receive from their buildings. Previously the company has offered
## monthly visits from a pest inspector as a solution to this problem.
## While this is the default solution of many property managers in NYC,
## the tenants are rarely home when the inspector visits, and so the manager
## reasons that this is a relatively expensive solution that is currently
## not very effective.

## One alternative to this problem is to deploy long term bait stations.
## In this alternative, child and pet safe bait stations are installed
## throughout the apartment building. Cockroaches obtain quick acting poison
## from these stations and distribute it throughout the colony.
## The manufacturer of these bait stations provides some indication of the
## space-to-bait efficacy, but the manager suspects that this guidance was
## not calculated with NYC roaches in mind. NYC roaches,
## the manager rationalizes, have more hustle than traditional roaches;
## and NYC buildings are built differently than other common residential
## buildings in the US. This is particularly important as the unit cost
## for each bait station per year is quite high.

pest_data <- readRDS('data/pest_data.RDS')
str(pest_data)

## ----describe-data-------------------------------------------------------
N_buildings <- length(unique(pest_data$building_id))
N_buildings

## ----data-plots----------------------------------------------------------
ggplot(pest_data, aes(x = complaints)) + 
  geom_bar()

ggplot(pest_data, aes(x = traps, y = complaints, color = live_in_super == TRUE)) + 
  geom_jitter()

## ---- data-plots-ts, fig.width = 6, fig.height = 8-----------------------
ggplot(pest_data, aes(x = date, y = complaints, color = live_in_super == TRUE)) + 
  geom_line(aes(linetype = "Number of complaints")) + 
  geom_point(color = "black") + 
  geom_line(aes(y = traps, linetype = "Number of traps"), color = "black", size = 0.25) + 
  facet_wrap(~building_id, scales = "free", ncol = 2, labeller = label_both) + 
  scale_x_date(name = "Month", date_labels = "%b") + 
  scale_y_continuous(name = "", limits = range(pest_data$complaints)) + 
  scale_linetype_discrete(name = "") + 
  scale_color_discrete(name = "Live-in super")


## ----stan-data-----------------------------------------------------------
stan_dat_simple <- list(
  N = nrow(pest_data), 
  complaints = pest_data$complaints,
  traps = pest_data$traps
)
str(stan_dat_simple)

## ----fit_P_real_data-----------------------------------------
## write the model in stan using the script "simple_poisson_regression.stan"
## use the following names for the data: N, complaints and traps
## in the generated quantities block compute y_rep for pppc
comp_model_P <- stan_model('stan_programs/simple_poisson_regression.stan')
fit_P_real_data <- sampling(comp_model_P, data = stan_dat_simple)

## We could also compile and fit the model (two lines above)
## with the 'stan' function

fit_P_real_data <- stan('stan_programs/simple_poisson_regression.stan',
                              data = stan_dat_simple)



## ----results_simple_P----------------------------------------------------
print(fit_P_real_data, pars = c('alpha','beta'))

## ----hist_simple_P-------------------------------------------------------
mcmc_hist(as.matrix(fit_P_real_data, pars = c('alpha','beta')))
mcmc_scatter(as.matrix(fit_P_real_data, pars = c('alpha','beta')), alpha = 0.2)

## ----extract_hypothetical_replications --------------------------------------------------------------------
y_rep <- as.matrix(fit_P_real_data, pars = "y_rep")

## ----marginal_PPC--------------------------------------------------------
ppc_dens_overlay(y = stan_dat_simple$complaints, y_rep[1:200,])

## ----proportion_of_zero_statistic--------------------------------------------------------------------
prop_zero <- function(x) mean(x == 0)
ppc_stat(y = stan_dat_simple$complaints, yrep = y_rep, stat = "prop_zero")

## ----standardized_residuals--------------------------------------------------------------------
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$complaints - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)

## ----predictive_intervals--------------------------------------------------------------------
ppc_intervals(
  y = stan_dat_simple$complaints, 
  yrep = y_rep,
  x = stan_dat_simple$traps
) + 
  labs(x = "Number of traps", y = "Number of complaints")

## ------------------------------------------------------------------------
ggplot(pest_data, aes(x = log(total_sq_foot), y = log1p(complaints))) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

## ----add_some_data--------------------------------------------------------------------
stan_dat_simple$log_sq_foot <- log(pest_data$total_sq_foot/1e4)
stan_dat_simple$live_in_super <- pest_data$live_in_super

## ----fit_mult_P_real_dat-------------------------------------------------
## extend the previous model including one more predictor and the offset
## call the extra predictor 'live_in_super' and the offset 'log_sq_foot'
fit_model_P_mult_real <- stan('stan_programs/multiple_poisson_regression.stan',
                              data = stan_dat_simple)
y_rep <- as.matrix(fit_model_P_mult_real, pars = "y_rep")
ppc_dens_overlay(stan_dat_simple$complaints, y_rep[1:200,])

## ------------------------------------------------------------------------
prop_zero <- function(x) mean(x == 0)
ppc_stat(y = stan_dat_simple$complaints, yrep = y_rep, stat = "prop_zero", binwidth = 0.01)

## ------------------------------------------------------------------------
ppc_intervals(
  y = stan_dat_simple$complaints, 
  yrep = y_rep,
  x = stan_dat_simple$traps
) + 
  labs(x = "Number of traps", y = "Number of complaints")

## ----runNB---------------------------------------------------------------
## write the negative binomial model

fitted_model_NB <- stan('stan_programs/multiple_NB_regression.stan',
                        data = stan_dat_simple)
samps_NB <- rstan::extract(fitted_model_NB)

mcmc_areas(fitted_model_NB, pars=c("alpha", "beta"),
          prob = 0.95)
## ----ppc-full------------------------------------------------------------
y_rep <- samps_NB$y_rep
ppc_dens_overlay(stan_dat_simple$complaints, y_rep[1:200,])

## ------------------------------------------------------------------------
ppc_stat(y = stan_dat_simple$complaints, yrep = y_rep, stat = "prop_zero")

## ----standardized_residuals--------------------------------------------------------------------
mean_inv_phi <- mean(samps_NB$inv_phi)
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$complaints - mean_y_rep) / sqrt(mean_y_rep + mean_y_rep^2*mean_inv_phi)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)

## ----predictive_intervals--------------------------------------------------------------------
ppc_intervals(
  y = stan_dat_simple$complaints, 
  yrep = y_rep,
  x = stan_dat_simple$traps
) + 
  labs(x = "Number of traps", y = "Number of complaints")

## ----ppc-group_means-----------------------------------------------------
ppc_stat_grouped(
  y = stan_dat_simple$complaints, 
  yrep = y_rep, 
  group = pest_data$building_id, 
  stat = 'mean',
  binwidth = 0.2
)
