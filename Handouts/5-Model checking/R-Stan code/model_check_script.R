library(rstan)
library(ggplot2)
library(rstanarm)
library(bayesplot)
library(dplyr)
library(lubridate)

########################
## Eight schools
########################


y <- c(28,8,-3,7,-1,1,18,12)
sigma <- c(15,10,16,11,9,11,10,18)
J <- 8

data <- list(y = y, sigma=sigma, J = J)
fit <- stan("8schools.stan", data = data, iter=200,
            cores = 4, chains =4)
sims <- extract(fit)
posterior <- as.matrix(fit)

theta_names <- paste0("theta[", 1:8, "]")
mu_names <- expression(mu)
tau_names <- expression(tau)
par_names <- c(mu_names, tau_names, theta_names)

mcmc_intervals(posterior, regex_pars = c("theta"))
pdf(file="8schools_areas.pdf", width=8, height =7)
mcmc_areas(posterior, regex_pars=c("mu", "tau", "theta"))+
  scale_y_discrete(labels = rev((parse(text= par_names))))+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)
dev.off()

# check 1: density

pdf(file="8schools_dens.pdf", width=8, height =7)
ppc_dens_overlay(y, sims$y_rep)+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

# check 2: ecdf

pdf(file="8schools_ecdf.pdf", width=8, height =7)
ppc_ecdf_overlay(y, sims$y_rep)+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

# check 3: 

pdf(file="8schools_intervals.pdf", width=8, height =7)
ppc_intervals(y, sims$y_rep)+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  labs(x="Schools")+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(4))
dev.off()

# check 4: stats

pdf(file="8schools_mean.pdf", width=5, height=5)
ppc_stat(y, sims$y_rep, stat="mean")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()
pdf(file="8schools_sd.pdf", width=5, height=5)
ppc_stat(y, sims$y_rep, stat="sd")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()
pdf(file="8schools_median.pdf", width=5, height=5)
ppc_stat(y, sims$y_rep, stat="median")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()
pdf(file="8schools_max.pdf", width=5, height=5)
ppc_stat(y, sims$y_rep, stat="max")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()

# check 5: bivariate stats

pdf(file="8schools_bivstats.pdf", width =5, height =5)
ppc_stat_2d(y, sims$y_rep)+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()
pdf(file="8schools_bivstats2.pdf", width =5, height =5)
ppc_stat_2d(y, sims$y_rep, c("median", "max"))+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()

# check 6: binned residuals


ppc_error_binned(y, sims$y_rep[1:8,])


################################
## PEST EXAMPLE
###############################

pest_data <- readRDS('pest_data.RDS')
str(pest_data)


N_buildings <- length(unique(pest_data$building_id))
N_buildings

#à preliminary plots
ggplot(pest_data, aes(x = complaints)) + 
  geom_bar()
#ggsave(file="hist_pest.pdf", width=8,height=6)

ggplot(pest_data, aes(x = traps, y = complaints, color = live_in_super == TRUE)) + 
  geom_jitter()

complaints <- pest_data$complaints
traps <- pest_data$traps

## Simple Poisson

stan_dat_simple <- list(
  N = nrow(pest_data), 
  complaints = pest_data$complaints,
  traps = pest_data$traps
)
comp_model_P <- stan_model('simple_poisson_regression.stan')
fit_model_P <- sampling(comp_model_P, data = stan_dat_simple, seed = 123)
print(fit_model_P, pars = c('alpha','beta'))

y_rep <- as.matrix(fit_model_P, pars = "y_rep")

#densities 

pdf(file="pest_dens.pdf", width =8, height =7)
ppc_dens_overlay(y = stan_dat_simple$complaints, y_rep[1:200,])+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

#proportion of zero

pdf(file="pest_zero.pdf", width =8, height =7)
prop_zero <- function(x) mean(x == 0)
ppc_stat(y = stan_dat_simple$complaints, yrep = y_rep, stat = "prop_zero")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

#residuals

pdf(file="pest_res.pdf", width =8, height =7)
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$complaints - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)+
  labs(x="Mean of y_rep", y= "Stand. residuals")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16))
dev.off()

# intervals

pdf(file="pest_intervals.pdf", width =8, height =7)
ppc_intervals(
  y = stan_dat_simple$complaints, 
  yrep = y_rep,
  x = stan_dat_simple$traps
) + 
  labs(x = "Number of traps", y = "Number of complaints")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
dev.off()


## Negative binomial

comp_model_NB <- stan_model('simple_NB_regression.stan')
fit_model_NB <- sampling(comp_model_NB, data = stan_dat_simple)
print(fit_model_NB, pars = c('alpha','beta'))


samps_NB <- rstan::extract(fit_model_NB)

y_rep <- samps_NB$y_rep
pdf(file="pest_dens_nb.pdf", width =8, height =7)
ppc_dens_overlay(stan_dat_simple$complaints, y_rep[1:200,])+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

#proportion of zero

pdf(file="pest_zero_nb.pdf", width =8, height =7)
prop_zero <- function(x) mean(x == 0)
ppc_stat(y = stan_dat_simple$complaints, yrep = y_rep, stat = "prop_zero")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

#residuals

pdf(file="pest_res_nb.pdf", width =8, height =7)
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_dat_simple$complaints - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)+
  labs(x="Mean of y_rep", y= "Stand. residuals")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16))
dev.off()

## Hierarchical model

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

# Make data list for Stan
stan_dat_hier <-
  with(pest_data,
       list(complaints = complaints,
            traps = traps,
            N = length(traps),
            J = N_buildings,
            M = N_months,
            log_sq_foot = log(pest_data$total_sq_foot/1e4),
            building_data = building_data[,-3],
            mo_idx = as.integer(as.factor(date)),
            K = 4,
            building_idx = building_idx
       )
  )



# centered parametrization

comp_model_NB_hier <- stan_model('hier_NB_regression.stan')
fitted_model_NB_hier <-
  sampling(
    comp_model_NB_hier,
    data = stan_dat_hier,
    chains = 4,
    cores = 4,
    iter = 4000
  )
print(fitted_model_NB_hier, pars = c('sigma_alpha','beta','mu','phi','alpha'))

pdf(file="pest_trace_cp.pdf", width =10, height =7)
mcmc_trace(
  as.array(fitted_model_NB_hier,pars = 'sigma_alpha'),
  np = nuts_params(fitted_model_NB_hier),
  window = c(500,1000)
)
dev.off()

pdf(file="pest_scatter_cp.pdf", width =10, height =7)
scatter_with_divs <- mcmc_scatter(
  as.array(fitted_model_NB_hier),
  pars = c("alpha[4]", 'sigma_alpha'),
  transform = list('sigma_alpha' = "log"),
  np = nuts_params(fitted_model_NB_hier)
)
scatter_with_divs
dev.off()

pdf(file="pest_parallel_cp.pdf", width =10, height =7)
parcoord_with_divs <- mcmc_parcoord(
  as.array(fitted_model_NB_hier, pars = c("sigma_alpha", "alpha")),
  np = nuts_params(fitted_model_NB_hier)
)
parcoord_with_divs
dev.off()

# non centered parametrization

comp_model_NB_hier_ncp <- stan_model('hier_NB_regression_ncp.stan')
fitted_model_NB_hier_ncp <- sampling(comp_model_NB_hier_ncp, data = stan_dat_hier, chains = 4, cores = 4)
print(fitted_model_NB_hier_ncp, pars = c('sigma_alpha','beta','mu','phi','alpha'))
samps_NB_hier_ncp <- rstan::extract(fitted_model_NB_hier_ncp, pars = c('y_rep','inv_phi'))

pdf(file="pest_scatter_ncp.pdf", width =12, height =7)
scatter_no_divs <- mcmc_scatter(
  as.array(fitted_model_NB_hier_ncp),
  pars = c("alpha[4]", 'sigma_alpha'),
  transform = list('sigma_alpha' = "log"),
  np = nuts_params(fitted_model_NB_hier_ncp)
)
bayesplot_grid(scatter_with_divs, scatter_no_divs,
               grid_args = list(ncol = 2), ylim = c(-11, 1))
dev.off()

pdf(file="pest_trace_ncp.pdf", width =10, height =7)
mcmc_trace(
  as.array(fitted_model_NB_hier_ncp,pars = 'sigma_alpha'),
  np = nuts_params(fitted_model_NB_hier),
  window = c(500,1000)
)
dev.off()

pdf(file="pest_parallel_ncp.pdf", width =10, height =7)
parcoord_no_divs <- mcmc_parcoord(
  as.array(fitted_model_NB_hier_ncp, pars = c("sigma_alpha", "alpha")),
  np = nuts_params(fitted_model_NB_hier_ncp)
)
bayesplot_grid(parcoord_with_divs, parcoord_no_divs,
               ylim = c(-3, 3))
dev.off()

# ppcheck: densities
y_rep <- as.matrix(fitted_model_NB_hier_ncp, pars = "y_rep")
color_scheme_set("blue")

pdf(file="pest_dens_hier_nb.pdf", width =8, height =7)
ppc_dens_overlay(stan_dat_hier$complaints, y_rep[1:200,])+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

# ppcheck: statistics

pdf(file="pest_stat_hier_nb.pdf", width=9, height=7)
ppc_stat_grouped(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  group = pest_data$building_id,
  stat = 'mean',
  binwidth = 0.5
)
dev.off()

# ppcheck: proportion of zeros

pdf(file="pest_zero_hier_nb.pdf", width=8, height =7)
prop_zero <- function(x) mean(x == 0)
ppc_stat(
  y = stan_dat_hier$complaints,
  yrep = y_rep,
  stat = prop_zero,
  binwidth = 0.025
)+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

dev.off()
