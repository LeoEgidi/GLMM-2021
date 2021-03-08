##############################
## 8 schools example
##############################

library(rstan)
library(bayesplot)
schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
fit_8schools <- stan(file = '8schools.stan', data = schools_dat)
print(fit_8schools, pars=c("mu", "tau", "theta"))
posterior <- as.array(fit_8schools)
color_scheme_set("red")

par_names <- c()
par_names[1]<-expression(mu)
par_names[2]<-expression(tau)

par_names[3]<-expression(theta[8])
par_names[4]<-expression(theta[7])
par_names[5]<-expression(theta[6])
par_names[6]<-expression(theta[5])
par_names[7]<-expression(theta[4])
par_names[8]<-expression(theta[3])
par_names[9]<-expression(theta[2])
par_names[10]<-expression(theta[1])

par_names[11]<-expression(eta[1])
par_names[12]<-expression(eta[2])
par_names[13]<-expression(eta[3])
par_names[14]<-expression(eta[4])
par_names[15]<-expression(eta[5])
par_names[16]<-expression(eta[6])
par_names[17]<-expression(eta[7])
par_names[18]<-expression(eta[8])



# posterior intervals
pdf(file="post_int_8schools.pdf", width =9, height=8.5)
mcmc_intervals(posterior,regex_pars=c("theta","tau", "mu" ))+
  scale_y_discrete(labels = rev((parse(text= par_names[1:10]))))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))+
  ggtitle("Posterior intervals")
dev.off()

# posterior areas
pdf(file="post_areas_8schools.pdf", width =9, height=8.5)
mcmc_areas(posterior, pars=c(  "theta[1]", "theta[2]", 
                               "theta[3]", "theta[4]", "theta[5]", "theta[6]",
                               "theta[7]", "theta[8]",
                               "tau", "mu" ))+
  scale_y_discrete(labels = rev((parse(text= par_names[1:10]))))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))+
  ggtitle("Posterior areas")
dev.off()


# marginal posterior
pdf(file="marg_post_8schools.pdf", width=12, height =8)
mcmc_dens(posterior)
dev.off()

# marginal posterior overlayed
pdf(file="marg_post_8schools_4chains.pdf", width=12.4, height =8)
mcmc_dens_overlay(posterior)
dev.off()

# bivariate plots
pdf(file="pairs_post_8schools.pdf", width=10, height=8)
mcmc_pairs(posterior, pars=c("mu", "tau"))
dev.off()

# trace plots
pdf(file="trace_post_8schools_4chains.pdf", width=12.4, height =8)
mcmc_trace(posterior)
dev.off()


