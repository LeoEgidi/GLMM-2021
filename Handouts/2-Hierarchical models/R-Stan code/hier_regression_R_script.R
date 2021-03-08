library (arm)
library(R2WinBUGS)
library(rstanarm)
library(ggplot2)
library(bayesplot)



#####################################
# HIERARCHICAL LINEAR MODELS
#####################################

## Radon data

srrs2 <- read.table ("srrs2.dat", header=T, sep=",")
mn <- srrs2$state=="MN"
radon <- srrs2$activity[mn]
log.radon <- log (ifelse (radon==0, .1, radon))
floor <- srrs2$floor[mn]       # 0 for basement, 1 for first floor
n <- length(radon)
y <- log.radon
x <- floor


## Partial pooling with no predictors

# get county index variable
county.name <- as.vector(srrs2$county[mn])
uniq <- unique(county.name)
J <- length(uniq)
county <- rep (NA, J)
for (i in 1:J){
  county[county.name==uniq[i]] <- i
}

# no predictors
ybarbar = mean(y)

sample.size <- as.vector (table (county))
sample.size.jittered <- sample.size*exp (runif (J, -.1, .1))
cty.mns = tapply(y,county,mean)
cty.vars = tapply(y,county,var)
cty.sds = mean(sqrt(cty.vars[!is.na(cty.vars)]))/sqrt(sample.size)
cty.sds.sep = sqrt(tapply(y,county,var)/sample.size)

# varying-intercept model, no predictors


mlm.radon.nopred <- stan_lmer(y ~ 1+ (1|county))

## Figure 12.1 (a)

pdf(file="radon_no_hier.pdf", width =11, height = 7)
par(mfrow=c(1,2))
plot (sample.size.jittered, cty.mns, cex.lab=1.6, cex.axis=1,
      xlab="sample size in county j",
      ylab="avg. log radon in county j",
      pch=20, log="x", cex=.8, mgp=c(1.5,.5,0),
      ylim=c(0,3.2), yaxt="n", xaxt="n", cex.main = 1.8)
axis (1, c(1,3,10,30,100), cex.axis=.9, mgp=c(1.5,.5,0))
axis (2, seq(0,3), cex.axis=.9, mgp=c(1.5,.5,0))
for (j in 1:J){
  lines (rep(sample.size.jittered[j],2),
         cty.mns[j] + c(-1,1)*cty.sds[j], lwd=.5)
  #         cty.mns[j] + c(-1,1)*mean(cty.sds[!is.na(cty.sds)]), lwd=.5)
}
abline(h=mlm.radon.nopred$coefficients[1])
title("No pooling",cex.main=1.8, line=1)
#abline(h=ybarbar)
points(sample.size.jittered[36],cty.mns[36],cex=4)

## Figure 12.1 (b)

plot (sample.size.jittered, 
      mlm.radon.nopred$coefficients[1]+ mlm.radon.nopred$coefficients[2:86], 
      cex.lab=1.6, cex.axis=1,
      xlab="sample size in county j",
      ylab="avg. log radon in county j",
      pch=20, log="x", cex=.8, mgp=c(1.5,.5,0),
      ylim=c(0,3.2), yaxt="n", xaxt="n")
axis (1, c(1,3,10,30,100), cex.axis=.9, mgp=c(1.5,.5,0))
axis (2, seq(0,3), cex.axis=.9, mgp=c(1.5,.5,0))
for (j in 1:J){
  lines (rep(sample.size.jittered[j],2),
         mlm.radon.nopred$coefficients[1]+mlm.radon.nopred$coefficients[j+1] + c(-1,1)*mlm.radon.nopred$ses[j+1],
         lwd=.5)
}
abline(h=mlm.radon.nopred$coefficients[1])
points(sample.size.jittered[36],
       mlm.radon.nopred$coefficients[1]+mlm.radon.nopred$coefficients[37],cex=4)#,col="red")
title("Multilevel model",cex.main=1.8, line=1)
dev.off()


## Partial pooling with predictors

## Complete pooling regression
lm.pooled <- lm (y ~ x)
display (lm.pooled)

## No pooling regression
lm.unpooled <- lm (y ~ x + factor(county) -1)
display (lm.unpooled)

## Comparing-complete pooling & no-pooling (Figure 12.2)
pdf(file="radon_8counties.pdf", height =6, width =12)
x.jitter <- x + runif(n,-.05,.05)
display8 <- c (36, 1, 35, 21, 14, 71, 61, 70)  # counties to be displayed
y.range <- range (y[!is.na(match(county,display8))])

par (mfrow=c(2,4), mar=c(3,2,3,1), oma=c(1,1,2,1))
for (j in display8){
  plot (x.jitter[county==j], y[county==j], xlim=c(-.05,1.05), ylim=y.range,
        xlab="floor", ylab="log radon level", cex.lab=1.8, cex.axis=1.1,
        pch=20, mgp=c(2,.7,0), xaxt="n", yaxt="n", cex.main=1.6,
        main=uniq[j], lwd=1.5)
  axis (1, c(0,1), mgp=c(2,.7,0), cex.axis=1.1)
  axis (2, seq(-1,3,2), mgp=c(2,.7,0), cex.axis=1.1)
  curve (coef(lm.pooled)[1] + coef(lm.pooled)[2]*x, lwd=3, lty=2, add=TRUE)
  curve (coef(lm.unpooled)[j+1] + coef(lm.unpooled)[1]*x, lwd=3, add=TRUE)
}
dev.off()

# varying-intercept model, with predictors


mlm.radon.pred <- stan_lmer(y ~ x+ (1|county))
print(mlm.radon.pred)
a.hat <- coefficients(mlm.radon.pred)$county[,1]
b.hat <- coefficients(mlm.radon.pred)$county[,2]


# plot wit the multilevel estimates (12.4)
pdf(file="radon_8counties_mlm.pdf", height =6, width =12)
x.jitter <- x + runif(n,-.05,.05)
display8 <- c (36, 1, 35, 21, 14, 71, 61, 70)  # counties to be displayed
y.range <- range (y[!is.na(match(county,display8))])

par (mfrow=c(2,4), mar=c(3,2,3,1), oma=c(1,1,2,1))
for (j in display8){
  plot (x.jitter[county==j], y[county==j], xlim=c(-.05,1.05), ylim=y.range,
        xlab="floor", ylab="log radon level", cex.lab=1.8, cex.axis=1.1,
        pch=20, mgp=c(2,.7,0), xaxt="n", yaxt="n", cex.main=1.6,
        main=uniq[j], lwd=1.5)
  axis (1, c(0,1), mgp=c(2,.7,0), cex.axis=1.1)
  axis (2, seq(-1,3,2), mgp=c(2,.7,0), cex.axis=1.1)
  curve (coef(lm.pooled)[1] + coef(lm.pooled)[2]*x, lwd=3, lty=2, add=TRUE)
  curve (coef(lm.unpooled)[j+1] + coef(lm.unpooled)[1]*x, lwd=3, add=TRUE)
  curve (a.hat[j] + b.hat[j]*x,  col="red", lwd =3, add=TRUE)
  }
dev.off()

## Eight schools


library(rstan)
library(invgamma)
library(ggplot2)
theme_set(theme_minimal())
library(gridExtra)
library(tidyr)
library(rprojroot)
#root<-has_dirname("BDA_R_demos-master/BDA_R_demos")$make_fix_file()



y <- c(28,8,-3,7,-1,1,18,12)
s <- c(15,10,16,11,9,11,10,18)

# separate analysis

pdf(file ="eightschool_sep.pdf", width =10, height=7)
x <- seq(-40, 60, length.out = 500)
df_sep <- mapply(function(y, s, x) dnorm(x, y, s), y, s, MoreArgs = list(x = x)) %>%
  as.data.frame() %>% setNames(LETTERS[1:8]) %>% cbind(x) %>% gather(school, p, -x)
labs1 <- c('Other Schools', 'School A')
plot_sep <- ggplot(data = df_sep) +
  geom_line(aes(x = x, y = p, color = (school=='A'), group = school)) +
  labs(x = 'Treatment effect', y = '', title = 'Separate model', color = '', size =rel(2)) +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))+
  theme(plot.title = element_text(hjust = 0.5, size =rel(2)),
        axis.title=element_text(size=22))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))
plot_sep
dev.off()

# complete pool

df_pool <- data.frame(x = x, p = dnorm(x, sum(y/s^2)/sum(1/s^2), sqrt(1/sum(1/s^2))))

pdf(file ="eightschool_compl.pdf", width =10, height=7)
plot_pool <- ggplot(data = df_pool) +
  geom_line(aes(x = x, y = p, color = '1')) +
  labs(x = 'Treatment effect', y = '', title = 'Pooled model', color = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = 'red', labels = 'All schools') +
  theme(legend.background = element_blank(), legend.position = c(0.7,0.9))+
  theme(plot.title = element_text(hjust = 0.5, size =rel(2)),
        axis.title=element_text(size=22))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))
plot_pool
dev.off()

# Load the pre-computed results for the hierarchical model Replace this with your own code in the related exercise

load("demo5_2.RData")
#  hierarchical model

pdf(file ="eightschool_hier.pdf", width =10, height=7)
df_hier <- as.data.frame(t(pxm)) %>% setNames(LETTERS[1:8]) %>%
  cbind(x) %>% gather(school, p, -x)
plot_hier <- ggplot(data = df_hier) +
  geom_line(aes(x = x, y = p, color = (school=='A'), group = school)) +
  labs(x = 'Treatment effect', y = '', title = 'Hierarchical model', color = '') +
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))+
  theme(plot.title = element_text(hjust = 0.5, size =rel(2)),
        axis.title=element_text(size=22))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))
plot_hier
dev.off()

#Plot separate, pooled, and hierarchical model

pdf("eightschools_all.pdf", width =11.5, height =8)
grid.arrange(plot_sep, plot_pool, plot_hier)
dev.off()

#Various marginal and conditional posterior summaries

pdf(file="eightschool_margtau.pdf", width =10, height =7)
df_margpost = data.frame(x = t(tt), p = t(tp))
title1 <- 'Marginal posterior density p(tau|y)'
plot_margpost <-
  ggplot(data = df_margpost) +
  geom_line(aes(x = x, y = p)) +
  labs(x = expression(tau), y = 'p(tau|y)', title = title1) +
  scale_y_continuous(breaks = NULL)+
  theme(plot.title = element_text(hjust = 0.5, size =rel(2)),
        axis.title=element_text(size=22))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9))
plot_margpost
dev.off()

df_condmeans <- as.data.frame(t(tm)) %>% setNames(LETTERS[1:8]) %>%
  cbind(x = t(tt)) %>% gather(school, p, -x)

pdf(file="eightschool_condmeans.pdf", width =10, height =7)
yl <- c(-5, 40)
title2 <- 'Conditional posterior means of effects E[theta_j|tau,y]'
plot_condmeans <- ggplot(data = df_condmeans) +
  geom_line(aes(x = x, y = p, color = (school=='A'), group = school)) +
  coord_cartesian(ylim = yl) +
  labs(x = expression(tau), y = 'E[theta_j|tau,y)', title = title2, color = '') +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))+
  theme(plot.title = element_text(hjust = 0.5, size =rel(2)),
        axis.title=element_text(size=22))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9)) 
plot_condmeans
dev.off()


df_condsds <- as.data.frame(t(tsd)) %>% setNames(LETTERS[1:8]) %>%
  cbind(x = t(tt)) %>% gather(school, p, -x)

pdf(file="eightschool_condsd.pdf", width =10, height =7)
yl <- c(0, 25)
title3 <- 'Conditional posterior standard deviations of effects sd[theta_j|tau,y]'
plot_condsds <- ggplot(data = df_condsds) +
  geom_line(aes(x = x, y = p, color = (school=='A'), group = school)) +
  coord_cartesian(ylim = yl) +
  labs(x = expression(tau), y = 'sd[theta_j|tau,y)', title = title3, color = '') +
  scale_color_manual(values = c('blue','red'), labels = labs1) +
  theme(legend.background = element_blank(), legend.position = c(0.8,0.9))+
  theme(plot.title = element_text(hjust = 0.5, size =rel(2)),
        axis.title=element_text(size=22))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.9)) 
plot_condsds
dev.off()



pdf(file="eightschool_postsummaries.pdf", width =11, height =9)
grid.arrange(plot_margpost, plot_condmeans, plot_condsds)
dev.off()

# hierarchical analysis

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
fit1 <- stan(file = '8schools.stan', data = schools_dat, 
             iter = 1000, chains = 4)
fit2 <- stan(file = '8schools_invgamma.stan', data = schools_dat, 
             iter = 1000, chains = 4)
fit3 <- stan(file = '8schools_halfcauchy.stan', data = schools_dat, 
             iter = 1000, chains = 4)

sims1 <- extract(fit1)
sims2 <- extract(fit2)
sims3 <- extract(fit3)


pdf(file="eightschools_3priors.pdf", width =12, height =5)
par(mfrow=c(1,3), mar =c(5,5,4,1))
hist(sims1$tau, probability = TRUE, xlab=expression(tau), main=
       "Uniform(0,100)", cex.main =3, breaks=30, cex.lab=2.5,
     col="gray", border ="black")
abline(h=0.01, col="red", lwd =3)
hist(sims2$tau, probability = TRUE, xlab=expression(tau), main= " InvGamma(0.01,0.01)" , 
     cex.main =3,  breaks=30, cex.lab=2.5, col="gray", border ="black")
curve(dinvgamma(x^2, 0.01,0.01)*2*x, add=TRUE, col="red", lwd =3)
hist(sims3$tau, probability = TRUE, xlab=expression(tau), main="HalfCauchy(0,2.5)",
     cex.main =3,  breaks=30, cex.lab=2.5, col="gray", border ="black")
curve(dcauchy(x, 0, 2.5), add=TRUE, col="red", lwd =3)
dev.off()


#####################################
## HIEARARCHICAL LOGISTIC REGRESSION
#####################################



# Set up the data for the election88 example

# Load in data for region indicators
# Use "state", an R data file (type ?state from the R command window for info)
#
# Regions:  1=northeast, 2=south, 3=north central, 4=west, 5=d.c.
# We have to insert d.c. (it is the 9th "state" in alphabetical order)

library (arm)
data (state)                  # "state" is an R data file
state.abbr <- c (state.abb[1:8], "DC", state.abb[9:50])
dc <- 9
not.dc <- c(1:8,10:51)
region <- c(3,4,4,3,4,4,1,1,5,3,3,4,4,2,2,2,2,3,3,1,1,1,2,2,3,2,4,2,4,1,1,4,1,3,2,2,3,4,1,1,3,2,3,3,4,1,3,4,1,2,4)

# Load in data from the CBS polls in 1988
# Data are at http://www.stat.columbia.edu/~gelman/arm/examples/election88
library (foreign)
polls <- read.dta ("polls.dta")
attach.all (polls)

# Select just the data from the last survey (#9158)
table (survey)                # look at the survey id's
ok <- survey==9158            # define the condition
polls.subset <- polls[ok,]    # select the subset of interest
attach.all (polls.subset)     # attach the subset
write.table (polls.subset, "polls.subset.dat")

print (polls.subset[1:5,])

# define other data summaries
y <- bush                  # 1 if support bush, 0 if support dukakis
n <- length(y)             # of survey respondents
n.age <- max(age)          # of age categories
n.edu <- max(edu)          # of education categories
n.state <- max(state)      # of states
n.region <- max(region)    # of regions

# compute unweighted and weighted averages for the U.S.
ok <- !is.na(y)                                    # remove the undecideds
cat ("national mean of raw data:", round (mean(y[ok]==1), 3), "\n")
cat ("national weighted mean of raw data:",
     round (sum((weight*y)[ok])/sum(weight[ok]), 3), "\n")

# compute weighted averages for the states
raw.weighted <- rep (NA, n.state)
names (raw.weighted) <- state.abbr
for (i in 1:n.state){
  ok <- !is.na(y) & state==i
  raw.weighted[i] <- sum ((weight*y)[ok])/sum(weight[ok])
}

# load in 1988 election data as a validation check
election88 <- read.dta ("election88.dta")
outcome <- election88$electionresult

# load in 1988 census data
census <- read.dta ("census88.dta")

# also include a measure of previous vote as a state-level predictor
presvote <- read.dta ("presvote.dta")
attach (presvote)
v.prev <- presvote$g76_84pr
not.dc <- c(1:8,10:51)
candidate.effects <- read.table ("candidate_effects.dat", header=T)
v.prev[not.dc] <- v.prev[not.dc] +
  (candidate.effects$X76 + candidate.effects$X80 + candidate.effects$X84)/3
# Data are at http://www.stat.columbia.edu/~gelman/arm/examples/election88

## Multilevel logistic regression: varying-intercept model

M1 <- glmer (y ~ black + female + (1 | state), family=binomial(link="logit"))
display (M1)

M1.rstanarm <- stan_glmer (y ~ black + female + (1 | state), family=binomial(link="logit"))
print(M1.rstanarm)


beta_names <- c(paste0("beta^", c( "black","female")), "gl.intercept")
alpha_names<-c()
for (i in 1:51){
  alpha_names[i] <- paste0(expression(alpha), "[", i,"]")
  }

posterior_M1 <- as.matrix(M1.rstanarm)

pdf("logistic_M1.pdf", width =12, height =7.4)
mcmc_intervals(posterior_M1, regex_pars=c( "black", 
                                          "female",
                                          "(Intercept)", "b"))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.4))+
  scale_y_discrete(labels = ((parse(text= c(beta_names, alpha_names)))))
  #ggtitle("MCMC estimates")+
  #theme(plot.title = element_text(hjust = 0.5, size =rel(2)))
dev.off()

# plot the random effects

pdf("random_effects_log.pdf", height =7, width =11)
int_ord <- sort(coef(M1.rstanarm)$state[,1], index.return=TRUE)$x
ord <- sort(coef(M1.rstanarm)$state[,1], index.return=TRUE)$ix
state.abbr.ord <- state.abbr[ord]
se_ord <- M1.rstanarm$ses[ord]
par(xaxt="n", mfrow=c(1,1), mar = c(5,2,2,1))
plot( 1:49, int_ord, ylim=c(-0.7,1.4), pch=19, bg=2, xlab="States", 
      ylab="Intercepts",  cex.main=1.9, cex.lab=1.9)
for (h in 1:49){
  segments(h, int_ord[h]-se_ord[h], h, int_ord[h]+se_ord[h], col="red")
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5) 
      abs(x - round(x)) < tol
  if (is.wholenumber(h/2)){
    text(h, int_ord[h]+se_ord[h]+0.1, state.abbr.ord[h], cex=1.1)}else{
      text(h, int_ord[h]-se_ord[h]-0.1, state.abbr.ord[h], cex=1.1)
    }
}
dev.off()


# varying-intercept and slope model

M2.rstanarm <- stan_glmer (y ~ black + female + (1+ female | state), 
             family=binomial(link="logit"))
print(M2.rstanarm)

pdf("random_effects_log2.pdf", height =7, width =11.7)
int_ord2 <- sort(coef(M2.rstanarm)$state[,1], index.return=TRUE)$x
ord2 <- sort(coef(M2.rstanarm)$state[,1], index.return=TRUE)$ix
state.abbr.ord2 <- state.abbr[ord2]
state_ind <- grep("b[(Intercept) state:",rownames(as.matrix(M2.rstanarm$ses)),fixed=TRUE)
se_ord2 <- M2.rstanarm$ses[state_ind][ord2]
par(xaxt="n", mfrow=c(1,2))
plot( 1:49, int_ord2, ylim=c(-0.7,1.6), pch=19, bg=2, xlab="States",
      ylab="Intercepts", main =expression(alpha[j]), cex.main =1.9, cex.lab =1.4)
for (h in 1:49){
  segments(h, int_ord2[h]-se_ord2[h], h, 
           int_ord2[h]+se_ord2[h], col="red")
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5) 
      abs(x - round(x)) < tol
  if (is.wholenumber(h/2)){
    text(h, int_ord2[h]+se_ord2[h]+0.1, 
         state.abbr.ord2[h], cex=0.6)}else{
           text(h, int_ord2[h]-se_ord2[h]-0.1, 
                state.abbr.ord2[h], cex=0.6)
         }
}

int_ord2.slope <- sort(coef(M2.rstanarm)$state[,3], index.return=TRUE)$x
ord2.slope <- sort(coef(M2.rstanarm)$state[,3], index.return=TRUE)$ix
state.abbr.ord2.slope <- state.abbr[ord2.slope]

female_ind <- grep("b[female state:",rownames(as.matrix(M2.rstanarm$ses)),fixed=TRUE)


se_ord2.slope <- M2.rstanarm$ses[female_ind][ord2.slope]
plot( 1:49, int_ord2.slope, ylim=c(-0.6,0.5), pch=19, bg=2, 
      xlab="States",
      ylab="Slopes", main =expression(beta[j]), cex.main =1.9, cex.lab =1.4)
for (h in 1:49){
  segments(h, int_ord2.slope[h]-se_ord2.slope[h], h,
           int_ord2.slope[h]+se_ord2.slope[h], col="blue")
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  
      abs(x - round(x)) < tol
  if (is.wholenumber(h/2)){
    text(h, int_ord2.slope[h]+se_ord2.slope[h]+0.1, 
         state.abbr.ord2.slope[h], cex=0.6)}else{
           text(h, int_ord2.slope[h]-se_ord2.slope[h]-0.1, 
                state.abbr.ord2.slope[h], cex=0.6)
         }
}
dev.off()


# Model comparison via looic

install.packages("loo")
library(loo)
lpd1 <- log_lik(M1.rstanarm)
loo1 <- loo(lpd1)
lpd2 <- log_lik(M2.rstanarm)
loo2 <- loo(lpd2)

c(loo1$looic, loo2$looic)

## A fuller model

# set up the predictors
age.edu <- n.edu*(age-1) + edu
region.full <- region[state]
v.prev.full <- v.prev[state]

# fit the model
M2 <- glmer (y ~ black + female + black:female + v.prev.full + (1 | age) + 
              (1 | edu) + (1 | age.edu) + (1 | state) + (1 | region.full), family=binomial(link="logit"))
display (M2)

### Fit the model in rstanarm
M2.rstanarm <- stan_glmer(y ~ black + female + black:female + v.prev.full + (1 | age) + 
                        (1 | edu) + (1 | age.edu) + (1 | state) + (1 | region.full), family=binomial(link="logit"))
## Plot Figure 14.1 
attach.bugs (M2.bugs)
par (mar=c(0,0,0,0))
summ <- M2.bugs$summary[c(2:28),3:7]
labels <- c("female","black","female x black",
            "18-29","30-44","45-64","65+",
            "no h.s.","high school","some college","college grad",
            "18-29 x no h.s.","18-29 x high school","18-29 x some college","18-29 x college grad",
            "30-44 x no h.s.","30-44 x high school","30-44 x some college","30-44 x college grad",
            "45-64 x no h.s.","45-64 x high school","45-64 x some college","45-64 x college grad",
            "65+ x no h.s.","65+ x high school","65+ x some college","65+ x college grad")

pos <- c (1:3, 5:8, 10:13, 15:18, 20:23, 25:28, 30:33)
bottom <- max(pos)+1

rng <- range(summ)
p.rng <- pretty(rng)
a <- -min(p.rng)/(max(p.rng)-min(p.rng))
b <- 1/(max(p.rng)-min(p.rng))
summ.adj <- a + b*summ
plot (c(-.25,1), c(2,-bottom-2), xlab="", ylab="", xaxt="n", yaxt="n",
      type="n", bty="n")
for (i in 1:nrow(summ)){
  text (-.25, -pos[i], labels[i], adj=0, cex=1.1)
  points (summ.adj[i,3], -pos[i], pch=20, cex=1.5)
  lines (summ.adj[i,c(2,4)], rep(-pos[i],2), lwd=4)
  lines (summ.adj[i,c(1,5)], rep(-pos[i],2), lwd=.5)
}
lines (rep(a,2), c(0,-bottom), lwd=.5)
lines (c(0,1), rep(0,2))
lines (c(0,1), rep(-bottom,2))
for (x in p.rng){
  text (a+b*x, 1, x, cex=1.2)
  lines (rep(a+b*x,2), c(0,-.2))
  text (a+b*x, -bottom-1, x, cex=1.2)
  lines (rep(a+b*x,2), -bottom+c(0,.2))
}

## Plot Figure 14.2 

# create linear predictors
attach.bugs (M2.bugs)
linpred <- rep (NA, n)
for (i in 1:n){
  linpred[i] <- mean (b.0 + b.female*female[i] + b.black*black[i] +
                        b.female.black*female[i]*black[i] + a.age[,age[i]] + a.edu[,edu[i]] +
                        a.age.edu[,age[i],edu[i]])
}

# plot the 8 states
par (mfrow=c(2,4))
y.jitter <- y + ifelse (y==0, runif (n, 0, .1), runif (n, -.1, 0))
state.name.all <- c(state.name[1:8], "District of Columbia", state.name[9:50])
for (j in c(2,3,4,8,6,7,5,9)) {
  plot (0, 0, xlim=range(linpred), ylim=c(0,1), yaxs="i", pch=20,
        xlab="linear predictor", ylab="Pr (support Bush)",
        main=state.name.all[j], type="n")
  for (s in 1:20){
    curve (invlogit (a.state[s,j] + x), lwd=.5, add=TRUE, col="gray20")}
  curve (invlogit (median (a.state[,j]) + x), lwd=2, add=TRUE)
  if (sum(state==j)>0) points (linpred[state==j], y.jitter[state==j])
}

## Using the model inferences to estimate avg opinion for each state

# construct the n.sims x 3264 matrix
L <- nrow (census)
y.pred <- array (NA, c(n.sims, L))
for (l in 1:L){
  y.pred[,l] <- invlogit(b.0 + b.female*census$female[l] +
                           b.black*census$black[l] + b.female.black*census$female[l]*census$black[l] +
                           a.age[,census$age[l]] + a.edu[,census$edu[l]] +
                           a.age.edu[,census$age[l],census$edu[l]] + a.state[,census$state[l]])
}

# average over strata within each state
y.pred.state <- array (NA, c(n.sims, n.state))
for (s in 1:n.sims){
  for (j in 1:n.state){
    ok <- census$state==j
    y.pred.state[s,j] <- sum(census$N[ok]*y.pred[s,ok])/sum(census$N[ok])
  }
}

# average over strata within each state
state.pred <- array (NA, c(n.state,3))
for (j in 1:n.state){
  state.pred[j,] <- quantile (y.pred.state[,j], c(.25,.5,.75))
}

## Plot Figure 14.3
attach (M2.bugs$sims.list)  ## ???
region.name <- c("Northeast", "Midwest", "South", "West", "D.C.")
par (mfrow=c(1,4), mar=c(4,4,3,1), oma=c(1,1,2,1))
for (k in 1:4){
  plot (range(v.prev[not.dc]), range(M2.bugs$median$a.state[not.dc]), cex.lab=1.2,
        cex.axis=1.2, cex.main=1.5, ylim=c(-.7,.7), xaxt="n", yaxt="n",
        xlab="R vote in prev elections",
        ylab="regression intercept", pch=20,
        main=region.name[k], type="n")
  axis (1, c(.5,.6,.7), cex.axis=1.2)
  axis (2, c(-.5,0,.5), cex.axis=1.2)
  for (j in (1:n.state)[region==k]){
    lines (rep(v.prev[j],2), quantile(a.state[,j], c(.25,.75)), lwd=.5, col="gray")
    text (v.prev[j], M2.bugs$median$a.state[j], state.abbr[j], cex=1.2)
  }
  curve (median(a.region[,k]) - .7 + median(v.prev)*x, lwd=.5, add=T)
}

## Plot Figure 14.5







