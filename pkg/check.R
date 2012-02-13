setwd("/home/federico/Documents/uni/dottorato_pd_2010/Rcode/simfms/pkg/R")
for (f in list.files()) source(f)

# debug(simulCRcf)
set.seed(90)
simulCRcf(nsim=30,
                trans = trans.cancer.reduced(),
                clock="f",
              # Frailty
                fdist="gamma", 
                ftheta=.5,
                nclus=NULL, 
                csize=NULL,
              # Covariates
                covs=list(
                  Treat=function(x) rbinom(x, 1, .5), 
                  Age=function(x) rnorm(x, mean=60, sd=7)),
                beta=list(Treat=log(c(.3, 1, rep(1.2, 3))), 
                          Age=log(c(.8, .9, rep(1.2, 3)))/10),
              # Marginals
                ctheta=1, 
                #   prev=NULL,
                marg="weibull", 
                  pars=rbind(lambda=rep(1, 5), rho=rep(1, 5)), 
                cens="gompertz",
                  cpars=rbind(lambda=rep(.7,3), rho=rep(1,3)),
                  adcens=72
  )

set.seed(10000000)
simulCRcf(nsim=30,
                trans = trans.cancer(),
                clock="f",
              # Frailty
                fdist="gamma", 
                ftheta=.5,
                nclus=NULL, 
                csize=NULL,
              # Covariates
                covs=list(
                  Treat=function(x) rbinom(x, 1, .5), 
                  Age=function(x) rnorm(x, mean=60, sd=7)),
                beta=list(Treat=log(c( .3, 1  , 1.2, 1  , 1.2,  .3, 1.2, 1.2)), 
                          Age  =log(c( .8,  .9, 1.2,  .9, 1.2,  .8, 1.2, 1.2))/10),
              # Marginals
                ctheta=1, 
                #   prev=NULL,
                marg="weibull", 
                  pars=rbind(lambda=rep(1, 8), rho=rep(1, 8)), 
                cens="gompertz",
                  cpars=rbind(lambda=rep(.7,4), rho=rep(1,4)),
                  adcens=72
  )


set.seed(10800)
# data = 
  simulCRcf(nsim=60,
                trans = trans.cancer.extended(),
                clock="f",
              # Frailty
                fdist="gamma", 
                ftheta=.5,
                nclus=NULL, 
                csize=NULL,
              # Covariates
                covs=list(
                  Treat=function(x) rbinom(x, 1, .5), 
                  Age=function(x) rnorm(x, mean=60, sd=7)),
                beta=list(Treat=log(c( .3, 1  ,  .6, 1.2, 1  , 1.2,  .3, 1.2, 1.2)), 
                          Age  =log(c( .8,  .9,  .85,1.2,  .9, 1.2,  .8, 1.2, 1.2))/10),
              # Marginals
                ctheta=1, 
                #   prev=NULL,
                marg="weibull", 
                  pars=rbind(lambda=rep(1, 9), rho=rep(1, 9)), 
                cens="gompertz",
                  cpars=rbind(lambda=rep(.7,4), rho=rep(1,4)),
                  adcens=72
  )
 data[(data$LRDM==1 &&data$LR==0 && data$DM==0), ]