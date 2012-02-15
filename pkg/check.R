rm(list=ls())
setwd("/home/federico/Documents/uni/dottorato_pd_2010/Rcode/simfms/pkg/R")
for (f in list.files()) 
  if (! f %in% c("simfms.R", "scan.tmat.R"))
    source(f)


nsim  = NULL
tmat=trans.cancer.reduced()
clock = "forward"
# Frailty
frailty = list(dist="gamma", par= .5)
nclus=5
csize=2
# Covariates
covs = list(age=function(x) rnorm(x, mean=60, sd=7),
            treat=function(x) rbinom(x, 1, .5))
beta = list(age=-2:2/100, treat=-2:2)
# Marginals
marg  = list(dist="weibull", lambda=1, rho=1:5)
cens  = list(dist="weibull", lambda=c(NED=1, DM=2, LR=3), rho=1, admin= 72)
# Copula
copula= list(name="clayton", par= 1)

source("simfms.R")
data

inTrans = NULL
subjs = 1:5
source("scan.tmat.R")
data

# simfms(tmat=trans.cancer.reduced(),
#        nclus=5, csize=2,
#        marg  = list(dist="weibull",
#                     lambda=1, rho=1:5),
#        frailty=list(dist="gamma", par=.5),
#        covs = list(age=function(x) rnorm(x, mean=60, sd=7),
#                    treat=function(x) rbinom(x, 1, .5)),
#        beta = list(age=-2:2/100, 
#                    treat=-2:2))
