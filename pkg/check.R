rm(list=ls())
setwd("/home/federico/Documents/uni/dottorato_pd_2010/Rcode/simfms/pkg/R")
for (f in list.files()) 
#    if (! f %in% c("simfms.R" #, "scan.tmat.R"
#                   ))
    source(f)

simfms(nsim  = NULL,
       tmat  = trans.cancer.reduced(),
       clock = "forward",
       frailty = list(dist="gamma",
                      par= .5),
       nclus = 5,
       csize = 2,
       covs = list(age=function(x) rnorm(x, mean=60, sd=7),
                   treat=function(x) rbinom(x, 1, .5)),
       beta = list(age=rep(.02,5), treat=rep(2,5)),
       marg  = list(dist="weibull",
                    lambda=1, rho=1), 
       cens  = list(dist="weibull", 
                    lambda=1, rho=1, 
                    admin= .5),
       copula= list(name="clayton",
                    par= 1)
       )


simfms(nsim  = NULL,
       tmat  = trans.cancer(),
       clock = "forward",
       frailty = list(dist="gamma",
                      par= .5),
       nclus = 5,
       csize = 2,
       covs = list(age=function(x) rnorm(x, mean=60, sd=7),
                   treat=function(x) rbinom(x, 1, .5)),
       beta = list(age=rep(.02, 8), treat=rep(2, 8)),
       marg  = list(dist="weibull",
                    lambda=1, rho=1), 
       cens  = list(dist="weibull", 
                    lambda=1, rho=1, 
                    admin= 72),
       copula= list(name="clayton",
                    par= 1)
       )


simfms(nsim  = NULL,
       tmat  = trans.cancer.extended(),
       clock = "forward",
       frailty = list(dist="gamma",
                      par= .5),
       nclus = 5,
       csize = 2,
       covs = list(age=function(x) rnorm(x, mean=60, sd=7),
                   treat=function(x) rbinom(x, 1, .5)),
       beta = list(age=rep(.02, 9), treat=rep(2, 9)),
       marg  = list(dist="weibull",
                    lambda=1, rho=1), 
       cens  = list(dist="weibull", 
                    lambda=1, rho=1, 
                    admin= 72),
       copula= list(name="clayton",
                    par= 1)
       )


nsim  = NULL
tmat=trans.cancer()
clock = "forward"
frailty = list(dist="gamma", par= .5)
nclus=4
csize=3
covs = list(age=function(x) rnorm(x, mean=60, sd=7),
            treat=function(x) rbinom(x, 1, .5))
beta = list(age=rep(.02, 8), treat=rep(2, 8))
marg  = list(dist="weibull", lambda=1, rho=1)
cens  = list(dist="weibull", lambda=.2, rho=1, admin= 72)
copula= list(name="clayton", par= 1)

source("simfms.R")
data

asd = data
table(asd$tr1.status, asd$tr4.status)
table(asd$tr2.status, asd$tr5.status)
table(asd$tr3.status, asd$tr4.status)
table(asd$tr3.status, asd$tr5.status)

data 
tmat





# simfms(tmat=trans.cancer.reduced(),
#        nclus=5, csize=2,
#        marg  = list(dist="weibull",
#                     lambda=1, rho=1:5),
#        frailty=list(dist="gamma", par=.5),
#        covs = list(age=function(x) rnorm(x, mean=60, sd=7),
#                    treat=function(x) rbinom(x, 1, .5)),
#        beta = list(age=-2:2/100, 
#                    treat=-2:2))
