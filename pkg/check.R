rm(list=ls())
setwd("/home/federico/Documents/uni/dottorato_pd_2010/Rcode/simfms/pkg/R")
for (f in list.files()) 
   if (! f %in% c("tune.simfms.R" , "scan.tmat.tune.R"
                  ))
    source(f)

# simfms(nsim  = NULL,
#        tmat  = trans.cancer.reduced(),
#        clock = "forward",
#        frailty = list(dist="gamma",
#                       par= .5),
#        nclus = 5,
#        csize = 2,
#        covs = list(age=function(x) rnorm(x, mean=60, sd=7),
#                    treat=function(x) rbinom(x, 1, .5)),
#        beta = list(age=rep(.02,5), treat=rep(2,5)),
#        marg  = list(dist="weibull",
#                     lambda=1, rho=1), 
#        cens  = list(dist="weibull", 
#                     lambda=1, rho=1, 
#                     admin= .5),
#        copula= list(name="clayton",
#                     par= 1)
#        )
# 
# 
# simfms(nsim  = NULL,
#        tmat  = trans.cancer(),
#        clock = "forward",
#        frailty = list(dist="gamma",
#                       par= .5),
#        nclus = 5,
#        csize = 2,
#        covs = list(age=function(x) rnorm(x, mean=60, sd=7),
#                    treat=function(x) rbinom(x, 1, .5)),
#        beta = list(age=rep(.02, 8), treat=rep(2, 8)),
#        marg  = list(dist="weibull",
#                     lambda=1, rho=1), 
#        cens  = list(dist="weibull", 
#                     lambda=1, rho=1, 
#                     admin= 72),
#        copula= list(name="clayton",
#                     par= 1)
#        )
# 
# 
# simfms(nsim  = NULL,
#        tmat  = trans.cancer.extended(),
#        clock = "forward",
#        frailty = list(dist="gamma",
#                       par= .5),
#        nclus = 5,
#        csize = 2,
#        covs = list(age=function(x) rnorm(x, mean=60, sd=7),
#                    treat=function(x) rbinom(x, 1, .5)),
#        beta = list(age=rep(.02, 9), treat=rep(2, 9)),
#        marg  = list(dist="weibull",
#                     lambda=1, rho=1), 
#        cens  = list(dist="weibull", 
#                     lambda=1, rho=1, 
#                     admin= 72),
#        copula= list(name="clayton",
#                     par= 1)
#        )
# 
# 
# 
# tune.simfms(nsim  = NULL,
#             tmat  = trans.cancer.reduced(),
#             clock = "forward",
#             frailty = list(dist="gamma",
#                            par= .5),
#             nclus = 10, 
#             csize = 2,
#             # Covariates
#             covs = list(age=function(x) rnorm(x, mean=60, sd=7),
#                         treat=function(x) rbinom(x, 1, .5)),
#             beta = list(age=c(log(.8), log(.9), rep(log(1.2), 3))/10, 
#                         treat=c(log(.3), 0 ,rep(1.2, 3))),
#             # Marginals
#             marg  = list(dist="weibull",
#                          lambda=1, rho=1), 
#             cens  = list(dist="weibull", 
#                          lambda=1, rho=1, 
#                          admin= 72),
#             # Copula
#             copula= list(name="clayton",
#                          par= 1),
#             # !!! - TARGET VALUES - !!!
#             target = list(prob = rbind(
#               NED= c(NED=NA, LR=0.34, DM=0.09, De=0.07),
#               LR=c(NA, NA, NA, 0.47),
#               DM=c(NA, NA, NA, 0.95),
#               De=NA),
#                           meds = rbind(
#                             NED= c(NED=NA, LR=6, DM=10, De=3),
#                             LR=c(NA, NA, NA, 3.25),
#                             DM=c(NA, NA, NA, 0.5),
#                             De=NA))
#             )



nsim  = NULL
       tmat  = trans.cancer.reduced()
       clock = "forward"
       frailty = list(dist="gamma",
                      par= .5)
       nclus = 5
       csize = 2
covs = list(age=function(x) rnorm(x, mean=60, sd=7),
            treat=function(x) rbinom(x, 1, .5))
beta = list(age=rep(.02, 5), treat=rep(2, 5))
marg  = list(dist="weibull", lambda=1, rho=1)
cens  = list(dist="weibull", lambda=.2, rho=1, admin= .72)
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
