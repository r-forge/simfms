rm(list=ls())
setwd("/home/federico/Documents/uni/dottorato_pd_2010/Rcode/simfms/pkg/R")
for (f in list.files()) source(f)


simfms(tmat=trans.cancer.reduced(),
       nclus=5, csize=2,
       marg  = list(dist="weibull",
                    lambda=1, rho=1:5),
       frailty=list(dist="gamma", par=.5),
       covs = list(age=function(x) rnorm(x, mean=60, sd=7),
                   treat=function(x) rbinom(x, 1, .5)),
       beta = list(age=-2:2/100, 
                   treat=-2:2))
