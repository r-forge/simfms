################################################################################
#  Tune simulation parameters for frailty multi-state data                     #
################################################################################
#                                                                              #
#  Looks for appropriate simulation parameters for multi-state data,           #
#  given target values of                                                      #
#   - competing events probabilities and                                       #
#   - medians of uncensored times                                              #
#                                                                              #
#  Its parameters are                                                          #
#   - nsim      : the number of subjects to simulate at each step              #
#   - tmat      : the transitions matrix                                       #
#   - clock     : either 'forward' or 'reset'                                  #
#   - frailty   : the frailty term specifications. A list with components      #
#                 dist: the name of the frailty distribution                   #
#                 par : the frailty parameter(s) value                         #
#   - nclus     : the number of clusters to simulate                           #
#   - csize     : the size(s) of cluster                                       #
#   - covs      : the covariates to simulate. A list with components           #
#                 nameOfCovariate = simulationfunction                         #
#   - beta      : the regression coefficients for covariates. A list of the    #
#                 same length as 'covs', with elements                         #
#                 nameOfCovariate = a vector of the same length as the number  #
#                 of transitions in 'tmat'                                     #
#   - marg      : the marginal baseline hazards. A list with components        #
#                 dist    : the name of the baseline hazard distribution       #
#                            (one value)                                       #
#                 eachpar : initial values of each baseline parameter          #
#                           (either one value or as many as the number         #
#                            of transitions in 'tmat')                         #
#   - cens      : the censoring time distributions. A list with components     #
#                 dist : the name of the censoring distributions (one value)   #
#                 eachpar : each censoring distribution parameter              #
#                           (either one value or as many as the number of      #
#                            possible starting states in 'tmat')               #
#                 admin: the time of administrative censoring                  #
#   - copula    : the copula model. A list with components                     #
#                 name : the name of the copula                                #
#                 par  : the copula parameter                                  #
#   - CEprobs   : the competing events probabilities. A matrix with the same   #
#                 structure as 'tmat', i.e. with NAs at the same places and    #
#                 with the line-specific probabilities associated to possible  #
#                 transitions in other cells                                   #
#   - copula    : the medians of uncensored times of possible events.          #
#                 A matrix with the same structure as 'tmat', i.e. with NAs    #
#                 at the same places and with the line-specific medians        # 
#                 in other cells                                               #
#   - target    : target values for probabilities of competing events and for  #
#                 medians of uncensored times of each transition.              #
#                 A list with elements 'prob' and 'meds', both vector of the   #
#                 same length as the number of transitions in 'tmat'           #
#                                                                              #
#                                                                              #
#   Date: February, 17, 2012                                                   #
#   Last modification on: February, 20, 2012                                   #
################################################################################

# tune.simfms <- function(nsim  = NULL,
#                         tmat  = NULL,
#                         clock = "forward",
#                         # Frailty
#                         frailty = list(dist="gamma",
#                                        par= .5),
#                         nclus = NULL, 
#                         csize = NULL,
#                         # Covariates
#                         covs = NULL,
#                         beta = NULL,
#                         # Marginals
#                         marg  = list(dist="weibull",
#                                      lambda=1, rho=1), 
#                         cens  = list(dist="weibull", 
#                                      lambda=1, rho=1, 
#                                      admin= 72),
#                         # Copula
#                         copula= list(name="clayton",
#                                      par= 1),
#                         # !!! - TARGET VALUES - !!!
#                         target = list(prob = NULL,
#                                       meds = NULL)
#                         ) {
  ### - CONTROLS - #############################################################
  checked <- checks(nsim=nsim,   tmat=tmat,   
                    clock=clock, frailty=frailty,
                    nclus=nclus, csize=csize,
                    covs=covs,   beta=beta,
                    marg=marg,   cens=cens,
                    copula=copula)
  nsim  <- checked$nsim
  nclus <- checked$nclus
  csize <- checked$csize
  clock <- checked$clock
  marg  <- checked$marg
  cens  <- checked$cens
  rm("checked")
  ###################################################### - END of CONTROLS - ###
  
  ### - INITIALIZATION - ########################################################
  initialized <- initialize.fms(nsim  = nsim,  tmat    = tmat, 
                                clock = clock, frailty = frailty, 
                                nclus = nclus, csize   = csize,
                                covs  = covs,  beta    = beta)
  data  <- initialized$data
  eta   <- initialized$eta
  rm("initialized")
  ################################################# - END of INITIALIZATION - ###
  
#   ### - TUNING of SIMULATION PARAMETERS - #######################################
#   # Detailed data for each transition
#   respars <- scan.tmat.tune(pars=NULL,
#                             data=data, atState=NULL, subjs=1:nrow(data),
#                             eta=eta,   tmat=tmat,    clock=clock,
#                             marg=marg, cens=cens,    copula=copula)
#   ########################################### - END of TUNING of PARAMETERS - ###
#   
#   return(respars)
# }