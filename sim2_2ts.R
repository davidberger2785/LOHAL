
# data structure: O=(L0, A0, L1, A1, L2)
# cc, the coefficient of instrumental variable

gen_data_fct <- function(n, cc=0, abar = c(1,1), tru_ind = FALSE) {

  L0 <- data.frame(L0.1 = runif(n, -1, 1), L0.2 = rbinom(n, 1, 0.50), L0.3 = runif(n, -1, 1),  L0.4 = runif(n))

  if (tru_ind== TRUE) {
    L1 <- data.frame(L1.1 = L0$L0.1*abar[1] + runif(n, 0, 1), L1.2 = rbinom(n, 1, 0.50), L1.3 = runif(n, -1, 1),  L1.4 = runif(n))
    L2 <- rbinom(n, 1, plogis(0.5*L0$L0.1+ 0.5*L0$L0.2 - 0.5*L0$L0.3 + 0.4*abar[1] + 0.3*L1$L1.1 + 0.7*L1$L1.2 -0.6*L1$L1.3 + 0.3*abar[2] ))
    return(mean(L2))

  } else {
    A0 <- rbinom(n, 1, plogis(0.5 * L0$L0.1 - L0$L0.2 + 2 * L0$L0.3 - cc*L0$L0.4))
    L1 <- data.frame(L1.1 = L0$L0.1*A0 + runif(n, 0, 1), L1.2 = rbinom(n, 1, 0.50), L1.3 = runif(n, -1, 1),  L1.4 = runif(n))
    A1 <- rbinom(n, 1, plogis(0.8 * L0$L0.3 + 0.5*A0 + 0.5 * L1$L1.1 - L1$L1.2 + 2 * L1$L1.3 - cc*L1$L1.4))
    L2 <- rbinom(n, 1, plogis(0.5*L0$L0.1+ 0.5*L0$L0.2 - 0.5*L0$L0.3 + 0.4*A0 + 0.3*L1$L1.1 + 0.7*L1$L1.2 -0.6*L1$L1.3 + 0.3*A1))
    return (list(A0=A0, A1=A1, L0=L0, L1=L1, L2=L2))
  }

}


truth <- gen_data_fct(n=1e7, tru_ind = TRUE);
truth # 0.7839227





library(hal9001)
library(SuperLearner)

n=1000
abar=c(1,1)





########### scenario 1: with instrumental variable  L0.4 and L1.4
# da_withI = gen_data_fct(n=1e3, cc = 2.5)
load( "./test/Rda/sim_data.Rda") # da_withI, da_noI
da_withI

cvFolds = 4



#------------ use drinf package

library(drinf)

drinf_withI <- drinf.tmle(
  L0 = da_withI$L0, L1 = da_withI$L1, L2 = da_withI$L2, A0 = da_withI$A0, A1 = da_withI$A1, abar = c(1,1), outcome_family = "binomial",
  SL.Q = "SL.hal9001", SL.g = "SL.hal9001", SL.Qr = "SL.hal9001", SL.gr = "SL.hal9001", flucOrd = c("targetQ2", "targetQ1"),
  universal = FALSE, universalStepSize = 1e-4, return.models = FALSE, verbose = TRUE, maxIter =30, return.ltmle = TRUE, allatonce = FALSE, tolg = 1e-2, tolQ = 1e-2
)

c(drinf_withI$est, drinf_withI$est.ltmle, drinf_withI$se, drinf_withI$se.ltmle)




#---------- obtain hal estimates for Qn and Ohal estimates for Gn
source("./test/ohal_nuisance_2ts_new.R")


fit_nuis <- ohal_nuisance_2ts(L0 = da_withI$L0, L1 = da_withI$L1, L2 = da_withI$L2, A0 = da_withI$A0, A1 = da_withI$A1, abar = abar, V = cvFolds,
                              outcome_family = "binomial", penalty_threshold = 1e-5, penalty_gamma = 1,
                              newdata = NULL, lambda_seq = exp(seq(-5, -13, length=50)), tolg = 1e-2)
head(fit_nuis)
colMeans(fit_nuis)
summary(fit_nuis)


Qn_list = list(Q2n = fit_nuis$Q2n_cv, Q1n = fit_nuis$Q1n_cv)
gn_list = list(g0n = fit_nuis$g0n_cv, g1n = fit_nuis$g1n_cv)

#--------- use hal to run reduced-dimension regression for G1n and G0n
source("./test/redReg_2ts.R")

folds <- rep(1:cvFolds, each = ceiling(n/cvFolds), length.out=n)
gnr_list <- sapply(1:cvFolds, redReg_gr_2ts, gn = gn_list, Qn = Qn_list, A0 = da_withI$A0, A1 = da_withI$A1, L2 = da_withI$L2, folds = folds, abar = abar,
                   verbose = TRUE, tolg = 1e-2, return.models = TRUE, simplify = FALSE)

length(gnr_list) # = cvFolds
# gnr_list[[1]]$g0nr, $g1nr, $h0nr, $h1nr, $hbarnar; gnr_list[[2]]$g0nr, $g1nr, $h0nr, $h1nr, $hbarnar
gnr_list[[4]]$g0nr  # length = n/cvFolds



#---------perform TMLE step using locally least favorable parametric submodels
maxIter = 30
# set up an empty vector of estimates
psin_trace <- rep(NA, maxIter)
se_trace <- rep(NA, maxIter)

#----------------------------------------
# target Q  according to flucOrd
#----------------------------------------
# initialize vectors in Qnstar and gnstar that will hold targeted nuisance parameters by setting equal to the initial values
Qnstar_list <- Qn_list
gnstar_list <- gn_list
full_Qnstar <- Qnstar_list
full_gnstar <- gnstar_list
tmp <-  do.call(Map, c(c, gnr_list))
# names(tmp): "g0nr","g1nr","h0nr","h1nr","hbarnr","g0rmod", "g1rmod","h0rmod","h1rmod","hbarrmod"
# length(tmp$g0nr) = n

full_gnr <- sapply(unique(names(tmp)), function(x) unname(unlist(tmp[names(tmp)==x])), simplify=FALSE)

est.naive <- mean(full_Qnstar$Q1n); est.naive

#-----------------------------------
# targeting loop for dr inference
#-----------------------------------

flucOrd = c("targetQ2","targetQ1")
meanif.dr <- Inf
n_minus_12 <- length(L2)^(-1/2) # 0.03162278
sqrt_n_max_block <- n_max_block <- sqrt_n_norm_block <- n_norm_block <- FALSE
sqrt_n_max_iter <- sqrt_n_norm_iter <- n_norm_iter <- n_max_iter <- maxIter
tolIF=1/(length(L2))
tolQ=tolg=1e-2
verbose = TRUE

iter <- 0
source("./test/targetQ_2ts.R")

while((max(abs(meanif.dr)) > tolIF | sqrt(sum(meanif.dr^2)) > tolIF) & iter < maxIter){

  iter <- iter + 1

  #--------------------
  # target Q
  #--------------------
  in_ct <- 0  ### ? what does it represent?
  for(ff in flucOrd){

    in_ct <- in_ct + 1
    if(verbose){ cat("Calling ", ff, " for targeting step. \n") }

    flucOut <- do.call(ff, args = list( A0 = da_withI$A0, A1 = da_withI$A1, L2 = da_withI$L2, Qn = full_Qnstar, gn = full_gnstar, gnr = full_gnr,
                                        tolg = tolg, tolQ = tolQ, abar = abar, return.models = FALSE))

    flucOutNames <- names(flucOut) ###

    # look for what names are in function output and assign values accordingly
    if("Q2nstar" %in% flucOutNames){
      if(verbose) cat("Q2n was targeted by ", ff,". \n")
      Qnstar_list$Q2n = as.numeric(flucOut$Q2nstar)
      full_Qnstar <- Qnstar_list
    }
    if("Q1nstar" %in% flucOutNames){
      if(verbose) cat("Q1n was targeted by ", ff,". \n")
      Qnstar_list$Q1n = as.numeric(flucOut$Q1nstar)
      full_Qnstar <- Qnstar_list
    }

  }

  #-------------------------
  # evaluate IF
  #-------------------------
  source("./test/evaluateIF_2ts.R")
  if.dr <- evaluateIF( A0 = da_withI$A0, A1 = da_withI$A1, L2 = da_withI$L2, Q2n = full_Qnstar$Q2n, Q1n = full_Qnstar$Q1n, g1n = full_gnstar$g1n, g0n = full_gnstar$g0n,
                       g0nr = full_gnr$g0nr, g1nr = full_gnr$g1nr, h0nr = full_gnr$h0nr, h1nr = full_gnr$h1nr, hbarnr = full_gnr$hbarnr, abar = abar )

  meanif.dr <- c(
    # original terms
    t(matrix(c(1,1,1))) %*% colMeans(Reduce("cbind",if.dr[1:3])),
    # extra terms targeting Q's
    t(matrix(c(1,1,1))) %*% colMeans(Reduce("cbind",if.dr[4:6]))
  )

  if(max(abs(meanif.dr)) < n_minus_12 & !sqrt_n_max_block){
    sqrt_n_max_iter <- iter
    sqrt_n_max_block <- TRUE # so it doesn't reset
  }

  if(sqrt(sum(meanif.dr^2)) < n_minus_12 & !sqrt_n_norm_block){
    sqrt_n_norm_iter <- iter
    sqrt_n_norm_block <- TRUE
  }

  if(sqrt(sum(meanif.dr^2)) < tolIF & !n_norm_block){
    n_norm_iter <- iter
    n_norm_block <- TRUE
  }

  if(max(abs(meanif.dr)) < tolIF & !n_max_block){
    n_max_iter <- iter
    n_max_block <- TRUE
  }

  # update iteration
  psin_trace[iter] <- mean(full_Qnstar$Q1n)
  se_trace[iter] <- sqrt(var( if.dr$Dstar0 + if.dr$Dstar1 + if.dr$Dstar2 - if.dr$DQ2.g1 - if.dr$DQ2.g0 - if.dr$DQ1.g0 )/length(A0))

  cat("\n mean ic = ", meanif.dr, "\n")
}


psin11 <- mean(full_Qnstar$Q1n)
se11 <- se_trace[iter]
c(psin11, se11)



####
res_withI = data.frame(est = c(psin11, drinf_withI$est, drinf_withI$est.ltmle), se = c(se11, drinf_withI$se, drinf_withI$se.ltmle))
rownames(res_withI) = c("OHAL", "Benk_dr", "Benk_ltmle")
res_withI$bias = truth-res_withI$est
res_withI
truth  # 0.7839227
#               est         se        bias
# OHAL       0.7934607 0.06368574 -0.00953802
# Benk_dr    0.8062385 0.04481287 -0.02231580
# Benk_ltmle 0.8119550 0.03796249 -0.02803227



########### scenario 2: no instrumental variable L0.4 and L1.4

# da_noI = gen_data_fct(n=1e3)


da_noI

#--------- use drinf package

library(drinf)

drinf_noI <- drinf.tmle(
  L0 = da_noI$L0, L1 = da_noI$L1, L2 = da_noI$L2, A0 = da_noI$A0, A1 = da_noI$A1, abar = c(1,1),
  SL.Q = "SL.hal9001", SL.g = "SL.hal9001", SL.Qr = "SL.hal9001", SL.gr = "SL.hal9001",
  flucOrd = c("targetQ2","targetQ1"),
  outcome_family = "binomial",
  universal = FALSE, universalStepSize = 1e-4, return.models = FALSE,
  verbose = TRUE, maxIter = 10,
  return.ltmle = TRUE, allatonce = FALSE,
  tolg = 1e-2, tolQ = 1e-2
)

c(drinf_noI$est, drinf_noI$est.ltmle, drinf_noI$se, drinf_noI$se.ltmle)



#---------- obtain hal estimates for Qn and Ohal estimates for Gn
source("./test/ohal_nuisance_2ts_new.R")
cvFolds = 4

fit_nuis_noI <- ohal_nuisance_2ts(L0 = da_noI$L0, L1 = da_noI$L1, L2 = da_noI$L2, A0 = da_noI$A0, A1 = da_noI$A1, abar = abar, V = cvFolds,
                              outcome_family = "binomial", penalty_threshold = 1e-5, penalty_gamma = 1,
                              newdata = NULL, lambda_seq = exp(seq(-5, -13, length=50)), tolg = 1e-2)
head(fit_nuis_noI)
colMeans(fit_nuis_noI)
summary(fit_nuis_noI)


Qn_list = list(Q2n = fit_nuis_noI$Q2n_cv, Q1n = fit_nuis_noI$Q1n_cv)
gn_list = list(g0n = fit_nuis_noI$g0n_cv, g1n = fit_nuis_noI$g1n_cv)

#--------- use hal to run reduced-dimension regression for G1n and G0n
source("./test/redReg_2ts.R")

folds <- rep(1:cvFolds, each = ceiling(n/cvFolds), length.out=n)
gnr_list <- sapply(1:cvFolds, redReg_gr_2ts, gn = gn_list, Qn = Qn_list, A0 = da_noI$A0, A1 = da_noI$A1, L2 = da_noI$L2, folds = folds, abar = abar,
                   verbose = TRUE, tolg = 1e-2, return.models = TRUE, simplify = FALSE)

length(gnr_list) # = cvFolds
# gnr_list[[1]]$g0nr, $g1nr, $h0nr, $h1nr, $hbarnar; gnr_list[[2]]$g0nr, $g1nr, $h0nr, $h1nr, $hbarnar
gnr_list[[4]]$g1nr  # length = n/cvFolds



#---------perform TMLE step using locally least favorable parametric submodels
maxIter = 30
# set up an empty vector of estimates
psin_trace <- rep(NA, maxIter)
se_trace <- rep(NA, maxIter)

#----------------------------------------
# target Q  according to flucOrd
#----------------------------------------
# initialize vectors in Qnstar and gnstar that will hold targeted nuisance parameters by setting equal to the initial values
Qnstar_list <- Qn_list
gnstar_list <- gn_list
full_Qnstar <- Qnstar_list
full_gnstar <- gnstar_list
tmp <-  do.call(Map, c(c, gnr_list))


full_gnr <- sapply(unique(names(tmp)), function(x) unname(unlist(tmp[names(tmp)==x])), simplify=FALSE)

est.naive <- mean(full_Qnstar$Q1n); est.naive

#-----------------------------------
# targeting loop for dr inference
#-----------------------------------

flucOrd = c("targetQ2","targetQ1")
meanif.dr <- Inf
n_minus_12 <- length(L2)^(-1/2) # 0.03162278
sqrt_n_max_block <- n_max_block <- sqrt_n_norm_block <- n_norm_block <- FALSE
sqrt_n_max_iter <- sqrt_n_norm_iter <- n_norm_iter <- n_max_iter <- maxIter
tolIF=1/(length(L2))
tolQ=tolg=1e-2
verbose = TRUE

iter <- 0
source("./test/targetQ_2ts.R")

while((max(abs(meanif.dr)) > tolIF | sqrt(sum(meanif.dr^2)) > tolIF) & iter < maxIter){

  iter <- iter + 1

  #--------------------
  # target Q
  #--------------------
  in_ct <- 0  ### ? what does it represent?
  for(ff in flucOrd){

    in_ct <- in_ct + 1
    if(verbose){ cat("Calling ", ff, " for targeting step. \n") }

    flucOut <- do.call(ff, args = list( A0 = da_noI$A0, A1 = da_noI$A1, L2 = da_noI$L2, Qn = full_Qnstar, gn = full_gnstar, gnr = full_gnr,
                                        tolg = tolg, tolQ = tolQ, abar = abar, return.models = FALSE))

    flucOutNames <- names(flucOut)

    # look for what names are in function output and assign values accordingly
    if("Q2nstar" %in% flucOutNames){
      if(verbose) cat("Q2n was targeted by ", ff,". \n")
      Qnstar_list$Q2n = as.numeric(flucOut$Q2nstar)
      full_Qnstar <- Qnstar_list
    }
    if("Q1nstar" %in% flucOutNames){
      if(verbose) cat("Q1n was targeted by ", ff,". \n")
      Qnstar_list$Q1n = as.numeric(flucOut$Q1nstar)
      full_Qnstar <- Qnstar_list
    }

  }

  #-------------------------
  # evaluate IF
  #-------------------------
  source("./test/evaluateIF_2ts.R")
  if.dr <- evaluateIF( A0 = da_noI$A0, A1 = da_noI$A1, L2 = da_noI$L2, Q2n = full_Qnstar$Q2n, Q1n = full_Qnstar$Q1n, g1n = full_gnstar$g1n, g0n = full_gnstar$g0n,
                       g0nr = full_gnr$g0nr, g1nr = full_gnr$g1nr, h0nr = full_gnr$h0nr, h1nr = full_gnr$h1nr, hbarnr = full_gnr$hbarnr, abar = abar )

  meanif.dr <- c(
    # original terms
    t(matrix(c(1,1,1))) %*% colMeans(Reduce("cbind",if.dr[1:3])),
    # extra terms targeting Q's
    t(matrix(c(1,1,1))) %*% colMeans(Reduce("cbind",if.dr[4:6]))
  )

  if(max(abs(meanif.dr)) < n_minus_12 & !sqrt_n_max_block){
    sqrt_n_max_iter <- iter
    sqrt_n_max_block <- TRUE # so it doesn't reset
  }

  if(sqrt(sum(meanif.dr^2)) < n_minus_12 & !sqrt_n_norm_block){
    sqrt_n_norm_iter <- iter
    sqrt_n_norm_block <- TRUE
  }

  if(sqrt(sum(meanif.dr^2)) < tolIF & !n_norm_block){
    n_norm_iter <- iter
    n_norm_block <- TRUE
  }

  if(max(abs(meanif.dr)) < tolIF & !n_max_block){
    n_max_iter <- iter
    n_max_block <- TRUE
  }

  # update iteration
  psin_trace[iter] <- mean(full_Qnstar$Q1n)
  se_trace[iter] <- sqrt(var( if.dr$Dstar0 + if.dr$Dstar1 + if.dr$Dstar2 - if.dr$DQ2.g1 - if.dr$DQ2.g0 - if.dr$DQ1.g0 )/length(A0))

  cat("\n mean ic = ", meanif.dr, "\n")
}


psin11_noI <- mean(full_Qnstar$Q1n)
se11_noI <- se_trace[iter]
c(psin11_noI, se11_noI)



###
res_noI = data.frame(est = c(psin11_noI, drinf_noI$est, drinf_noI$est.ltmle), se = c(se11_noI, drinf_noI$se, drinf_noI$se.ltmle))
rownames(res_noI) = c("OHAL", "Benk_dr", "Benk_ltmle")
res_noI$bias = truth-res_noI$est
res_noI
#               est         se       bias
# OHAL       0.7707707 0.03022389 0.01315197
# Benk_dr    0.7686500 0.02880545 0.01527272
# Benk_ltmle 0.7664925 0.02741904 0.01743021














