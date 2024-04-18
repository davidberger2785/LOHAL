# small test


# fix this


library(SuperLearner)
library(glmnet)
library(gam)
library(lqa)
library(tmle)
library(drtmle)
library(hal9001)

# this setting contains instrumental variable W4
dat_gen_instru <- function(n, trueInd = FALSE, seed=sample(1:100000,1)) {

  set.seed(seed)
  if ( !trueInd) {
    W <- data.frame(W1 = runif(n, -1, 1), W2 = rbinom(n, 1, 0.50), W3 = runif(n, -1, 1),  W4 = runif(n))
    A <- rbinom(n, 1, plogis(0.5 * W$W1 - W$W3 + 2 * W$W3 * W$W2 - 2.5*W$W4))
    Qbar0AW <- plogis(-2 * W$W1*as.numeric(W$W1 > -1/2) - W$W3 + 2 * W$W3 * W$W2 + A)
    Y <- rbinom(n, 1, Qbar0AW)
    dat = data.frame(cbind(W, A, Y))

    return(dat)

  } else {

    big_n <- 1e6
    big_W <- data.frame(W1 = runif(big_n, -1, 1), W2 = rbinom(big_n, 1, 0.50), W3 = runif(big_n, -1, 1))
    big_Qbar01W <- plogis(-2 * big_W$W1*as.numeric(big_W$W1 > -1/2) - big_W$W3 + 2 * big_W$W3 * big_W$W2 + 1)
    big_Qbar00W <- plogis(-2 * big_W$W1*as.numeric(big_W$W1 > -1/2) - big_W$W3 + 2 * big_W$W3 * big_W$W2 + 0)
    truth <- mean(big_Qbar01W - big_Qbar00W)

    return(c(mean(big_Qbar01W), mean(big_Qbar00W), truth))
  }
}

truth = dat_gen_instru(n=1000, trueInd = TRUE)
# 0.6288655 0.4251338 0.2037317


# simulation for one time step from Ju's paper

n <- 1000
datt = dat_gen_instru(n=n)
W = datt[c("W1", "W2", "W3", "W4")]
A = datt$A
Y = datt$Y
# regression formulas for iptw and gcomp
correct_gform <- "A ~ W1 + W3 + I(W3*W2) + W4"
correct_Qform <- "Y ~ I(W1*as.numeric(W1 > 0)) + W3 + I(W3*W2) + A"
sim_fam <- binomial()

# call Ju`s ohal_nuisance function to get nuisance estimates
source("ohal_functions.R")
fit <- ohal_nuisance(W = W, A = A, Y = Y, V = 5,
                     outcome_family = sim_fam$family,
                     lambda_seq = exp(seq(-1, -13, length=100)))
head(fit)
#         QAW       Q1W       Q0W  G1W_OHAL   G1W_HAL     cv_QAW    cv_Q1W     cv_Q0W cv_G1W_OHAL cv_G1W_HAL fold_vec
# 1 0.5459684 0.7681405 0.5459684 0.2668404 0.1522684 0.4657354 0.7079995 0.4657354   0.1730536  0.1611381        1
# 2 0.4781575 0.7162674 0.4781575 0.2091937 0.1920035 0.4896961 0.7274533 0.4896961   0.1543022  0.1794550        1
# 3 0.7940551 0.7940551 0.5832422 0.5434328 0.4972074 0.7993286 0.7993286 0.5888327   0.4853740  0.3889695        1
# 4 0.5243819 0.7523252 0.5243819 0.1779877 0.1427289 0.4515866 0.6960798 0.4515866   0.1454095  0.1569921        1
# 5 0.3297739 0.5754784 0.3297739 0.1370402 0.2301067 0.3255393 0.5731051 0.3255393   0.1164464  0.2360555        1
# 6 0.7127190 0.7127190 0.4738189 0.2091937 0.2709809 0.7052025 0.7052025 0.4623800   0.1543022  0.1759232        1

# call drtmle with OHAL
drtmle_ohal_fit <- drtmle(W = W, A = A, Y = Y,
                          a_0 = c(0,1), tolg = 1e-4,
                          Qn = list(fit$Q0W, fit$Q1W),
                          gn = list(1 - fit$G1W_OHAL, fit$G1W_OHAL),
                          SL_gr = "SL.hal9001",
                          guard = "g")
# $est
#
# 0 0.4263246  # mean(drtmle_ohal_fit$nuisance_drtmle$QnStar[[1]])
# 1 0.6917282  # mean(drtmle_ohal_fit$nuisance_drtmle$QnStar[[2]])
#
# $cov
# 0            1
# 0 3.142364e-04 3.552782e-05
# 1 3.552782e-05 9.471821e-04

# get drtmle for OHAL ATE
drtmle_ohal_ci <- ci(drtmle_ohal_fit, contrast = c(-1, 1))  # default of est value is "drtmle"
se_drtmle <- (drtmle_ohal_ci$drtmle[3] - drtmle_ohal_ci$drtmle[2]) / 2 / 1.96




#########################################################################
## use functions from drinf.tmle

a0 = 1  # estimate the fixed marginal expectation


Qn = list(Q0n = fit$Q0W, Q1n = fit$Q1W)
gn = list(g0n = 1 - fit$G1W_OHAL, g1n = fit$G1W_OHAL)
n <- length(Y)
cvFolds = 1 # only use one cvfolds so far
folds <- rep(1:cvFolds, each = ceiling(n/cvFolds), length.out=n)
SL.Qr = SL.gr = "SL.hal9001"
verbose = TRUE  # default value in drinf.tmle
return.models = FALSE  # default value in drinf.tmle
tolg=1e-8; tolQ=1e-8  # default value in drinf.tmle
tolg=0.01 # use the default value in drtmle

##----------------------------------------
## run reduced-dimension regression of g
source("sim1_functions.R")
Qnr.gnr_list <- sapply(1:cvFolds, redReg1,
                       gn = gn, Qn = Qn, A = A, Y = Y, folds = folds, a0 = a0,
                       verbose = verbose, tolg = tolg, SL.Qr = SL.Qr, SL.gr = SL.gr,
                       return.models = return.models, simplify = FALSE)



#----------------------------------------
# targeting according to locally least favorable models

#----------------------------------------
# target Q according to flucOrd
#----------------------------------------
# initialize vectors in Qnstar and gnstar that will hold
# targeted nuisance parameters by setting equal to the initial values
full_Qnstar <- Qnstar_list<- Qn
full_gnstar <- gnstar_list<- gn
tmp <- do.call(Map, c(c, Qnr.gnr_list))
full_Qnr.gnr <- vector(mode = "list")
full_Qnr.gnr$Qnr <- sapply(unique(names(tmp$Qnr)), function(x) unname(unlist(tmp$Qnr[names(tmp$Qnr)==x])), simplify=FALSE)
full_Qnr.gnr$gnr <- sapply(unique(names(tmp$gnr)), function(x) unname(unlist(tmp$gnr[names(tmp$gnr)==x])), simplify=FALSE)

# est.naive <- mean(full_Qnstar$Q1n)
# 0.6308019
# est.naive <- mean(full_Qnstar$Q0n)

#-----------------------------------
# targeting loop for dr inference
#-----------------------------------
iter <- 0
meanif.dr <- Inf
n_minus_12 <- length(Y)^(-1/2)
maxIter = 20  # default value in drinf.tmle # default value in drtmle: maxIter = 3
# set up an empty vector of estimates
psin_trace <- rep(NA, maxIter)
se_trace <- rep(NA, maxIter)
tolIF=1/(length(Y))
sqrt_n_max_block <- n_max_block <- FALSE
sqrt_n_norm_block <- n_norm_block <- FALSE
sqrt_n_max_iter <- sqrt_n_norm_iter <- n_norm_iter <- n_max_iter <- maxIter
flucOrd = c("targetQ")  # only consider for targeting Q now


while((max(abs(meanif.dr)) > tolIF | sqrt(sum(meanif.dr^2)) > tolIF) & iter < maxIter){

  iter <- iter + 1
  #--------------------
  # target Q
  #--------------------
  in_ct <- 0
  for(ff in flucOrd){
    in_ct <- in_ct + 1
    if(verbose){
      cat("Calling ", ff, " for targeting step. \n")
    }
    # call targetQ(), for Qnr.gnr only need reduce gnr, Qnr.gnr$gnr
    flucOut <- do.call(ff, args = list( A=A, Y = Y, Qn = full_Qnstar, gn = full_gnstar, Qnr.gnr = full_Qnr.gnr, tolg = tolg, tolQ = tolQ, a0 = a0,
                                        return.models = return.models, SL.Qr = SL.Qr, SL.gr = SL.gr, verbose = verbose))  ### output: $Q2nstar (n values); $flucmod
    flucOutNames <- names(flucOut) ### "Q1nstar" "flucmod"


    # look for what names are in function output and assign values accordingly
    if("Qnstar" %in% flucOutNames){
      if(verbose) cat("Qn was targeted by ", ff,". \n")
      split_Qnstar <- split(flucOut$Qnstar, folds)
      if (a0==1) Qnstar_list$Q1n <- as.numeric(split_Qnstar[[1]]  ) else if (a0==0) Qnstar_list$Q0n <- as.numeric(split_Qnstar[[1]]  )
      # Qnstar_list <- mapply(x = Qnstar_list, values = split_Q1nstar[[1]], FUN = list_replace, MoreArgs = list(list = "Q1n"), SIMPLIFY = FALSE)
      full_Qnstar <- Qnstar_list
    }

  }

  #-------------------------
  # evaluate IF
  #-------------------------
  if (a0==1) {Qestn = full_Qnstar$Q1n; gestn = full_gnstar$g1n} else if (a0==0) {Qestn = full_Qnstar$Q0n; gestn = full_gnstar$g0n}

  if.dr <- evaluateIF1( A = A, Y = Y, Q1n = Qestn, g1n = gestn, g1nr = full_Qnr.gnr$gnr$g1nr, g2nr = full_Qnr.gnr$gnr$g2nr, a0 = a0)

  meanif.dr <- c(
    # original terms
    t(matrix(c(1,1))) %*% colMeans(Reduce("cbind",if.dr[1:2])),
    # extra terms targeting Q's
    mean(Reduce("cbind",if.dr[3]))
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
  se_trace[iter] <- sqrt(var( (if.dr$Dstar0 + if.dr$Dstar1 - if.dr$DQ1.g1) )/length(A))

  cat("\n mean ic = ", meanif.dr, "\n")
} # end of while loop




#------------------------------------------
# evaluate parameter and compute std err
#------------------------------------------
# if set a0=1
psin1 <- mean(full_Qnstar$Q1n); if.dr1 = if.dr
se1 <- sqrt(var( (if.dr1$Dstar0 + if.dr1$Dstar1 - if.dr1$DQ1.g1) )/length(A))
esta1 = c(psin1, psin1 +1.96*c(-1,1)*se1)

# # if set a0=0
# psin0 <- mean(full_Qnstar$Q0n); if.dr0 = if.dr
# se0 <- sqrt(var( (if.dr0$Dstar0 + if.dr0$Dstar1 - if.dr0$DQ1.g1) )/length(A))
# esta0 = c(psin0, psin0 +1.96*c(-1,1)*se0)


# results using modified drinf.tmle
res = rbind(esta0, esta1)
colnames(res) = c("est", "cil", "ciu")
#           est       cil       ciu
# esta0 0.3967765 0.3616956 0.4318574
# esta1 0.6825548 0.6201867 0.7449228



# results using drtmle
ci(drtmle_ohal_fit)$drtmle
#       est       cil       ciu
# 0 0.3967846 0.3616750 0.4318942
# 1 0.6825548 0.6201878 0.7449217






########### G-com ###########

gcomp <- function(W, A, Y, formula = "Y ~ .", family = binomail()){
  out_fit <- glm(as.formula(formula), data = data.frame(W, A = A, Y), family = family)
  Q1hat <- predict(out_fit, newdata = data.frame(W, A = 1, Y), type = "response")
  Q0hat <- predict(out_fit, newdata = data.frame(W, A = 0, Y), type = "response")
  return(c(mean(Q1hat), mean(Q0hat), mean(Q1hat - Q0hat)))
}

gcomp_one_boot <- function(W, A, Y, formula = "Y ~ .", family = binomial()){
  resamp_idx <- sample(1:length(Y), replace = TRUE)
  return(gcomp(W[resamp_idx, , drop = FALSE], A = A[resamp_idx], Y = Y[resamp_idx],
               formula = formula, family = family))
}

gcomp_boot_ci <- function(W, A, Y, formula = "Y ~ .", nboot = 5e2, family = binomial()){
  ate_star <- replicate(nboot, gcomp_one_boot(W, A, Y, formula, family = family))
  ci = sapply(1:3, function(x) as.numeric(quantile(ate_star[x,], p = c(0.025, 0.975))))
  colnames(ci) = c("E[Q1hat]", "E[Q0hat]", "ATE")
  return(ci) # 1st col for Q1hat, 2nd col for Q0hat, 3rd col for ATE
  # return(as.numeric(quantile(ate_star, p = c(0.025, 0.975))))
}


est_gcomp_right <- gcomp(W, A, Y, formula = correct_Qform, family = sim_fam)
ci_gcomp_right <- gcomp_boot_ci(W, A, Y, formula = correct_Qform, family = sim_fam)
gcomp_right <- list(gcomp_est = est_gcomp_right, gcomp_ci = ci_gcomp_right)
# $gcomp_est
# [1] 0.6493520 0.4004721 0.2488800
#
# $gcomp_ci
# E[Q1hat]  E[Q0hat]       ATE
# [1,] 0.5878254 0.3657378 0.1833215
# [2,] 0.7034692 0.4372652 0.3097588



res_gcom = cbind(est_gcomp_right[2:1],ci_gcomp_right[,2], ci_gcomp_right[,1])
colnames(res_gcom) = c("est", "cil", "ciu"); rownames(res_gcom) = c("0", "1")
#       est       cil       ciu
# 0 0.4004721 0.3655522 0.5962230
# 1 0.6493520 0.4317182 0.7041806





############ IPTW ##########

iptw <- function(W, A, Y, formula = "A ~ .", family = binomial()){
  ps_fit <- glm(as.formula(formula), data = data.frame(W, A = A), family = family)
  g1hat <- predict(ps_fit, type = "response")
  est <- sum(A * Y /(g1hat))/sum(A/g1hat) - sum((1-A)*Y/(1- g1hat))/sum((1-A)/(1-g1hat))
  return(c(sum(A * Y /(g1hat))/sum(A/g1hat), sum((1-A)*Y/(1- g1hat))/sum((1-A)/(1-g1hat)), est))
}

iptw_one_boot <- function(W, A, Y, formula = "A ~ .", family = binomial()){
  resamp_idx <- sample(1:length(A), replace = TRUE)
  return(iptw(W[resamp_idx, , drop = FALSE], A = A[resamp_idx], Y = Y[resamp_idx],
              formula = formula, family = family))
}

iptw_boot_ci <- function(W, A, Y, formula = "Y ~ .", nboot = 5e2, family = binomail()){
  ate_star <- replicate(nboot, iptw_one_boot(W, A, Y, formula, family = family))
  ci = sapply(1:3, function(x) as.numeric(quantile(ate_star[x,], p = c(0.025, 0.975))))
  colnames(ci) = c("E[Q1hat]", "E[Q0hat]", "ATE")
  return(ci) # 1st col for Q1hat, 2nd col for Q0hat, 3rd col for ATE
  # return(as.numeric(quantile(ate_star, p = c(0.025, 0.975))))
}


est_iptw_right <- iptw(W, A, Y, formula = correct_gform, family = sim_fam)
ci_iptw_right <- iptw_boot_ci(W, A, Y, formula = correct_gform, family = sim_fam)
iptw_right <- list(iptw_est = est_iptw_right, iptw_ci = ci_iptw_right)
# iptw_est
# [1] 0.6870389 0.3959046 0.2911343
#
# $iptw_ci
#       E[Q1hat]  E[Q0hat]       ATE
# [1,] 0.6002881 0.3642765 0.1918264
# [2,] 0.7745947 0.4328871 0.3790200



res_iptw = cbind(est_iptw_right[2:1],ci_iptw_right[,2], ci_iptw_right[,1])
colnames(res_iptw) = c("est", "cil", "ciu"); rownames(res_iptw) = c("0", "1")
#       est       cil       ciu
# 0 0.3959046 0.3643497 0.6064910
# 1 0.6870389 0.4325825 0.7663278



