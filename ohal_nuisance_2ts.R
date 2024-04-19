
# use HAL to estimate Qn
# use OHAL based on Qn to estimate Gn
# for two-time step
# extension on Ju's function ohal_nuisance() for one-time step




#--------------------------- main function ------------------------
# main function to estimate Qn and Gn
ohal_nuisance_2ts <- function(L0, A0, L1, A1, L2, abar, V = 10, outcome_family = "binomial", penalty_threshold = 1e-10,
                              penalty_gamma = 1, newdata = NULL, lambda_seq = exp(seq(0, -13, length=2000)), tolg = tolg){

  n <- length(L2)
  # make a vector of cv-folds

  fold_idx <- chunk2(seq_len(n), V) # output V chunks: i.e. $`9` 801-900
  folds <- rep(seq_len(V), unlist(lapply(fold_idx, length)))  # 1...1 (100 1s), 2...2 (100 2s),...
  # fold_idx <- unlist(fold_idx, use.names = FALSE)  # 1,2,..., 1000


  ################### outcome
  # cross-validation for OR
  cv_out <- sapply(seq_len(V), one_hal_or_m2ts,  # use one_hal_or(), output $QAW, $riskQ (MSE), $Q1W, $Q0W, $or_coef (coefs from hal of Q)
                   folds = folds, L0 = L0, A0 = A0, L1 = L1, A1 = A1, L2 = L2,
                   lambda_seq = lambda_seq, outcome_family = outcome_family, simplify = FALSE)

  risk_all = esum_list(lapply(cv_out, "[[", "risk") )
  lambda_or_idx <- which(risk_all == min(risk_all), arr.ind = TRUE)# take smallest lambda that has smallest risk


  # re-fit on full data using the optimal lambda
  full_hal_or = one_hal_or_m2ts(validFold = NULL, folds = NULL, L0 = L0, A0 = A0, L1 = L1, A1 = A1, L2 = L2,
                               lambda_seq = lambda_seq, outcome_family = outcome_family)

  Q2n_select <- full_hal_or$Qn$Q2n[, lambda_or_idx[1]]
  Q1n_select <- full_hal_or$Qn$Q1n[[lambda_or_idx[1]]][, lambda_or_idx[2]]


  # get cv selected OR coefficients for each fold
  Q2mod_coef_list <- lapply(cv_out, function(x){ x$coefs$Q2mod[, lambda_or_idx[1]] })  # length(cv_out)=length(or_coef_list)=5, different length for each or_coef_list[[a]] for a = 1,...,5
  Q1mod_coef_list <- lapply(cv_out, function(x){ x$coefs$Q1mod[[lambda_or_idx[1]]][, lambda_or_idx[2]] })
  # sapply(1:5, function(x) summary(Q1mod_coef_list[[x]]))
  # sapply(1:5, function(x) length(Q2mod_coef_list[[x]])) #  2387 2389 2386 2386 2387
  # sapply(1:5, function(x) length(Q1mod_coef_list[[x]]))  # 801 801 801 801 801


  ################### PS

  cv_out_ps <- mapply(seq_len(V), Q2mod_coef_list = Q2mod_coef_list, Q1mod_coef_list = Q1mod_coef_list, FUN = one_hal_ps_2ts,  # use one_hal_ps(), output
                      MoreArgs = list(folds = folds, L0 = L0, A0 = A0, L1 = L1, A1 = A1, L2 = L2, penalty_gamma = penalty_gamma,  # should define gamma value
                                      lambda_seq = lambda_seq, penalty_threshold = penalty_threshold, tolg = tolg),
                      SIMPLIFY = FALSE)

  ohal_risks1 <- colMeans(Reduce("rbind",lapply(cv_out_ps, "[[", "G1n_OHAL_risk")))
  ohal_risks0 <- colMeans(Reduce("rbind",lapply(cv_out_ps, "[[", "G0n_OHAL_risk")))
  # take smallest lambda that has smallest risk
  ohal_lambda_ps1_idx <- which.min(ohal_risks1)[1]
  ohal_lambda_ps0_idx <- which.min(ohal_risks0)[1]

  # re-fit on full data
  full_hal_ps <- one_hal_ps_2ts(validFold = NULL, folds = NULL, L0 = L0, A0 = A0, L1 = L1, A1 = A1, L2 = L2,
                                penalty_gamma = penalty_gamma, lambda_seq = lambda_seq, newdata = newdata,
                                Q2mod_coef_list = full_hal_or$coefs$Q2mod[, lambda_or_idx[1]] ,  # !!! using the lambda correspongding to min(risk of Q) in full Q model
                                Q1mod_coef_list = full_hal_or$coefs$Q1mod[[lambda_or_idx[1]]][, lambda_or_idx[2]],
                                penalty_threshold = penalty_threshold, tolg = tolg)

  # G1W at cv selected lambda
  G1n_OHAL_select <- full_hal_ps$G1n_OHAL[, ohal_lambda_ps1_idx]
  G0n_OHAL_select <- full_hal_ps$G0n_OHAL[, ohal_lambda_ps0_idx]


  ###################  format outcome
  out <- data.frame(Q2n = Q2n_select, Q1n = Q1n_select, G1n = G1n_OHAL_select, G0n = G0n_OHAL_select, folds = folds)

  return(out)
}


#--------------------- helper functions ---------------------
chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))



esum_list <- function(list) {
  ll = length(list); suml = 0
  for (i in 1:ll){  suml = suml + as.matrix(list[[i]]) }
  return(suml)
}


# A helper function for fitting a single fold of HAL on outcome. Can be used both for fitting cross-validated HAL and fitting HAL to full data.

one_hal_or_m2ts <- function(validFold, folds, L0, A0, L1, A1, L2, outcome_family = "binomial",
                            newdata = NULL, lambda_seq = exp(seq(0, -13, length=2000))){
  n <- length(L2)
  n_valid <- sum(validFold == folds)
  n_train <- n - n_valid

  # if called during cross-validation
  if(is.null(folds)){
    train_L0 <- valid_L0 <- L0; train_L1 <- valid_L1 <- L1; train_A0 <- valid_A0 <- A0
    train_A1 <- valid_A1 <- A1; train_L2 <- valid_L2 <- L2
  }else{
    # training data
    train_L0 <- L0[folds != validFold]; train_L1 <- L1[folds != validFold]
    train_A0 <- A0[folds != validFold]; train_A1 <- A1[folds != validFold]; train_L2 <- L2[folds != validFold]

    # validation data
    valid_L0 <- L0[folds == validFold]; valid_L1 <- L1[folds == validFold]
    valid_A0 <- A0[folds == validFold]; valid_A1 <- A1[folds == validFold]; valid_L2 <- L2[folds == validFold]
  }


  # for outcome regression
  #--------
  # Q2n
  #--------
  # x2_fit = cbind(train_L0,train_L1)[train_A0==abar[1] & train_A1==abar[2],]
  # y2_fit = train_L2[train_A0==abar[1] & train_A1==abar[2]]
  # new_x2 = cbind(L0,L1)
  x2_fit = cbind(train_L0,train_L1)
  y2_fit = as.matrix(train_L2)
  new_x2 = cbind(L0,L1)
  basis_list <- enumerate_basis(x2_fit, NULL) # length(basis_list): 9608
  x2_basis <- hal9001:::make_design_matrix(x2_fit, basis_list)  # dim(x2_basis): 226 678
  copy_map <- hal9001:::make_copy_map(x2_basis) # length(copy_map) 652, output a list of numeric vectors indicating indices of basis functions that are identical in the training set
  unique_columns <- as.numeric(names(copy_map)) # length(unique_columns) 652
  # subset to non-duplicated columns
  x2_basis_fit <- x2_basis[, unique_columns]
  dim(x2_basis_fit)
  x2_basis_fit = cbind(train_A0, train_A1, x2_basis_fit) # dim(x2_basis_fit)

  ##### fit outcome regression
  hal_lasso2 <- glmnet::glmnet(x = x2_basis_fit, y = y2_fit, family = outcome_family, lambda = lambda_seq, standardize = FALSE)

  new_x_basis2 <- hal9001:::make_design_matrix(new_x2, basis_list)
  new_x_basis2 <- as.matrix(new_x_basis2[, unique_columns])  # dim(new_x_basis2): n_valid * 652

  or_pred_matrix2 <- cbind(rep(1, ifelse(is.null(validFold), n, n_valid)), cbind(A0=abar[1],A1=abar[2], new_x_basis2)) %*% rbind(hal_lasso2$a0, hal_lasso2$beta) # dim(or_pred_matrix2)  n * length(lambda_seq)

  if(outcome_family == "binomial"){ or_pred_matrix2 <- apply(or_pred_matrix2, 2, plogis) } # dim(or_pred_matrix2)  n * length(lambda_seq)

  #--------
  # Q1n
  #--------
  # x1_fit = as.matrix(train_L0[train_A0==abar[1]]);
  # new_x1 = as.matrix(valid_L0)
  x1_fit = as.matrix(train_L0);
  new_x1 = as.matrix(valid_L0)
  basis_list <- enumerate_basis(x1_fit, NULL)
  x1_basis <- hal9001:::make_design_matrix(x1_fit, basis_list)  # dim: n_train * 445
  copy_map <- hal9001:::make_copy_map(x1_basis) # length(copy_map) 445, output a list of numeric vectors indicating indices of basis functions that are identical in the training set
  unique_columns <- as.numeric(names(copy_map)) # length(unique_columns) 445
  # subset to non-duplicated columns
  x1_basis_fit <- x1_basis[, unique_columns] # dim(x1_basis_fit)
  x1_basis_fit = cbind(train_A0, x1_basis_fit) # dim(x2_basis_fit)

  risks_allamb = NULL
  or_pred_list1 = Q1mod_coef_list = list()
  new_x_basis1 <- hal9001:::make_design_matrix(new_x1, basis_list)
  new_x_basis1 <- as.matrix(new_x_basis1[, unique_columns])  # dim(new_x_basis1): n_valid * 445

  for (j in 1:length(lambda_seq)){

    if(!is.null(folds)){ y1_fit <-  or_pred_matrix2[,j][folds != validFold] }else{ y1_fit <- or_pred_matrix2[,j] }

    hal_lasso1 <- glmnet::glmnet(x = x1_basis_fit, y = y1_fit, lambda = lambda_seq, standardize = FALSE)

    or_pred_matrix1 <- cbind(rep(1, ifelse(is.null(validFold), n, n_valid)), cbind(A0=abar[1],new_x_basis1)) %*% rbind(hal_lasso1$a0, hal_lasso1$beta) # dim(or_pred_matrix2)  n * length(lambda_seq)
    # if(outcome_family == "binomial"){ or_pred_matrix1 <- apply(or_pred_matrix1, 2, plogis) } # dim(or_pred_matrix2)

    # compute MSE/loglik
    if(is.null(newdata)){ or_risk <- apply(or_pred_matrix1, 2, function(x){mean((L2[folds == validFold] - x)^2)}) } else { or_risk <- NULL }
    risks_allamb = rbind(risks_allamb, or_risk) # dim(risks_allamb)  length(lambda_seq) X length(lambda_seq) where row refers to the lambda in Q2n, column referes to the lambda in Q1n

    or_pred_list1[[j]] = or_pred_matrix1
    Q1mod_coef_list[[j]] = hal_lasso1$beta

  }


  ###### format output
  out <- list()
  out$Qn$Q2n <- or_pred_matrix2; out$Qn$Q1n <- or_pred_list1
  out$risk = risks_allamb
  out$coefs$Q2mod = hal_lasso2$beta; out$coefs$Q1mod = Q1mod_coef_list


  return(out)

}



# A helper function for fitting a single fold of OHAL on treatment. Can be used both for fitting cross-validated HAL and fitting HAL to full data.

one_hal_ps_2ts <- function(validFold, folds, L0, A0, L1, A1, L2, Q2mod_coef_list, Q1mod_coef_list, # passed in from HAL OR
                           penalty_gamma = 1, penalty_threshold = 1e-5, newdata = NULL, lambda_seq = exp(seq(0, -13, length=2000)), tolg = 1e-5){
  n <- length(L2)
  n_valid <- sum(validFold == folds) # 200 if V=5, n=1000
  n_train <- n - n_valid

  # if called during cross-validation
  if(is.null(folds)){
    train_L0 <- valid_L0 <- L0; train_L1 <- valid_L1 <- L1; train_A0 <- valid_A0 <- A0
    train_A1 <- valid_A1 <- A1; train_L2 <- valid_L2 <- L2
  }else{
    # training data
    train_L0 <- L0[folds != validFold]; train_L1 <- L1[folds != validFold]
    train_A0 <- A0[folds != validFold]; train_A1 <- A1[folds != validFold]; train_L2 <- L2[folds != validFold]

    # validation data
    valid_L0 <- L0[folds == validFold]; valid_L1 <- L1[folds == validFold]
    valid_A0 <- A0[folds == validFold]; valid_A1 <- A1[folds == validFold]; valid_L2 <- L2[folds == validFold]
  }


  # for PS
  #--------
  # g0n
  #--------
  x0_fit = as.matrix(train_L0)
  a0_fit = train_A0
  new_x0 = as.matrix(L0)

  basis_list <- enumerate_basis(x0_fit, NULL)
  x0_basis <- hal9001:::make_design_matrix(x0_fit, basis_list)  # dim: n_train * 800
  copy_map <- hal9001:::make_copy_map(x0_basis) # length(copy_map) 800, output a list of numeric vectors indicating indices of basis functions that are identical in the training set
  unique_columns <- as.numeric(names(copy_map)) # length(unique_columns) 800
  # subset to non-duplicated columns
  x0_basis_fit <- x0_basis[, unique_columns] # dim(x0_basis_fit)

  # Q coefficients
  # !!! the number of coefs from Q should equal to the dimension of covariates included into the model
  or_coef_tmp0 <- abs(Q1mod_coef_list[-1])
  penalty.factor <- or_coef_tmp0^(-penalty_gamma) # length(penalty.factor)= 451, could be different across different folds
  penalty.factor[or_coef_tmp0 < penalty_threshold] <- Inf # sum(penalty.factor== Inf)=449

  if(!all(or_coef_tmp0 == 0)){
    # fit ohal
    ohal_lasso0 <- glmnet::glmnet(x = x0_basis_fit, y = a0_fit, lambda = lambda_seq, family = 'binomial',
                                  standardize = FALSE, penalty.factor = penalty.factor)
  }else{
    ohal_lasso0 <- list(beta = matrix(0, nrow = dim(x0_basis_fit)[2], ncol = length(lambda_seq)), a0 = qlogis(mean(a0_fit)))
  }

  new_x_basis0 <- hal9001:::make_design_matrix(new_x0, basis_list)
  new_x_basis0 <- as.matrix(new_x_basis0[, unique_columns])

  ps_pred_matrix0 <- cbind(rep(1, ifelse(is.null(validFold), n, n_valid)), new_x_basis0) %*% rbind(ohal_lasso0$a0, ohal_lasso0$beta) # dim(ps_pred_matrix0)  n * length(lambda_seq)
  ps_pred_matrix0 <- apply(ps_pred_matrix0, 2, plogis)  # dim(or_pred_matrix2)  n * length(lambda_seq)
  ps_pred_matrix0[ps_pred_matrix0 < tolg] <- tolg
  ps_pred_matrix0[ps_pred_matrix0 > 1 - tolg] <- 1 - tolg

  #--------
  # g1n
  #--------
  x1_fit = cbind(train_L0, train_L1)
  a1_fit = train_A1
  new_x1 = cbind(L0, L1)

  basis_list <- enumerate_basis(x1_fit, NULL)
  x1_basis <- hal9001:::make_design_matrix(x1_fit, basis_list)  # dim: n_train * 800
  copy_map <- hal9001:::make_copy_map(x1_basis) # length(copy_map) 800, output a list of numeric vectors indicating indices of basis functions that are identical in the training set
  unique_columns <- as.numeric(names(copy_map)) # length(unique_columns) 800
  # subset to non-duplicated columns
  x1_basis_fit <- x1_basis[, unique_columns] # dim(x1_basis_fit)
  x1_basis_fit = cbind(train_A0, x1_basis_fit)

  # Q coefficients
  # !!! the number of coefs from Q should equal to the dimension of covariates included into the model
  or_coef_tmp2 <- abs(Q2mod_coef_list[-2])
  penalty.factor <- or_coef_tmp2^(-penalty_gamma) # length(penalty.factor)= 2388, could be different across different folds
  penalty.factor[or_coef_tmp2 < penalty_threshold] <- Inf # sum(penalty.factor== Inf)=1918


  if(!all(or_coef_tmp2 == 0)){
    # fit ohal
    ohal_lasso1 <- glmnet::glmnet(x = x1_basis_fit, y = a1_fit, family = 'binomial', lambda = lambda_seq,
                                  standardize = FALSE, penalty.factor = penalty.factor)
  }else{
    ohal_lasso1 <- list(beta = matrix(0, nrow = dim(x_basis_fit)[2], ncol = length(lambda_seq)), # dim(beta) 8220 *length(lambda_seq)
                        a0 = qlogis(mean(a1_fit)))
  }


  # predictions on validation sample
  new_x_basis1 <- hal9001:::make_design_matrix(new_x1, basis_list)
  new_x_basis1 <- as.matrix(new_x_basis1[, unique_columns])

  ps_pred_matrix1 <- cbind(rep(1, ifelse(is.null(validFold), n, n_valid)), cbind(A0=abar[1],new_x_basis1)) %*% rbind(ohal_lasso1$a0, ohal_lasso1$beta) # dim(ps_pred_matrix0)  n * length(lambda_seq)
  ps_pred_matrix1 <- apply(ps_pred_matrix1, 2, plogis)  # dim(or_pred_matrix2)  n * length(lambda_seq)

  ps_pred_matrix1[ps_pred_matrix1 < tolg] <- tolg
  ps_pred_matrix1[ps_pred_matrix1 > 1 - tolg] <- 1 - tolg


  # compute MSE/loglik
  if(is.null(newdata)){
    ohal_ps_risk0 <- apply(ps_pred_matrix0, 2, function(x){mean(ifelse(A0[folds == validFold] == 1, -log(x), -log(1 - x))) })
    ohal_ps_risk1 <- apply(ps_pred_matrix1, 2, function(x){mean(ifelse(A1[folds == validFold] == 1, -log(x), -log(1 - x))) })
  }else{ ohal_ps_risk1 <- ohal_ps_risk0 <- NULL }


  ###### format output
  out <- list()
  out$G1n_OHAL <- ps_pred_matrix1
  out$G0n_OHAL <- ps_pred_matrix0

  if(!is.null(validFold)){
    out$G1n_OHAL_risk <- ohal_ps_risk1
    out$G0n_OHAL_risk <- ohal_ps_risk0
  }
  return(out)
}




########################### testing using simple setting

# O=(L0, A0, L1, A1, L2)  use the same data structure with drinf.tmle()

library(hal9001)
library(SuperLearner)

n=1000

L0 <- runif(n,-1,1)
A0 <- rbinom(n, 1, plogis(L0^2))
L1 <- L0^2*A0 + runif(n)
A1 <- rbinom(n, 1, plogis(L0*L1))
L2 <- rbinom(n, 1, plogis(L0^2*A0*A1 + L1))

abar=c(1,1)


fit_nuis <- ohal_nuisance_2ts(L0, A0, L1, A1, L2, abar = abar, V = 2,
                              outcome_family = "binomial", penalty_threshold = 1e-5, penalty_gamma = 1,
                              newdata = NULL, lambda_seq = exp(seq(-3.5, -13, length=50)), tolg = 1e-5)
head(fit_nuis)
summary(fit_nuis)












#' #' A helper function for fitting a single fold of HAL. Can be used both for fitting cross-validated HAL and fitting HAL to full data.
#' one_hal_or_2ts <- function(validFold, folds, L0, A0, L1, A1, L2, outcome_family = "binomial",
#'                            newdata = NULL, lambda_seq = exp(seq(-1,-13,length=1000))){
#'   n <- length(L2)
#'   # n_valid <- sum(validFold == folds)
#'   # n_train <- n - n_valid
#'
#'     # if called during cross-validation
#'     if(is.null(folds)){
#'       train_L0 <- valid_L0 <- L0; train_L1 <- valid_L1 <- L1; train_A0 <- valid_A0 <- A0
#'       train_A1 <- valid_A1 <- A1; train_L2 <- valid_L2 <- L2
#'     }else{
#'       # training data
#'       train_L0 <- L0[folds != validFold]; train_L1 <- L1[folds != validFold]
#'       train_A0 <- A0[folds != validFold]; train_A1 <- A1[folds != validFold]; train_L2 <- L2[folds != validFold]
#'
#'       # validation data
#'       valid_L0 <- L0[folds == validFold]; valid_L1 <- L1[folds == validFold]
#'       valid_A0 <- A0[folds == validFold]; valid_A1 <- A1[folds == validFold]; valid_L2 <- L2[folds == validFold]
#'     }
#'
#'
#'     # for outcome regression
#'     #------------------------
#'     Q2mod <- do.call(SL.Q, args=c(list(
#'       Y=eval(parse(text=paste0( ifelse(stratify, "train_L2[train_A0==abar[1] & train_A1==abar[2]]", "train_L2"  ) ))),
#'       X=eval(parse(text=paste0( ifelse(stratify, "cbind(train_L0,train_L1)[train_A0==abar[1] & train_A1==abar[2],]",
#'                "cbind(train_L0,train_L1,train_A0,train_A1)" ) ))),
#'       newX = eval(parse(text=paste0( ifelse(stratify,  "cbind(L0,L1)", "cbind(L0,L1,A0=abar[1],A1=abar[2])" ) ))),
#'       SL.library=SL.Q,
#'       obsWeights = rep(1, ifelse(stratify, length(train_L2[train_A0==abar[1] & train_A1==abar[2]]),length(train_L2))),
#'       lambda = lambda_seq, verbose=verbose), SL.Q.options
#'     ))
#'
#'
#'     Q2n_coefs = Q2mod$fit$object$coefs   # length(Q2n_coefs)  228
#'     Q2n_all <- as.numeric(Q2mod$pred)  # length(Q2n_all)  n
#'     # replace extrapolated predictions with smallest/largest value
#'     Q2n_all[Q2n_all < min(L2)] <- min(L2)
#'     Q2n_all[Q2n_all > max(L2)] <- max(L2)
#'
#'     if(!is.null(folds)){ train_Q2n <- Q2n_all[folds != validFold] }else{ train_Q2n <- Q2n_all }
#'
#'     Q1mod <- do.call(SL.Q, args = c(list(
#'       Y=eval(parse(text=paste0( ifelse(stratify, "train_Q2n[train_A0==abar[1]]",  "train_Q2n" ) ))),
#'       X=eval(parse(text=paste0( ifelse(stratify, "train_L0[train_A0==abar[1]]", "cbind(train_L0,train_A0)" ) ))),
#'       newX = eval(parse(text=paste0( ifelse(stratify, "valid_L0", "cbind(valid_L0,A0=abar[1])" ) ))),
#'       SL.library = SL.Q,
#'       obsWeights = rep(1, ifelse(stratify, length(train_Q2n[train_A0==abar[1]]),length(train_Q2n))),
#'       lambda = lambda_seq, verbose = verbose), SL.Q.options
#'     ))
#'
#'     Q1n_coefs = Q1mod$fit$object$coefs  # length(Q1n_coefs)  26
#'     Q1n <- as.numeric(Q1mod$pred) # length(Q1n)  n/V
#'     # replace extrapolated predictions with  smallest/largest value
#'     Q1n[Q1n < min(Q2n)] <- min(Q2n)
#'     Q1n[Q1n > max(Q2n)] <- max(Q2n)
#'
#'
#'     # compute MSE/loglik
#'     if(is.null(newdata)){
#'       if(outcome_family == "gaussian"){ or_risk <- mean((L2[folds == validFold] - Q1n)^2)
#'       } else { or_risk <-  mean(ifelse(L2[folds == validFold] == 1, -log(Q1n), -log(1 - Q1n)))  }
#'     }else{ or_risk <- NULL }
#'
#'     ###### format output
#'     out <- list()
#'     out$Qn$Q2n <- Q2n; out$Qn$Q1n <- Q1n
#'     out$risk = or_risk
#'     out$Lambda_opt <- data.frame(Q2mod_lambda = Q2mod$fit$object$lambda_star, Q1mod_lambda = Q1mod$fit$object$lambda_star)
#'     out$coefs$Q2mod = Q2n_coefs; out$coefs$Q1mod = Q1n_coefs
#'     out$Q2mod = Q2mod; out$Q1mod = Q1mod
#'
#'   return(out)
#'
#' }




# for ps estimation
#' A helper function for fitting a single fold of HAL. Can be used both for fitting cross-validated HAL and fitting HAL to full data.



