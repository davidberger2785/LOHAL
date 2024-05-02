# target Q at two time steps, extended on drinf.tmle package

#------------------------------------------------------------------
#' targetQ2: Function that targets Q2n in two steps. First solving original equations, then solving sum of two new equations.
#'
#' @param A0 A \code{vector} treatment delivered at baseline.
#' @param A1 A \code{vector} treatment deliver after \code{L1} is measured.
#' @param L2 A \code{vector} outcome of interest.
#' @param Qn A \code{list} of current estimates of Q2n and Q1n
#' @param gn A \code{list} of current estimates of g1n and g0n
#' @param gnr A \code{list} of current estimates of reduced dim. regressions
#' @param abar A \code{vector} of length 2 indicating the treatment assignment that is of interest.
#' @param tolg A \code{numeric} indicating the truncation level for conditional treatment probabilities.
#' @param tolQ A \code{numeric}
#' @param return.models  A \code{boolean} indicating whether the fluctuation model should be returned with the output.
#' @param tol.coef A \code{numeric} indicating the coefficient above which the minimization along the
#' submodel using \code{glm} is deemed to be unreasonable. In these cases \code{optim} is used
#' instead to perform the fluctuation along the same submodel.
#' @importFrom SuperLearner trimLogit
#'
#' @return A list with named entries corresponding to the estimators of the fluctuated nuisance parameters evaluated at the observed data values. If
#' \code{return.models = TRUE} output also includes the fitted fluctuation model.

targetQ2 <- function( A0, A1, L2, Qn, gn, gnr, abar, tolg, tolQ, return.models, tol.coef = 1e2, ... ){
  #-------------------------------------------
  # making outcomes for logistic fluctuation
  #-------------------------------------------
  # length of output
  n <- length(gn$g0n)

  # scale L2, Q2n, Q1n to be in (0,1)
  L2.min <- min(L2); L2.max <- max(L2)
  # scale L2
  L2s <- (L2 - L2.min)/(L2.max - L2.min)
  # scale Q2n,Q1n
  Q2ns <- (Qn$Q2n - L2.min)/(L2.max - L2.min)

  #-------------------------------------------
  # making offsets for logistic fluctuation
  #-------------------------------------------
  flucOff <- c( SuperLearner::trimLogit(Q2ns, trim = tolQ) )

  #-------------------------------------------
  # making covariates for fluctuation
  #-------------------------------------------
  # the original "clever covariates"
  combinedPropensity <- gn$g0n * gn$g1n
  combinedPropensity[combinedPropensity < tolg] <- tolg

  flucCov1 <- c( as.numeric(A0==abar[1] & A1==abar[2])/(combinedPropensity) ) # the usual guy

  # the new "clever covariates" for Q
  flucCov2 <- c(as.numeric(A0==abar[1] & A1==abar[2]) * ((gnr$hbarnr +gnr$h1nr)/gnr$g1nr) )

  #-------------------------------------------
  # making covariates for prediction
  #-------------------------------------------
  # getting the values of the clever covariates evaluated at \bar{A} = abar
  predCov1 <- c( 1/combinedPropensity  ) # all c(A0,A1) = abar

  predCov2 <- c( (gnr$hbarnr + gnr$h1nr)/gnr$g1nr ) # all c(A0,A1) = abar

  #-------------------------------------------
  # fitting fluctuation submodel

  #-------- second fluctuation submodel to solve new equations
  flucmod2 <- suppressWarnings(glm( formula = "out ~ -1 + offset(fo) + fc1",
                                    data = data.frame(out = L2s, fo = flucOff, fc1 = flucCov2), family = binomial(), start = 0 ))

  if(abs(flucmod2$coefficients) < tol.coef){
    # get predictions
    tmp <- predict(  flucmod2, newdata = data.frame(out = 0, fo = flucOff, fc1 = predCov2), type = "response" )
  }else{
    # use optim to try the minimization along intercept only submodel if glm looks wonky
    flucmod2 <- optim( par = 0, fn = wnegloglik, gr = gradient.wnegloglik, method = "L-BFGS-B", lower = -tol.coef,
                       upper = tol.coef, control = list(maxit = 10000), Y = L2s, offset = flucOff, weight = flucCov2 )
    epsilon <- flucmod2$par
    tmp <- plogis(flucOff +  epsilon)
  }

  #---------- first fluctuation submodel to solve original equations
  # new offset
  flucOff <- c( SuperLearner::trimLogit(tmp, trim = tolQ) )

  flucmod1 <- suppressWarnings(glm( formula = "out ~ -1 + offset(fo) + fc1",
                                    data = data.frame(out = L2s, fo = flucOff,  fc1 = flucCov1), family = binomial(), start = 0 ))
  # see if the fluctuation coefficient is reasonable
  if(abs(flucmod1$coefficients) < tol.coef){
    # get predictions
    Q2nstar <- predict( flucmod1, newdata = data.frame(out = 0, fo = flucOff, fc1 = predCov1), type = "response" )*(L2.max - L2.min) + L2.min
  }else{
    # use optim to try the minimization along submodel if glm looks wonky
    flucmod1 <- optim( par = 0, fn = offnegloglik, gr = gradient.offnegloglik, method = "L-BFGS-B",
                       lower = -tol.coef, upper = tol.coef, control = list(maxit = 10000), Y = L2s, offset = flucOff, weight = flucCov1 )
    epsilon <- flucmod1$par
    Q2nstar <- stats::plogis(flucOff + predCov1 * epsilon)*(L2.max - L2.min) + L2.min
  }

  #--------------
  # output
  #-------------
  out <- list( Q2nstar = Q2nstar, flucmod = NULL )
  if(return.models){ out$flucmod = list(flucmod1, flucmod2) }

  return(out)
}





#' targetQ1: Function that targets Q1n in two steps. First solving original equations, then solving sum of two new equations.
#'
targetQ1 <- function( A0, A1, L2, Qn, gn, gnr, abar, tolg, tolQ, return.models, tol.coef = 1e2, ... ){
  #-------------------------------------------
  # making outcomes for logistic fluctuation
  #-------------------------------------------
  # length of output
  n <- length(gn$g0n)

  # scale L2, Q2n, Q1n to be in (0,1)
  L2.min <- min(L2); L2.max <- max(L2)
  # scale Q2n,Q1n
  Q2ns <- (Qn$Q2n - L2.min)/(L2.max - L2.min)
  Q1ns <- (Qn$Q1n - L2.min)/(L2.max - L2.min)
  #-------------------------------------------
  # making offsets for logistic fluctuation
  #-------------------------------------------
  flucOff <- c( SuperLearner::trimLogit(Q1ns, trim = tolQ) )

  #-------------------------------------------
  # making covariates for fluctuation
  #-------------------------------------------
  # the original "clever covariates"
  flucCov1 <- c( as.numeric(A0 == abar[1])/gn$g0n  ) # the usual guy

  ### --- the new "clever covariates" for Q
  flucCov2 <- c(as.numeric(A0 == abar[1]) * (gnr$h0nr/gnr$g0nr) )

  #-------------------------------------------
  # making covariates for prediction
  #-------------------------------------------
  # getting the values of the clever covariates evaluated at \bar{A} = abar
  predCov1 <- c( 1/gn$g0n ) # all A0 == abar[1]

  ### --- the new prediction for Q
  predCov2 <- c( gnr$h0nr/gnr$g0nr) # all A0 == abar[1]

  #-------------------------------------------
  # fitting fluctuation submodel

  #---------- second fluctuation submodel to solve new equations
  flucmod2 <- suppressWarnings( glm( formula = "out ~ -1 + offset(fo) + fc1",
                                     data = data.frame(out = Q2ns, fo = flucOff, fc1 = flucCov2), family = binomial(), start = 0 ))

  if(abs(flucmod2$coefficients) < tol.coef){
    # get predictions
    tmp <- predict( flucmod2, newdata = data.frame(out = 0, fo = flucOff, fc1 = predCov2), type = "response" )
  }else{
    # use optim to try the minimization along intercept only submodel if glm looks wonky
    flucmod2 <- optim(
      par = 0, fn = offnegloglik, gr = gradient.offnegloglik, method = "L-BFGS-B", lower = -tol.coef,
      upper = tol.coef,control = list(maxit = 10000), Y = Q2ns, offset = flucOff, weight = flucCov2 )
    epsilon <- flucmod2$par
    tmp <- plogis(flucOff + predCov2 * epsilon)*(L2.max - L2.min) + L2.min
  }


  #------------ first fluctuation submodel to solve original equations
  # new offset
  flucOff <- c( SuperLearner::trimLogit(tmp, trim = tolQ) )

  flucmod1 <- suppressWarnings(glm( formula = "out ~ -1 + offset(fo) + fc1",
                                    data = data.frame(out = Q2ns, fo = flucOff, fc1 = flucCov1), family = binomial(), start = 0 ))

  if(abs(flucmod1$coefficients) < tol.coef){
    # get predictions
    Q1nstar <- predict( flucmod1, newdata = data.frame(out = 0, fo = flucOff, fc1 = predCov1), type = "response" )*(L2.max - L2.min) + L2.min
  }else{
    # use optim to try the minimization along submodel if glm looks wonky
    flucmod1 <- optim( par = 0, fn = offnegloglik, gr = gradient.offnegloglik, method = "L-BFGS-B", lower = -tol.coef,
                       upper = tol.coef, control = list(maxit = 10000), Y = Q2ns, offset = flucOff, weight = flucCov1)
    epsilon <- flucmod1$par
    Q1nstar <- plogis(flucOff + predCov1 * epsilon)*(L2.max - L2.min) + L2.min
  }

  #--------------
  # output
  #-------------
  out <- list( Q1nstar = Q1nstar, flucmod = NULL )
  if(return.models){ out$flucmod = list(flucmod1, flucmod2) }

  return(out)
}



#' wnegloglik: A function that computes the weighted negative log-likelihood loss of theintercept only fluctuation model.
#'
#' @param epsilon The scalar parameter of the fluctuation submodel.
#' @param weight The \code{vector} of weights, i.e., the clever covariates.
#' @param Y The \code{vector} of regression outcomes.
#' @param offset The \code{vector} of offsets.
#' @return A \code{numeric} value of the negative log-likelihood loss

wnegloglik <- function(epsilon, weight, Y, offset){
  X <- as.matrix(cbind(offset, rep(1, length(offset))))
  mu <- plogis(X%*%as.matrix(c(1,epsilon)))
  mu[mu == 0] <- .Machine$double.neg.eps
  mu[mu == 1] <- 1 - .Machine$double.neg.eps
  wloglik <- sum(weight * (Y * log(mu) + (1-Y)*log(1-mu)))
  return(-wloglik)
}


offnegloglik <- function(epsilon, weight, Y, offset){
  X <- as.matrix(cbind(offset, weight))
  mu <- plogis(X%*%as.matrix(c(1,epsilon)))
  mu[mu == 0] <- .Machine$double.neg.eps
  mu[mu == 1] <- 1 - .Machine$double.neg.eps
  ologlik <- sum((Y * log(mu) + (1-Y)*log(1-mu)))
  return(-ologlik)
}

gradient.offnegloglik <- function(epsilon, weight, Y, offset){
  X <- as.matrix(cbind(offset, weight))
  mu <- plogis(X%*%matrix(c(1,epsilon)))
  mu[mu == 0] <- .Machine$double.neg.eps
  mu[mu == 1] <- 1 - .Machine$double.neg.eps
  grad <- crossprod(weight, Y - mu)
  return(-grad)
}


#' gradient.wnegloglik
#' @param epsilon The scalar parameter of the fluctuation submodel.
#' @param weight The \code{vector} of weights, i.e., the clever covariates.
#' @param Y The \code{vector} of regression outcomes.
#' @param offset The \code{vector} of offsets.
#'
#' @return A \code{vector} of the gradient of the loss evaluated at the data inputs.
gradient.wnegloglik <- function(epsilon, weight, Y, offset){
  X <- cbind(offset, rep(1,length(offset)))
  mu <- plogis(X%*%matrix(c(1,epsilon)))
  mu[mu == 0] <- .Machine$double.neg.eps
  mu[mu == 1] <- 1 - .Machine$double.neg.eps
  grad <- crossprod(weight, Y - mu)
  return(-grad)
}


#' Helper function to replace values in a list
#' @param x a list
#' @param list indices of list
#' @param values what values to replace with
#' @return x with list values replaced
#' @export
list_replace <- function(x, list, values){
  x[[list]] <- values; x
}
