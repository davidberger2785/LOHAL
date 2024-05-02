

#' This function computes the residuals that are used as outcomes in the reduced dimension regressions for g.

residG <- function( A0, A1, g0n, g1n, abar, ... ){
  rg0 <- (as.numeric(A0==abar[1])-g0n)/g0n
  rg1 <- as.numeric(A0==abar[1])/g0n * (as.numeric(A1==abar[2]) - g1n)/g1n
  return(list( rg0 = rg0, rg1 = rg1  ))
}




#' estimategr: A function used to estimate the reduced dimension regressions for g.

estimategr <- function( rg0, rg1, A0, A1, folds, validFold, Q2n, Q1n, SL.gr="SL.hal9001", abar, return.models, tolg, verbose, ... ){
  all1 <- all(folds ==1)
  if(all1){
    train_Q1n <- valid_Q1n <- Q1n
    train_Q2n <- valid_Q2n <- Q2n
    train_A0 <- valid_A0 <- A0
    train_A1 <- valid_A1 <- A1
  }else{
    # training data
    train_Q1n <- Q1n[folds != validFold]
    train_Q2n <- Q2n[folds != validFold]
    train_A0 <- A0[folds != validFold]
    train_A1 <- A1[folds != validFold]

    # validation data
    valid_A0 <- A0[folds == validFold]
    valid_A1 <- A1[folds == validFold]
    valid_Q1n <- Q1n[folds == validFold]
    valid_Q2n <- Q2n[folds == validFold]
  }

  #-------------
  # g0nr: A0 ~ Q1n
  #-------------
  g0rmod <- do.call(SL.gr, args=list(
      Y=as.numeric(train_A0==abar[1]),
      X=data.frame(Q1n = train_Q1n),
      newX=data.frame(Q1n = Q1n), # need predictions on full data
      SL.library= SL.gr,
      obsWeights = rep(1, length(train_Q1n)),
      family = binomial(),
      method = "method.CC_nloglik_mod",
      cvControl = list(V = 2),
      verbose=verbose))

  g0nr <- g0rmod$pred[folds == validFold]
  if(!all1){ train_g0nr <- g0rmod$pred[folds != validFold] }else{ train_g0nr <- g0nr }

  g0nr[g0nr < tolg] <- tolg
  train_g0nr[train_g0nr < tolg] <- tolg

  #----------------
  # g1nr: A0*A1 ~ Q2n
  #----------------
  g1rmod <- do.call( SL.gr,args=list(
      Y=as.numeric(train_A0==abar[1] & train_A1==abar[2]),
      X=data.frame(Q2n = train_Q2n),
      newX = data.frame(Q2n = valid_Q2n),
      SL.library=SL.gr,
      obsWeights = rep(1, length(train_Q2n)),
      family = binomial(),
      method = "method.CC_nloglik_mod",
      cvControl = list(V = 2),
      verbose=verbose))

  g1nr <- g1rmod$pred

  # trim small values
  g1nr[g1nr < tolg] <- tolg

  #-------------------------------
  # h0nr: rg0 [= (A0 - g0n)/g0n] ~ Q1n
  #-------------------------------
  h0rmod <- do.call(SL.gr,args=list(
      Y = rg0,
      X = data.frame(Q1n = train_Q1n),
      newX = data.frame(Q1n = Q1n), # need predictions on full data here
      SL.library = SL.gr,
      obsWeights = rep(1, length(train_Q1n)),
      family = gaussian(),
      method = "method.CC_LS_mod",
      cvControl = list(V = 2),
      verbose = verbose))

  h0nr <- h0rmod$pred[folds == validFold]
  if(!all1){train_h0nr <- h0rmod$pred[folds != validFold] }else{ train_h0nr <- h0nr }

  #---------------------------------------
  # h1rn: rg1 [= A0/g0n * (A1 - g1n)/g1n] ~ Q2n
  #---------------------------------------
  h1rmod <- do.call( SL.gr, args=list(
      Y=rg1,
      X=data.frame(Q2n = train_Q2n),
      newX = data.frame(Q2n = valid_Q2n),
      SL.library=SL.gr,
      obsWeights = rep(1, length(train_Q2n)),
      family = gaussian(),
      method = "method.CC_LS_mod",
      cvControl = list(V = 2),
      verbose=verbose))

  h1nr <- h1rmod$pred

  #--------------------------------------
  # hbarnr: A0/g0nr * h0nr ~ Q2n
  #--------------------------------------
  hbarrmod <- do.call( SL.gr, args=list(
      Y = train_A0 / train_g0nr * train_h0nr,
      X = data.frame(Q2n = train_Q2n),
      newX = data.frame(Q2n = valid_Q2n),
      SL.library = SL.gr,
      obsWeights = rep(1, length(train_Q2n)),
      family = gaussian(),
      method = "method.CC_LS_mod",
      cvControl = list(V = 2),
      verbose = verbose))

  hbarnr <- hbarrmod$pred


  #--------
  # return
  #--------
  out <- list(g0nr = g0nr, g1nr = g1nr, h0nr = h0nr, h1nr = h1nr, hbarnr = hbarnr, g0rmod = NULL, g1rmod = NULL,
              h0rmod = NULL, h1rmod = NULL, hbarrmod = NULL)

  if(return.models){
    out$g0rmod <- g0rmod; out$g1rmod <- g1rmod
    out$h0rmod <- h0rmod; out$h1rmod <- h1rmod; out$hbarrmod <- hbarrmod
  }

  return(out)
}



#' reduced-dimension regression on G
#' Function that computes the reduced dimension regressions needed for the extra targeting steps of the tmle.

redReg_gr_2ts <- function(L2, A0, A1, Qn, gn, folds, validFold, abar, verbose, return.models = FALSE,
                          SL.gr="SL.hal9001", tolg, update = FALSE, gnr = NULL, ...){


    if(!update){  ## --- update = FALSE

        full_gn <- gn; full_Qn <- Qn

        if(all(folds == 1)){
          train_L2 <- valid_L2 <- L2
          train_A0 <- valid_A0 <- A0
          train_A1 <- valid_A1 <- A1
          train_Q1n <- valid_Q1n <- full_Qn$Q1n
          train_Q2n <- valid_Q2n <- full_Qn$Q2n
          train_g0n <- valid_g0n <- full_gn$g0n
          train_g1n <- valid_g1n <- full_gn$g1n
        }else{
          # training data
          train_L2 <- L2[folds != validFold]
          train_A0 <- A0[folds != validFold]
          train_A1 <- A1[folds != validFold]
          train_Q1n <- full_Qn$Q1n[folds != validFold]
          train_Q2n <- full_Qn$Q2n[folds != validFold]
          train_g0n <- full_gn$g0n[folds != validFold]
          train_g1n <- full_gn$g1n[folds != validFold]

          # validation data
          valid_L2 <- L2[folds == validFold]
          valid_A0 <- A0[folds == validFold]
          valid_A1 <- A1[folds == validFold]
          valid_Q1n <- full_Qn$Q1n[folds == validFold]
          valid_Q2n <- full_Qn$Q2n[folds == validFold]
          valid_g0n <- full_gn$g0n[folds == validFold]
          valid_g1n <- full_gn$g1n[folds == validFold]
        }
        #---------------------------------------
        # residual outcomes for gr regressions
        #---------------------------------------
        rg <- residG( A0 = train_A0, A1 = train_A1, g0n = train_g0n, g1n = train_g1n, abar = abar )

        #---------------------------
        # estimate gr regressions
        #---------------------------
        gnr <- estimategr(
            rg0 = rg$rg0, rg1 = rg$rg1, folds = folds, validFold = validFold,
            A0 = A0, A1 = A1, Q2n = full_Qn$Q2n, Q1n = full_Qn$Q1n, SL.gr = SL.gr, abar = abar,
            return.models = return.models, tolg = tolg, verbose = verbose )

        return( gnr )

    }else{ ## !--- update = TRUE
        out <- predict_gnr(g0n = gn$g0n, g1n = gn$g1n, A0 = A0, A1 = A1, L2 = L2,
                           SL.gr = SL.gr, abar = abar, verbose = verbose, gnr = gnr, tolg = tolg)
        return(out)
    }
}
