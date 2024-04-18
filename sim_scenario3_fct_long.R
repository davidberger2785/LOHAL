#all functions for the longitudinal scenario with baseline covariates only
library(gtools)

aslongish<-function(odat,var.list,K){ #this will need to be modified for var.LIST
  n=dim(odat)[1]
  odat$id=1:n
  longdat=odat
  for (t in 2:K){
  longdat=as.data.frame(rbind(longdat,odat)) #stack a dataset per time point
  }
  longdat$t=rep(0:(K-1), each=n)
  longdat[longdat$t!=(K-1),names(longdat) %in% "Y"]=NA #sets the outcome to NA except at the last time point
  longdat$A=rep(NA,K*n)
  longdat$A[longdat$t==0]=odat$A0
  for(t in 1:(K-1)){
    longdat$A[longdat$t==t]=odat[,length(var.list)+t+1]
  }
  #vector of all A node names
  #allA=paste(rep("A",K),0:(K-1),sep="")
  #longdat=longdat[,!(names(longdat) %in% allA)]
  return(longdat)
}

getQ<-function(odat,var.LIST, K){
  n=dim(odat)[1]
  allA=paste(rep("A",K),0:(K-1),sep="")
  #enumerate all possible combinations of a-bar
  abarmat=permutations(n=2,K,v=0:1,repeats=T)
  abarbigmat=data.frame(abarmat)
  abarbigmat= abarbigmat[rep(seq_len(nrow(abarmat)), each = n), ]
  names(abarbigmat)=allA
  Q_stack=matrix(ncol=K,nrow=nrow(abarmat)*n)
  #construct a data frame, then give names
  for (a in 1:dim(abarmat)[1]){
    Qk=odat$Y
    for (k in K:1){
      Qk.full.form = formula(paste("Qk~",paste(c(var.LIST[[k]],allA[1:k]),collapse="+")))
      Qkmod<-lm(Qk.full.form,data=odat) 
      Qk=predict(Qkmod, newdata=as.data.frame(cbind(odat[,var.LIST[[k]]],abarbigmat[((a-1)*n+1):((a)*n),1:k, drop = FALSE]))) 
      Q_stack[((a-1)*n+1):((a)*n),k]=Qk
    }
  }
  return(list(Q_stack=Q_stack,abarbigmat=abarbigmat))
}

Gcomp_est<-function(L0,Q_stack,abarbigmat){
  L0_stack=rep(L0,nrow(Q_stack)/length(L0))
  cumA_stack=rowSums(abarbigmat)
  MSM_mod=lm(Q_stack[,1]~L0_stack+cumA_stack)
  return(coefficients(MSM_mod))
}

get_betas<-function(odat,Q_stack,var.LIST,abarbigmat,K){
  n=dim(odat)[1]
  pa=nrow(abarbigmat)/n
  Vars_s=NULL
  #Vars_s=do.call(rbind,replicate(pa,odat[,c(var.list)],simplify=FALSE))
  beta=NULL
  #look through time points
  for (k in K:2){
    Vars_s=do.call(rbind,replicate(pa,odat[,c(var.LIST[[k]])],simplify=FALSE))
    Vars_s=cbind(Vars_s,abarbigmat)
    outcomek = Q_stack[,k]
    beta.full.form = formula(paste("outcomek~",paste(c(var.LIST[[k]],names(abarbigmat)[1:(k-1)]),collapse="+")))
    modlm=lm(beta.full.form,data=Vars_s)
    beta=c(coefficients(modlm),beta)
  }
  #for K=1
  Vars_s=do.call(rbind,replicate(pa,odat[,c(var.LIST[[1]])],simplify=FALSE))
  Vars_s=cbind(Vars_s,abarbigmat)
  beta.full.form = formula(paste("outcomek~",paste(c(var.LIST[[1]]),collapse="+")))
  modlm=lm(beta.full.form,data=Vars_s)
  beta=c(coefficients(modlm),beta)
  return(beta)
}

#next is piprod
#generic compute propensity score function
#output is vector of length n
#fuse is a Graph or null
library(matrixStats)
piprod<-function(alpha,odat,var.LIST,K){
  #length variable set (baseline only)
  pikA=matrix(nrow=length(odat$Y),ncol=K)
    plast=0
    for (k in 1:K){
      if (k>1){cond_treats=paste("A",0:(k-2),sep="")}else{cond_treats=NULL}
      pk=length(var.LIST[[k]])
      alphak=alpha[(plast+1):(plast+1+pk+length(cond_treats))] #the alphas for each time point including intercept and past tx
      Ak=odat[,paste("A",k-1,sep="")]
      pik=t(plogis(alphak%*%t(cbind(1,odat[,c(var.LIST[[k]],cond_treats)]))))
      pikA[,k]=Ak*pik+(1-Ak)*(1-pik)
      plast=plast+1+pk+length(cond_treats)
    }
  pikAprod=rowCumprods(pikA)
  return(pikAprod)
}

#weighted mean to use in the balance criterion
#V is a matrix
#outputs a vector
wmd<-function(tx,V,w,sub=rep(1,dim(V)[1])){
  wmd_vect=function(x) {abs(sum(tx*sub*x*w)/sum(tx*sub*w)-sum((1-tx)*sub*x*w)/sum((1-tx)*sub*w))}
  apply(V,2,wmd_vect)
}

##balance criterion
#sigmas are standard errors of beta estimates
#w is an n*K matrix with columns corresponding to time points k=1,...,K
BC<-function(odat,w,betas,sigmas=rep(1,length(betas)),var.LIST,K){
  #calculate the weighted differences for every (baseline) variable at each time-point
  #betas and sigmas should only be for the covariates
  wmdk=NULL
  for (k in 1:K){
      wk=w[,k]
      #NEED TO CHECK FOR NO VARIABLES AT A TIME POINT for general case
      wmdk=c(wmdk, wmd(odat[,paste("A",k-1,sep="")],odat[,var.LIST[[k]]],wk) )
    }
  return( sum(abs(betas)*wmdk/sigmas) )
}

#generic IPTW function
IPTW_est<-function(y,w,cuma,l0){
  MSM_mod=lm(y~l0+cuma,weights=w)
  return(coefficients(MSM_mod))
}

IPTW_split<-function(odat,var.LIST,K){
  alpha=NULL
  Amat=NULL
  for (k in 1:K){
    outcome=paste(paste("A",k-1,sep=""),"~",sep="")
    if (k>1){cond_treats=paste("A",0:(k-2),sep="")}else{cond_treats=NULL}
    g.full.form=formula(paste(outcome,paste(c(var.LIST[[k]],cond_treats),collapse="+")))
    pmod<-glm(g.full.form,data=odat,family=binomial())
    alpha=c(alpha,coef(pmod))
    Amat=cbind(Amat,odat[,paste("A",k-1,sep="")])
  }
  pi=piprod(alpha,odat,var.LIST=var.LIST,K=K)[,K]
  cumA<-rowSums(Amat)
  return(IPTW_est(y=odat$Y,w=1/pi,cuma=cumA,l0=odat[,1]))
  #MSM_mod=lm(Y~L0*cumA,weights=1/(pi0*pi1),data=odat)
}

IPTW_pool<-function(longdat,odat,var.LIST,K){
  #pooled model
  Amat=NULL
  longdat[is.na(longdat)]=0
  for (k in 1:K){
    Amat=cbind(Amat,odat[,paste("A",k-1,sep="")]) #construct matrix of treatments for later
    if (k>1){
      treats=paste("A",0:(k-2),sep="")
      #cond_treats=paste(treats,paste(c("*(t==",k-1,")"),collapse=""),sep="")
      }else{treats=NULL}
    if(!identical(var.LIST[[k]], character(0))) { ##NEED TO CHECK WHAT HAPPENS WHEN NOTHING IS SELECTED. Can't be null
      if (k==1){ expk=paste("I(1*(t==0))+",paste(paste("I(",c(var.LIST[[k]]),paste(c("*(t==",k-1,")"),collapse=""),")",collapse="+"),collapse="+"),collapse="")
        }else{
          expk=paste(expk,paste("+I(1*(t==",k-1,"))+",collapse=""),paste("I(",c(var.LIST[[k]],treats),paste(c("*(t==",k-1,")"),collapse=""),")",collapse="+"),collapse="+")
        }
    } else{ #if no variables are in the list for this time point
      if (k==1){ expk="I(1*(t==0))"
      } else{ #k=2,...,K
        expk=paste(expk,paste("+I(1*(t==",k-1,"))+",collapse=""),paste("I(",treats,paste(c("*(t==",k-1,")"),collapse=""),")",collapse="+"),collapse="+")
        }
    }
  }
  g.full.form = formula(paste("A~-1+",expk,collapse="+"))

  pimod=glm(g.full.form,family=binomial(),data=longdat)
  alpha<-coefficients(pimod)
  pi=piprod(alpha=alpha,odat=odat,var.LIST=var.LIST,K=K)
  cumA<-rowSums(Amat)
  return(list(est=IPTW_est(y=odat$Y,w=1/pi[,K],cuma=cumA,l0=odat[,1]), alpha=alpha, pi=pi)) 
  #MSM_mod=lm(Y~L0*cumA,weights=1/(pi0*pi1),data=odat)
}


IPTW_pool_fuse<-function(longdat,odat,var.LIST,fuse.graph,K){
  #pooled model
  Amat=NULL
  longdat[is.na(longdat)]=0
  NC=length(fuse.graph)
  for (k in 1:K){
    Amat=cbind(Amat,odat[,paste("A",k-1,sep="")]) #construct matrix of treatments for later
    if (k>1){
      treats=paste("A",0:(k-2),sep="")
      #cond_treats=paste(treats,paste(c("*(t==",k-1,")"),collapse=""),sep="")
    }else{treats=NULL}
    if(!identical(var.LIST[[k]], character(0))) { 
      #find the variables to remove because they're in a clique
      removevar=NULL
      for (c in 1:NC){ 
        ind=which(fuse.graph[[c]][,1]==k) ;
        if(!identical(var.LIST[[k]], character(0))){ removevar=c(removevar,fuse.graph[[c]][ind,2])  } 
      }
      if (k==1){ 
        if(!identical(var.LIST[[k]][-removevar], character(0))){ #if there are variables left
        expk=paste("I(1*(t==0))+",paste(paste("I(",c(var.LIST[[k]][-removevar]),paste(c("*(t==",k-1,")"),collapse=""),")",collapse="+"),collapse="+"),collapse="")
          } else{ #if there are no variables left after removing the ones in cliques
            expk="I(1*(t==0))"
          } 
      } else{ #if k>1, need to include treatments (here, the remaining var list can be empty)
        expk=paste(expk,paste("+I(1*(t==",k-1,"))+",collapse=""),paste("I(",c(var.LIST[[k]][-removevar],treats),paste(c("*(t==",k-1,")"),collapse=""),")",collapse="+"),collapse="+")
      }
    } else{ #if no variables are in the list for this time point
      if (k==1){ expk="I(1*(t==0))"
      } else{ #k=2,...,K
        expk=paste(expk,paste("+I(1*(t==",k-1,"))+",collapse=""),paste("I(",treats,paste(c("*(t==",k-1,")"),collapse=""),")",collapse="+"),collapse="+")
      }
    }
  }
  #now we need to add the cliques back in with an indicator for any of the corresponding times
  expc=NULL
  for (c in 1:NC){
    expc = paste(expc, paste( paste("I( (", paste("t==",fuse.graph[[c]][,1]-1,collapse="|"), ")*",collapse=""), var.LIST[[fuse.graph[[c]][1,1]]][fuse.graph[[c]][1,2]] , ")", collapse="" ), sep="+")
    }
  
  g.full.form = formula(paste("A~-1+",expk,expc,collapse="+"))
  
  pimod=glm(g.full.form,family=binomial(),data=longdat)
  #do it directly because the alpha names have changed
  alpha<-coefficients(pimod)  # !!! ly add this line 
  pi=predict(pimod,type="response")
  pipr=rep(1,length=length(odat$Y))
  for (k in 1:K){
    treat=paste("A",k-1,sep="")
    pipr=pipr*( pi[longdat$t==(k-1)]*longdat[longdat$t==(k-1),treat] + (1-pi[longdat$t==(k-1)])*(1-longdat[longdat$t==(k-1),treat]) )
    
  }
  cumA<-rowSums(Amat)
  return(list(est=IPTW_est(y=odat$Y,w=1/pipr,cuma=cumA,l0=odat[,1]), alpha=alpha, pi=pipr))  # should remove 'alpha=alpha'
  #MSM_mod=lm(Y~L0*cumA,weights=1/(pi0*pi1),data=odat)
}


#Input betas
#n is sample size
IPTW_OALASSO<-function(longdat,odat,betas,var.LIST,sigmas=rep(1,length(betas)), K){

  est=rep(NULL,2)
  bcvals=rep(NULL,2) #size?
  saved.coefs=rep(NULL,length(betas))
  saved.coefvals=rep(NULL,length(betas))
  saved.alphas=NULL
  

  y=longdat$A
  n=length(y)/K
  x=NULL
  pk=rep(NA,K) # this saves length of each component
  nopen=list() #this saves indices that we don't penalize (intercepts and treatments and could add user-input here)
  ##construct design matrix to use in glmnet
  #g.full.form = formula(paste("A~-1+",paste(c("I(1-t)",paste("I(",var.list0,"*(1-t)",")",collapse="+"),"t",paste("I(",var.list0,"*t",")",collapse="+"),paste("I(",var.list1,"*t",")",collapse="+"),"I(A0*t)"),collapse="+")))
  #x=cbind(1-longdat$t,longdat[,var.list0]*(1-longdat$t),longdat$t,longdat[,var.list0]*(longdat$t),longdat[,var.list1]*(longdat$t),longdat$A0)
  for (k in 1:K){
    if (k>1){
      treats=paste("A",0:(k-2),sep="")
      pk[k]= length(c(var.LIST[[k]],treats)) +1
      nopen[[k]]=c(sum(pk[1:k-1])+1,(sum(pk[1:k-1])+1+length(var.LIST[[k]])+1):(sum(pk[1:k-1])+1+length(var.LIST[[k]])+length(treats)) )
    }else{ #k=1
      treats=NULL
      pk[1]= length(c(var.LIST[[1]])) +1
      nopen[[1]]=1
    }
    x=cbind(cbind( x, rep(1,n*K)*(longdat$t==(k-1)) ) ,longdat[,c(var.LIST[[k]],treats)]*(longdat$t==(k-1)) )
  }
  
  gamma=2.5
  
  gam=gamma
  AL_w=1/abs(betas)^gam
  AL_w[unlist(nopen)]=0 #no penalties on the selected elements #THERE IS AN ERROR HERE
  #logit_oal = glmnet( x=x,y=y,alpha=1,standardize=TRUE,  penalty.factor=AL_w, family='binomial',intercept=FALSE) #selects a vector of lambdas
  vectorlam = exp(seq(0,20,length=50))
  logit_oal = glmnet( x=x,y=y,alpha=1,standardize=TRUE,  penalty.factor=AL_w, family='binomial',intercept=FALSE,lambda=vectorlam)
  coefvals=coef(logit_oal)
  lambdas=logit_oal$lambda
  
  selectedcovs<-coef(logit_oal,s=logit_oal$lambda)!=0
  
  #if (k>1){cond_treats=paste("A",0:(k-2),sep="")}else{cond_treats=NULL}
  #pk=length(var.LIST[[k]])
  #alphak=alpha[(plast+1):(plast+1+pk+length(cond_treats))] #the alphas for each time point including intercept and past tx
  #Ak=odat[,paste("A",k-1,sep="")]
  #pik=t(plogis(alphak%*%t(cbind(1,odat[,c(var.LIST[[k]],cond_treats)]))))
  #pikA[,k]=Ak*pik+(1-Ak)*(1-pik)
  #plast=plast+1+pk+length(cond_treats)
  
  for(m in 1:(dim(selectedcovs)[2])){
    lam=lambdas[m]
    var.LIST_star=list()
    plast=1
    for (k in 1:K){
      pk=length(var.LIST[[k]])
      var.LIST_star[[k]]=var.LIST[[k]][selectedcovs[(plast+2):(plast+2+pk-1),m]]
      plast=plast+2+pk-1+k-1
    }
    selIPTW=IPTW_pool(longdat,odat,var.LIST=var.LIST_star,K=K)
    #sel.mod<-glm.fit(y=y,x=x[,selectedcovs[-1,m]],family=binomial())
    alpha=rep(0,dim(x)[2])
    alpha[selectedcovs[-1,m]]= selIPTW$alpha
    saved.alphas=rbind(saved.alphas,alpha);
    pi=selIPTW$pi #matrix of cumulative propensities
    w=1/pi #for now, could stabilize later #this is a matrix of weights
    #calculate the balance criterion given the weights at these tuning parameters
    #only send betas for the covariates -- use complete list of covariates, not the selected list
    beta_vars=betas[unlist(var.LIST)]
    sigma_vars=sigmas[unlist(var.LIST)]
    bcvals=rbind(bcvals,c(lam,BC(odat,w=w,betas=beta_vars,sigmas=sigma_vars,var.LIST=var.LIST,K=K)))
    ##output the estimated parameter values for each value of lambda
    est=rbind(est,c(lam,gam,selIPTW$est))
    
  }
  #save the selected coefficients
  saved.coefs=cbind(saved.coefs,selectedcovs[-1,]) #concatenates the binary matrix of which covariates were selected
  #plot(logit_oal)
  #plot(bcvals)
  
  #pick the value at the minimum of the balance criterion 
  ind=which.min(bcvals[,2]) 
  a.alphas=saved.alphas[ind,]; names(a.alphas)=names(x)
  return(list(est=est[ind,],selcovs=saved.coefs[,ind],covvals=a.alphas))
}


getBIC<-function(x,y,alphas,pen=log(nrow(x))){
  df=length(unique(alphas[!is.na(alphas)])) #how many non-zero coefficients = num of parameters
  lin_pred=as.matrix(x)[,!is.na(alphas)]%*%alphas[!is.na(alphas)]
  prob_pred=1/(1+exp(-lin_pred))
  Dev <- -2*sum(y*log(prob_pred) +(1-y)*log(1-prob_pred))
  BICc<- Dev + pen*df
  return(BICc)
}


library(FusedLasso)
source("C://Users/alain/Documents/Liuyan_20180731/UdeM/NDIT OALASSO/Mireille Code/fusedlasso_writeover.R") #this just replaced cBind with cbind

#alpha_pre comes from OALASSO procedure, is zero for the unselected variables, with common names for the elements that can be fused
#the inserted list.vars is the original list (not the selected one)

IPTW_fOALASSO_2steps<-function(longdat,odat, betas, sigmas, alpha_pre,var.LIST,fuse.graph,K,pk){
  est=rep(NULL,2)
  bicvals=rep(NULL,2)
  saved.coefs=rep(NULL,length(betas))
  saved.fused=rep(NULL,3)
  saved.alphas=rep(NULL,length(betas))
  Amat=NULL
  x=NULL
  nopen=list()
  for (k in 1:K){
    Amat=cbind(Amat,odat[,paste("A",k-1,sep="")]) #construct matrix of treatments for later
    if (k>1){
      treats=paste("A",0:(k-2),sep="")
      pk[k]= length(c(var.LIST[[k]],treats)) +1
      nopen[[k]]=c(sum(pk[1:k-1])+1,(sum(pk[1:k-1])+1+length(var.LIST[[k]])+1):(sum(pk[1:k-1])+1+length(var.LIST[[k]])+length(treats)) )
    }else{ #k=1
      treats=NULL
      pk[1]= length(c(var.LIST[[1]])) +1
      nopen[[1]]=1
    }
    x=cbind(cbind( x, rep(1,n*K)*(longdat$t==(k-1)) ) ,longdat[,c(var.LIST[[k]],treats)]*(longdat$t==(k-1)) )
  }
  y=longdat$A
  
  gam=2.5
  #AL_w=1/(abs(alpha_pre)^gam) #NOW WE USE STANDARD ADAPTIVE WEIGHTS -- these came from running OALASSO
  #AL_w[c(1,4,9)]=0 #no penalty placed on the two intercepts, alpha0 and alpha3, or on a0, alpha8
  AL_w=rep(0,length(alpha_pre)) #No more penalization of individual main terms
  
  sel<-alpha_pre!=0
  p=sum(sel)
  graphConn=vector("list", p) #initialize the clique (or pairs) graph
  graphWeights=vector("list", p)
  alpha_reduced=alpha_pre[sel]
  #find the cliques by common name but exclude intercepts and treatments
  varnames=unique(unlist(var.LIST))
  #commonvars=intersect(names(x[1:3])[sel0],names(x[4:9])[sel1])
  for (j in varnames){
    indices=which((names(x)[sel])==j) #indices post selection to use in the graph
    if (length(indices)>1){
      for (m in 1:length(indices)){ #for each index make a node, then connect it to others in the clique (reciprocally)
        graphConn[[indices[m]]]=as.integer(indices[-m]) #this is its clique
        for(o in indices[-m]){
          graphWeights[[indices[m]]]=c(graphWeights[[indices[m]]],1/abs(alpha_reduced[indices[m]]-alpha_reduced[o])^gam) #add on a weight for each node in clique
        }
      }
    }
  }
  
  
  Graph=list(graphConn,graphWeights)
  
  #lambda1=exp(rep(seq(0,-30,length=10),10))
  lambda1=rep(0,20)
  lambda2=sort(exp(seq(0,-20,length=20)),decreasing=T)
  
  AL_w=AL_w[sel]
  obj=fusedlasso(X=as.matrix(x[,sel]),y=y,family="binomial",
                 wLambda1=AL_w,graph=Graph,addIntercept=F,lambda2=lambda2,lambda1=lambda1)
  
  selectedcovs<-matrix(F,ncol=length(betas),nrow=length(lambda1))
  selectedcovs[,sel]= as.logical(t(obj$beta!=0)) #what is selected in the second phase
  
  #cbind(lambda1,lambda2,t(obj$beta))[1:20,]
  
  
  for (k in 1:length(lambda1)){
    
    #####use BIC
    
    #given fused coefficients and selected variables, fit the pooled tx model
    
    ind=(selectedcovs)[k,]
    
    alpha=rep(0,length(betas)); names(alpha)=names(x)
    #alpha[ind]=obj$beta[obj$beta[,k]!=0,k] #selected from the second phase
    
    x2=x
    #run through all nodes in the graph to see if fused
    f_ind=rep(0,length(graphConn))
    for (j in 1:length(graphConn)){
      if(!is.null(graphConn[[j]])){
        for (m in graphConn[[j]]){
          if(j>m){
            f_ind[j]=((obj$beta[j,k]==obj$beta[m,k])&(obj$beta[j,k]!=0)) #we indicate the latter covariate(s), to be removed; we keep the minimum index
            x2[,sel][,m]=x2[,sel][,j]*f_ind[j]+x2[,sel][,m] #IF FUSED we fuse the two covariates or intercepts and put them in the minimum index
          } 
        }
      }
    }
    f_ind_exp=rep(F,length(betas)); f_ind_exp[sel]=f_ind
    to_remove=!ind | f_ind_exp==1 #we'll remove all covariates not selected from each model #and also those that were fused to one in an earlier model
    
    sel.mod<-glm.fit(y=y,x=x2[,!to_remove],family=binomial()) 
    alpha[!to_remove]=coef(sel.mod)
    
    #now need to transform alpha back by replacing the second (& later) value(s) with the first
    if(sum(f_ind)!=0){
      for (j in which(f_ind==1)){
        alpha[sel][j]=alpha[sel][min(graphConn[[j]])]
      }}
    
    saved.alphas=rbind(saved.alphas,alpha); 
    pi=piprod(alpha,odat,var.LIST=var.LIST,K=K) #don't need to indicate which fused
    w=1/pi #for now, could stabilize later #is matrix
    #calculate the balance criterion given the weights at these tuning parameters
    bicvals=rbind(bicvals,c(lambda1[k],lambda2[k],getBIC(x=x,y=y,alphas=alpha)))
    #bicvals=rbind(bicvals,c(lambda1[k],lambda2[k],BC(odat=odat,w=w,betas=betas,sigmas=sigmas)))
    ##output the estimated parameter values for each value of lambda
    est=rbind(est,c(lambda1[k],lambda2[k],gam, IPTW_est(y=odat$Y,w=w[,K],cuma=rowSums(Amat),l0=odat[,1])))
    #saved.fused=rbind(saved.fused,f_ind)
  }
  #save the selected coefficients
  saved.coefs=rbind(saved.coefs,selectedcovs) #concatenates the binary matrix of which covariates were selected
  
  
  #pick the value at the minimum of the BIC criterion 
  ind=which.min(bicvals[,3])
  
  return(list(est=est[ind,],selcovs=saved.coefs[ind,],savedalphas=saved.alphas[ind,]))    
  
}

###run simulations

seed=read.table("C://Users/alain/Documents/Liuyan_20180731/UdeM/NDIT OALASSO/Mireille Code/seeds.txt",header=F)$V1
source("C://Users/alain/Documents/Liuyan_20180731/UdeM/NDIT OALASSO/Mireille Code/bootstrap.vect.R")
get_betas_forbs<-function(odat,var.LIST){ 
  Q=getQ(odat,var.LIST,K=K); 
  betas=get_betas(odat=odat,Q$Q_stack,var.LIST=var.LIST,abarbigmat=Q$abarbigmat,K=K); 
  return(betas) 
  }

library(glmnet)
library(data.table)
library(dplyr)
library(glmnet)
library(sandwich)
library(FusedLasso)
library(MASS)


#fuse.graph_star #which non-zero variables SHOULD be fused
var.list=gendat_long5K(n=1000,seed=seed[19])$var.list # !! changed to give initial value for var.list, modified on 09-28
p=length(var.list)
K=5
fuse.graph_star=vector("list", 2) #2 or 4 is the number of cliques. Each matrix represents the time and index of each member wrt var.LIST
fuse.graph_star[[1]]=matrix(c(1,1,2,1,3,1,4,1,5,1),nrow=5,byrow=T)
fuse.graph_star[[2]]=matrix(c(1,2,2,2,3,2,4,2,5,2),nrow=5,byrow=T)
#fuse.graph_star[[3]]=matrix(c(1,3,2,3,3,3,4,3,5,3),nrow=5,byrow=T)
#fuse.graph_star[[4]]=matrix(c(1,4,2,4,3,4,4,4,5,4),nrow=5,byrow=T)


var.LIST=list()
var.list_star=c( "Xc0_1","Xc0_2","Xp0_1","Xp0_2")
var.LIST_star=list()
pk=vector(length=K) #this is the length of each k-specific alpha vector (related to the propensity score model at k)
for (k in 1:K){var.LIST[[k]]=var.list; var.LIST_star[[k]]=var.list_star; pk[k]=length(var.LIST[[k]])+length(1:k)}



ests=matrix(nrow=1000,ncol=18)
which_covs=matrix( nrow=1000,ncol=2*(length(var.list)*K+sum(0:(K-1))+K) )
covsvals=matrix(nrow=1000,ncol=2*(length(var.list)*K+sum(0:(K-1))+K) )
whichfused=matrix(nrow=1000,ncol=20) #ignore intercepts and tx
# truevals<-c(0,0.6,0.5)
n=200
for (l in 1:1000){
    Datastuff=gendat_long5K(n,seed=seed[l])
    odat<-Datastuff$Data
    
    longdat=aslongish(odat,var.list=var.list,K=5)
    longdat[is.na(longdat)]=0
    
    Q=getQ(odat,var.LIST=var.LIST,K=K)
    ests[l,1:3]=Gcomp_est(odat[,1],Q_stack=Q$Q_stack,abarbigmat = Q$abarbigmat)
    ests[l,4:6]=IPTW_split(odat,var.LIST=var.LIST,K=K)
    
    #oracle variable select
    ests[l,7:9]=IPTW_pool(longdat,odat,var.LIST=var.LIST_star,K=K)$est #target 1 
    #oracle variable select + fuse
    ests[l,10:12]=IPTW_pool_fuse(longdat=longdat,odat=odat,var.LIST=var.LIST_star,fuse.graph=fuse.graph_star,K=K)$est
    
    #print(l)
    
    bs_betastars<-bootstrap.vect(get_betas_forbs,odat,Kk=100,var.LIST=var.LIST)
    sigs<-sqrt(bs_betastars$VAR)
    
    betas=get_betas(odat,Q$Q_stack,var.LIST,abarbigmat=Q$abarbigmat,K=K)
    iptwmod<-IPTW_OALASSO(longdat,odat,betas,var.LIST=var.LIST,sigmas=sigs,K=K)
    ests[l,13:15]=iptwmod$est[3:5]
    which_covs[l,1:115]=as.vector(iptwmod$selcovs)
    
    fiptwmod<-IPTW_fOALASSO_2steps(longdat,odat,betas=betas,sigmas=sigs,alpha_pre=as.vector(iptwmod$covvals),var.LIST=var.LIST,K=K,pk=pk)
    ests[l,16:18]=fiptwmod$est[4:6]
    which_covs[l,116:230]=as.vector(fiptwmod$selcovs)
    alphas=as.vector(fiptwmod$savedalphas)
    
    whichfused[l,]=c( #to update when you have the final design
      (alphas[1+1]==alphas[1+pk[1]+1])&(alphas[1+1]!=0), #1-2
      (alphas[1+2]==alphas[2+pk[1]+1])&(alphas[2+1]!=0),
      
      (alphas[1+1]==alphas[1+pk[1]+pk[2]+1])&(alphas[1+1]!=0), #1-3
      (alphas[1+2]==alphas[2+pk[1]+pk[2]+1])&(alphas[2+1]!=0),
      
      (alphas[1+1]==alphas[1+pk[1]+pk[2]+pk[3]+1])&(alphas[1+1]!=0), #1-4
      (alphas[1+2]==alphas[2+pk[1]+pk[2]+pk[3]+1])&(alphas[2+1]!=0),
      
      (alphas[1+1]==alphas[1+pk[1]+pk[2]+pk[3]+pk[4]+1])&(alphas[1+1]!=0), #1-5
      (alphas[1+2]==alphas[2+pk[1]+pk[2]+pk[3]+pk[4]+1])&(alphas[2+1]!=0),
     
      (alphas[1+pk[1]+pk[2]+1]==alphas[1+pk[1]+1])&(alphas[1+pk[1]+1]!=0), #2-3
      (alphas[2+pk[1]+pk[2]+1]==alphas[2+pk[1]+1])&(alphas[2+pk[1]+1]!=0),
      
      (alphas[1+pk[1]+pk[2]+pk[3]+1]==alphas[1+pk[1]+1])&(alphas[1+pk[1]+1]!=0), #2-4
      (alphas[2+pk[1]+pk[2]+pk[3]+1]==alphas[2+pk[1]+1])&(alphas[2+pk[1]+1]!=0),
     
      (alphas[1+pk[1]+pk[2]+pk[3]+pk[4]+1]==alphas[1+pk[1]+1])&(alphas[1+pk[1]+1]!=0), #2-5
      (alphas[2+pk[1]+pk[2]+pk[3]+pk[4]+1]==alphas[2+pk[1]+1])&(alphas[2+pk[1]+1]!=0),
      
      (alphas[1+pk[1]+pk[2]+pk[3]+1]==alphas[1+pk[1]+pk[2]+1])&(alphas[1+pk[1]+pk[2]+1]!=0), #3-4
      (alphas[2+pk[1]+pk[2]+pk[3]+1]==alphas[2+pk[1]+pk[2]+1])&(alphas[2+pk[1]+pk[2]+1]!=0),
      
      (alphas[1+pk[1]+pk[2]+pk[3]+pk[4]+1]==alphas[1+pk[1]+pk[2]+1])&(alphas[1+pk[1]+pk[2]+1]!=0), #3-5
      (alphas[2+pk[1]+pk[2]+pk[3]+pk[4]+1]==alphas[2+pk[1]+pk[2]+1])&(alphas[2+pk[1]+pk[2]+1]!=0),
      
      (alphas[1+pk[1]+pk[2]+pk[3]+pk[4]+1]==alphas[1+pk[1]+pk[2]+pk[3]+1])&(alphas[1+pk[1]+pk[2]+pk[3]+1]!=0), #4-5
      (alphas[2+pk[1]+pk[2]+pk[3]+pk[4]+1]==alphas[2+pk[1]+pk[2]+pk[3]+1])&(alphas[2+pk[1]+pk[2]+pk[3]+1]!=0)
      
    )
    
    write(c(l, iptwmod$est[1:2],fiptwmod$est[2:3],ests[l,]),paste0("./simulation_txt/long_fOAL_gam2_5_ests_n",n,"_upd.txt"),ncol=450,append=T)
    write(c(which_covs[l,],whichfused[l,]),paste0("./simulation_txt/long_fOAL_gam2_5_select_n",n,"_upd.txt"),ncol=450,append=T)
    print(l)

}

cm=colMeans(ests-matrix(rep(truevals,6*1000),nrow=1000,byrow=T),na.rm=T) #bias

colMeans((ests-matrix(rep(truevals,6*1000),nrow=1000,byrow=T))^2,na.rm=T) #MSE

colMeans((ests-cm)^2,na.rm=T) #var

colMeans(which_covs,na.rm=T)



# all results using 5 time points data were computed by the computer in lab
n=200 # or n=500
datest<-read.table(paste0("./simulation_txt/long_fOAL_gam2_5_ests_n",n,"_upd.txt"), header=F)[,6:23]
# n=1000  
# datest<-read.table(paste0("./simulation_txt/long_fOAL_gam2_5_ests_n",n,".txt"))[,6:23] # when n =1000, remove header=F due to different format of txt file
nr=dim(datest)[1]
truevals<-c(0, 1.14, 0.5)

rnbias=abs(round(sqrt(n)*matrix(colMeans(datest-matrix(rep(truevals,6*nr),nrow=nr,byrow=T)),nrow=6,byrow=T),1))
colm<-colMeans(datest,na.rm=T)
nmse=round(n*matrix(colMeans((datest-matrix(rep(truevals,6*nr),nrow=nr,byrow=T))^2,na.rm=T),nrow=6,byrow=T),0)

for(k in 1:6){
  print(paste(paste(paste(rnbias[k,],nmse[k,],sep="("),collapse=")&"),")",sep=""))
}

# # n=200, 
# 1] "-0.1(10)&-0.1(4)&0(1)"
# [1] "-4.6(107)&-2.1(35)&1.9(13)"
# [1] "-1.2(41)&-1.3(17)&0.5(6)"
# [1] "-1.3(41)&-1.2(16)&0.5(5)"
# [1] "-2.5(38)&-1.9(16)&1.1(5)"
# [1] "-2.5(35)&-1.8(14)&1.1(5)"
# # n=500, 
# [1] "-0.1(9)&-0.1(4)&0(1)"
# [1] "-5.5(154)&-2.5(45)&2.4(21)"
# [1] "-1(57)&-1.4(21)&0.5(8)"
# [1] "-1.1(52)&-1.4(20)&0.5(7)"
# [1] "-2.4(48)&-2.1(18)&1.1(7)"
# [1] "-2.4(47)&-2.1(18)&1.1(7)"
# # n=1000, 
# [1] "0.1(9)&0(4)&0(1)"
# [1] "-6.8(225)&-2.4(62)&2.9(32)"
# [1] "-0.5(74)&-1.6(28)&0.3(11)"
# [1] "-0.5(74)&-1.6(28)&0.3(11)"
# [1] "-2.4(55)&-2.4(24)&1.1(8)"
# [1] "-2.3(56)&-2.3(25)&1.1(9)"




#------------------------------------ figures -----------------------------------

library(ggplot2)
n=1000
selest<-read.table(paste0("./simulation_txt/long_fOAL_gam2_5_select_n",n,"_upd.txt"),header=F)

ressel=data.frame(rep(c(rep("C",2),rep("P",2),rep("I",2),rep("S",14)),5))
names(ressel)="type"
ressel$val=colMeans(selest)[c(2:21,23:42,45:64,68:87,92:111)]
ressel$name=paste(c(rep(0,length(var.list)),rep(1,length(var.list)),rep(2,length(var.list)),rep(3,length(var.list)),rep(4,length(var.list))),rep(var.list,5),sep="")

ggplot(ressel, aes(x=name, y=val, fill=type)) +
  geom_bar(stat="identity")+theme_minimal()+ 
  scale_fill_manual(values=c("C"="#4cbb17","P"= "#56B4E9","I"="red","S"= "orange"))+
  xlab("covariate") + geom_hline(yintercept=1,linetype="dashed", color = "black")+
  theme(axis.text=element_text(size=10,angle=90), axis.title=element_text(size=16),  plot.title=element_text(hjust=0.5,size=19),legend.position="right") +
  ggtitle(paste0("Proportion selection, n=", n))


# ggsave(
#   paste0("./figures/Scenlong_Select", n,"n_upd.jpeg"),
#   plot = last_plot(),
#   device = "jpeg",
#   path = NULL,
#   scale = 1,
#   width = 8,
#   height = 6,
#   units = c("in", "cm", "mm", "px"),
#   dpi = 300,
#   limitsize = TRUE
# )


resfuse=data.frame(
  c(
    paste("1-2",var.list[1:2],sep=""),
    paste("1-3",var.list[1:2],sep=""),
    paste("1-4",var.list[1:2],sep=""),
    paste("1-5",var.list[1:2],sep=""),
    paste("2-3",var.list[1:2],sep=""),
    paste("2-4",var.list[1:2],sep=""),
    paste("2-5",var.list[1:2],sep=""),
    paste("3-4",var.list[1:2],sep=""),
    paste("3-5",var.list[1:2],sep=""),
    paste("4-5",var.list[1:2],sep="")
  )
)
names(resfuse)="covariate"
resfuse$val=colMeans( selest[,231:250] )
resfuse$type=rep(c("C","C","C","C"), 5)
ggplot(resfuse, aes(x=covariate, y=val, fill=type)) +
  geom_bar(stat="identity")+theme_minimal()+ 
  scale_fill_manual(values=c("C"="#4cbb17","P"= "#56B4E9","I"="red","S"= "orange"))+
  xlab("covariate") + geom_hline(yintercept=1,linetype="dashed", color = "black")+
  theme(axis.text=element_text(size=10,angle=90), axis.title=element_text(size=16),  plot.title=element_text(hjust=0.5,size=19),legend.position="right") +
  ggtitle(paste0("Proportion fused and non-zero, n=",n))

# ggsave(
#   paste0("./figures/Scenlong_Fuse",n,"n_upd.jpeg"),
#   plot = last_plot(),
#   device = "jpeg",
#   path = NULL,
#   scale = 1,
#   width = 8,
#   height = 6,
#   units = c("in", "cm", "mm", "px"),
#   dpi = 300,
#   limitsize = TRUE
# )

