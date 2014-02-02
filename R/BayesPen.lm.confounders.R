BayesPen.lm.confounders <-
function(y,x,u,prior, nIter=500, burnIn=100, thin=1, update, force, max.steps=NULL, max.refit, saveAt="", include.me=FALSE, z.score=FALSE){
    
  if(missing(prior)) prior<- NULL
  if(missing(update)) update <- round(500/10)
  if(missing(force)) force <- NULL
  if(missing(max.refit)) max.refit <- NULL
  
  if (nIter < ncol(x)){
    stop("The number of MCMC draws must be set to be larger than the dimension for cronfounder selection.")
  }
  
  x <- as.matrix(x)
  fit.out<-BLR2(y, XR=cbind(x,u), prior=prior, nIter = nIter, burnIn = burnIn, thin = thin, update=update, saveAt=saveAt)
  if(is.null(dim(x)) | ncol(x)==1){
    confounder.weights <- c(0,abs(BLR2(x, XR=u, prior=prior, nIter = nIter, burnIn = burnIn, thin = thin, update=update, saveAt=saveAt)$bR))
    force <- c(1,force+1)
  }else{
    confounder.weights<-rep(0,ncol(u)+ncol(x))
    for(i in 1:ncol(x)){
      if(include.me){
        fit.temp <- BLR2(x[,i], XR=cbind(x[,-i],u), prior=prior, nIter = nIter, burnIn = burnIn, thin = thin, update=update, saveAt=saveAt)
        weight.temp <- abs(fit.temp$bR)
        if(z.score) weight.temp <- weight.temp/sqrt(diag(fit.temp$COV.bR))*sqrt(diag(fit.out$COV.bR))[-1]
        confounder.weights <- confounder.weights+c(rep(0,ncol(x)),weight.temp[-c(1:c(ncol(x)-1))])
      }else{
        fit.temp <- BLR2(x[,i], XR=u, prior=prior, nIter = nIter, burnIn = burnIn, thin = thin, update=update, saveAt=saveAt)
        weight.temp <- abs(fit.temp$bR)
        if(z.score) weight.temp <- weight.temp/sqrt(diag(fit.temp$COV.bR))*sqrt(diag(fit.out$COV.bR)[-c(1:3)])
        confounder.weights <- confounder.weights+c(rep(0,ncol(x)),weight.temp)
      }
    }
    force <- c(1:ncol(x),force+ncol(x))
  }
  print("Model fitting complete, start post processing.")
  print("------------------------------------------------------------")

  fit <- BayesPen(beta=fit.out$bR, beta_cov=fit.out$COV.bR,joint=TRUE,confounder.weights=confounder.weights,force=force,max.steps=max.steps)
  
  print("Refit Model")
  print("------------------------------------------------------------")
  refit <- BayesPen.refit(y,cbind(x,u),fit=fit, max.refit=max.refit)

  print("Complete")
  print("------------------------------------------------------------")
  out <- c(fit,refit[names(refit)[-grep("joint",names(refit))]])
  out$lm <- fit.out
  out$confounder.weights <- confounder.weights
  return(out)
  
}
