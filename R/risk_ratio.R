#' riskratio
#' @param x either data.frame of class mpti_weib_dat or parameter list of class mpti_weib_parm
#' @param ... additional arguments to class specific function calls
#' @export

riskratio <- function(x,...){
 UseMethod("riskratio")
}

#' Estimate survival risk ratios and mediation proportion for binary mediator using Cox proportional hazards
#  
#  VanderWeele T. Explanation in causal inference: methods for mediation and interaction. Oxford University Press, 2015.
#  
#  Vandenberghe S, Duchateau L, Slaets L et al. Surrogate marker analysis in cancer clinical trials through time-
#  to-event mediation techniques. Statistical methods in medical research 2018; 27(11): 3367-3385.
#  
#  Zhou J, Jiang X, Xia HA, Wei P, Hobbs BP (2021). A Survival Mediation Model with Bayesian Model Averaging, 
#  Statistical Methods in Medical Research, 30(11):2413-2427. doi:10.1177/09622802211037069
#
#  Zhou J, Jiang X, Xia HA, Wei P, Hobbs BP (2022). Predicting outcomes of phase III oncology trials with Bayesian 
#  mediation modeling of tumor response, Statistics in Medicine, 41(4):751-768. doi:10.1002/sim.9268
# 
#'
#' @param x mstate_data object
#' @param tseq times to evaluate risk ratios
#' @importFrom stats as.formula glm formula model.matrix model.frame stepfun coef
#' @importFrom survival Surv survfit coxph
#' @return list of risk ratio calculations at observed event times
#' @export

riskratio.mstate_data <- function(x, tseq=NULL){

 d <- x@data
 class(d) <- "data.frame"
 covars <- NULL
 if( !(length(x@covar)==1 && is.na(x@covar)) ){
  covars <- x@covar
 }

 out <- riskratio(x=d, y="OS", event="OSevent", trt="trt", mediator="resp", covars=covars, tseq=tseq)

 return(out)
}

#' Estimate survival risk ratios and mediation proportion for binary mediator using Cox proportional hazards
#  
#  VanderWeele T. Explanation in causal inference: methods for mediation and interaction. Oxford University Press, 2015.
#  
#  Vandenberghe S, Duchateau L, Slaets L et al. Surrogate marker analysis in cancer clinical trials through time-
#  to-event mediation techniques. Statistical methods in medical research 2018; 27(11): 3367-3385.
#  
#  Zhou J, Jiang X, Xia HA, Wei P, Hobbs BP (2021). A Survival Mediation Model with Bayesian Model Averaging, 
#  Statistical Methods in Medical Research, 30(11):2413-2427. doi:10.1177/09622802211037069
#
#  Zhou J, Jiang X, Xia HA, Wei P, Hobbs BP (2022). Predicting outcomes of phase III oncology trials with Bayesian 
#  mediation modeling of tumor response, Statistics in Medicine, 41(4):751-768. doi:10.1002/sim.9268
# 
#'
#' @param x imstate_data object
#' @param tseq times to evaluate risk ratios
#' @importFrom stats as.formula glm formula model.matrix model.frame stepfun coef
#' @importFrom survival Surv survfit coxph
#' @return list of risk ratio calculations at observed event times
#' @export

### check ###
riskratio.imstate_data <- function(x, tseq=NULL){

 d <- x@data
 class(d) <- "data.frame"
 covars <- NULL
 if( !(length(x@covar)==1 && is.na(x@covar)) ){
  covars <- x@covar
 }

 out <- riskratio(x=d, y="OS", event="OSevent", trt="trt", mediator="resp", covars=covars, tseq=tseq)

 return(out)
}


#' Estimate survival risk ratios and mediation proportion for binary mediator using Cox proportional hazards
#  
#  VanderWeele T. Explanation in causal inference: methods for mediation and interaction. Oxford University Press, 2015.
#  
#  Vandenberghe S, Duchateau L, Slaets L et al. Surrogate marker analysis in cancer clinical trials through time-
#  to-event mediation techniques. Statistical methods in medical research 2018; 27(11): 3367-3385.
#  
#  Zhou J, Jiang X, Xia HA, Wei P, Hobbs BP (2021). A Survival Mediation Model with Bayesian Model Averaging, 
#  Statistical Methods in Medical Research, 30(11):2413-2427. doi:10.1177/09622802211037069
#
#  Zhou J, Jiang X, Xia HA, Wei P, Hobbs BP (2022). Predicting outcomes of phase III oncology trials with Bayesian 
#  mediation modeling of tumor response, Statistics in Medicine, 41(4):751-768. doi:10.1002/sim.9268
# 
#'
#' @param x data frame
#' @param y name of TTF duration
#' @param event name of TTF event
#' @param trt treatment name
#' @param mediator mediator var name (surrogate)
#' @param covars names of prognostic covariates
#' @param tseq times to evaluate risk ratios
#' @importFrom stats as.formula glm formula model.matrix model.frame stepfun coef
#' @importFrom survival Surv survfit coxph
#' @return list of risk ratio calculations at observed event times
#' @export

riskratio.data.frame <-function(x, y="OS", event="OSevent", trt="trt", mediator="resp", covars=NULL, tseq=NULL ){

 if( !is.data.frame(x) ){ stop("x must be a data.frame") }
 if( !(y %in% names(x)) ){ stop("y not found in x") }
 if( !(event %in% names(x)) ){ stop("event not found in x") }
 if( !(trt %in% names(x)) ){ stop("trt not found in x") }
 if( !(mediator %in% names(x)) ){ stop("mediator not found in x") }
 if( length(covars)>0 ){ 
  if( sum(!(covars %in% names(x)))>0 ){ stop("covariate not found in x") }
 }

 if( !is.factor(x[,trt]) ){ x[,trt] <- factor(x[,trt]) }
 if( !is.factor(x[,mediator]) ){ x[,mediator] <- factor(x[,mediator]) }

 bsx <- base::subset( x, c(TRUE,rep(FALSE,nrow(x)-1)))[,c(trt,mediator,covars)]
 bsx[1,trt] <- levels(x[,trt])[[1]]
 bsx[1,mediator] <- levels(x[,mediator])[[1]]

 formu <- paste0("survival::Surv(",y,",",event,") ~ ",trt," + ",mediator)
 if(length(covars) > 0){
  for(g in covars){
   formu <- paste(formu, " + ", g, sep = "")
   if( is.factor(x[,g]) ){
    bsx[1,g] <- levels(x[,g])[[1]]
   } else{
    bsx[1,g] <- mean(x[,g],na.rm=TRUE)
   }
  }
 }

 coxfit <- survival::coxph( stats::as.formula(formu), data=x)

 if(is.null(tseq)){ tseq<-survival::survfit(coxfit)$time }
 
 rawdat<-moddat<-stats::model.frame(coxfit)
 mod.mat<-stats::model.matrix(coxfit)

 newdat<-rawdat[1,-1]
 newdat[1,]<-bsx

 var.ord<-c(grep(trt,colnames(newdat)),grep(mediator,colnames(newdat)))
 if( length(covars)>0 ){ for(g in covars){
  var.ord <- c(var.ord, grep(g,colnames(newdat)))
 }}

 colnames(newdat)[var.ord]<-c(trt,mediator,covars)
 colnames(rawdat)[-1][var.ord]<-c(trt,mediator,covars)

 bs.surv<-survival::survfit(coxfit,newdata=newdat)
 Ht<-stats::stepfun(bs.surv$time,c(0,-log(bs.surv$surv)))(tseq)

 St.pred<-exp(-Ht%o%as.numeric(exp(mod.mat%*%coef(coxfit))))
  
 St.a0m0<-apply(St.pred[,rawdat[,trt]==levels(x[,trt])[[1]]],1,mean)
 St.a1m1<-apply(St.pred[,rawdat[,trt]==levels(x[,trt])[[2]]],1,mean)
  
 moddat[,grep(trt,colnames(moddat))]<-1
 mod.mat1<-stats::model.matrix(stats::formula(coxfit),moddat)[,-1]
 St.med<-exp(-Ht%o%as.numeric(exp(mod.mat1%*%stats::coef(coxfit))))
 St.a1m0<-apply(St.med[,rawdat[,trt]==levels(x[,trt])[[1]]],1,mean)

 #### Alternative calculation to adjust for effects of prognostic covariates ####
 if( FALSE ){
  datr<-rawdat
  datr[,1] <- NULL

  St.a1m0.LE<-numeric()
  formu <- paste0("St ~ 1 ")
  if(length(covars) > 0){
   for(g in covars){ formu <- paste(formu, " + ", g, sep = "") }
  }

  for(j in 1:length(tseq)){
    datr$St<-St.med[j,]
    w <- datr[,trt]==levels(x[,trt])[[1]]
    datrs <- base::subset(datr,w)
    fit.logit<-stats::glm(stats::as.formula(formu),data=datrs,family=quasibinomial(link = "logit"))
    St.a1m0.LE[j]<-mean(predict(fit.logit,datr,type="response"))
  }
 }
 ################################################################################

 # total treatment effect
 RR_tot<-St.a1m1/(St.a0m0*as.numeric(St.a0m0>0)+(1-as.numeric(St.a0m0>0)))
 DIF_tot<-log(St.a1m1)-log(St.a0m0) 

 # direct effect of Trt on risk ratio
 RR_d<-St.a1m0/(St.a0m0*as.numeric(St.a0m0>0)+(1-as.numeric(St.a0m0>0)))
 #RR_d.LE<-St.a1m0.LE/(St.a0m0*as.numeric(St.a0m0>0)+(1-as.numeric(St.a0m0>0)))
 DIF_d<-log(St.a1m0)-log(St.a0m0)

 # indirect effect of Trt on risk ratio
 RR_m<-St.a1m1/(St.a1m0*as.numeric(St.a1m0>0)+(1-as.numeric(St.a1m0>0)))
 #RR_m.LE<-St.a1m1/(St.a1m0.LE*as.numeric(St.a1m0.LE>0)+(1-as.numeric(St.a1m0.LE>0)))
 DIF_m<-log(St.a1m1)-log(St.a1m0)

 # mediation proportion
 PROP_m<-(St.a1m1-St.a1m0)/(St.a1m1-St.a0m0)
 #PROP_m.LE<-(St.a1m1-St.a1m0.LE)/(St.a1m1-St.a0m0)

 #out <- as.list(data.frame(tseq,St.a0m0,St.a1m1,St.a1m0,St.a1m0.LE,RR_d,RR_d.LE,RR_m,RR_m.LE,RR_tot,DIF_d,DIF_m,DIF_tot,PROP_m,PROP_m.LE))
 out <- as.list(data.frame(tseq,St.a0m0,St.a1m1,St.a1m0,RR_d,RR_m,RR_tot,DIF_d,DIF_m,DIF_tot,PROP_m))

 attributes(out)[["scenario"]] <- attributes(x)$scenario

 lctr <- x[,trt]==levels(x[,trt])[[1]]
 x_ctr <- base::subset(x,lctr)
 x_trt <- base::subset(x,!lctr)
 x_ctr <- base::subset(x_ctr,!is.na(x_ctr[,mediator]))
 x_trt <- base::subset(x_trt,!is.na(x_trt[,mediator]))

 attributes(out)[["beta.ctr"]] <- c(1,1)
 attributes(out)[["beta.trt"]] <- c(1,1)

 if( nrow(x_ctr)>0 ){
  attributes(out)[["beta.ctr"]] <- c(1+sum(x_ctr[,mediator]==1), 1+nrow(x_ctr)-sum(x_ctr[,mediator]==1))
 }
 if( nrow(x_trt)>0 ){
  attributes(out)[["beta.trt"]] <- c(1+sum(x_trt[,mediator]==1), 1+nrow(x_trt)-sum(x_trt[,mediator]==1))
 }

 attributes(out)[["summary"]] <- c( "median.ctr"=round(tseq[which.min(abs(St.a0m0-0.5))],2),
                                    "FUP.ctr"=median( x_ctr[,y] ), 
                                    "FUP.ctr_lb"=summary( x_ctr[,y] )[[2]],
                                    "FUP.ctr_ub"=summary( x_ctr[,y] )[[5]],
                                    "median.trt"=round(tseq[which.min(abs(St.a1m1-0.5))],2),
                                    "FUP.trt"=median( x_trt[,y] ), 
                                    "FUP.trt_lb"=summary( x_trt[,y] )[[2]],
                                    "FUP.trt_ub"=summary( x_trt[,y] )[[5]] )
 
 class(out) <- c("riskratio_est",class(out))
 
 return(out)
}


#' Estimate survival risk ratios and mediation proportion for binary mediator using Cox proportional hazards
#  
#'
#' @param x data frame
#' @param y name of TTF duration
#' @param event name of TTF event
#' @param trt treatment name
#' @param mediator mediator var name (surrogate)
#' @param covars names of prognostic covariates
#' @param tseq times to evaluate risk ratios
#' @importFrom stats as.formula glm formula model.matrix model.frame stepfun coef
#' @importFrom survival Surv survfit coxph
#' @return list of risk ratio calculations at observed event times
#' @export

riskratio.mpti_weib_dat <-function(x, y="OS", event="OSevent", trt="trt", mediator="resp", covars=NULL, tseq=NULL ){

 if( !is.data.frame(x) ){ stop("x must be a data.frame") }
 if( !(y %in% names(x)) ){ stop("y not found in x") }
 if( !(event %in% names(x)) ){ stop("event not found in x") }
 if( !(trt %in% names(x)) ){ stop("trt not found in x") }
 if( !(mediator %in% names(x)) ){ stop("mediator not found in x") }
 if( length(covars)>0 ){ 
  if( sum(!(covars %in% names(x)))>0 ){ stop("covariate not found in x") }
 }

 if( !is.factor(x[,trt]) ){ x[,trt] <- factor(x[,trt]) }
 if( !is.factor(x[,mediator]) ){ x[,mediator] <- factor(x[,mediator]) }

 bsx <- base::subset( x, c(TRUE,rep(FALSE,nrow(x)-1)))[,c(trt,mediator,covars)]
 bsx[1,trt] <- levels(x[,trt])[[1]]
 bsx[1,mediator] <- levels(x[,mediator])[[1]]

 formu <- paste0("survival::Surv(",y,",",event,") ~ ",trt," + ",mediator)
 if(length(covars) > 0){
  for(g in covars){
   formu <- paste(formu, " + ", g, sep = "")
   if( is.factor(x[,g]) ){
    bsx[1,g] <- levels(x[,g])[[1]]
   } else{
    bsx[1,g] <- mean(x[,g],na.rm=TRUE)
   }
  }
 }

 coxfit <- survival::coxph( stats::as.formula(formu), data=x)

 if(is.null(tseq)){ tseq<-survival::survfit(coxfit)$time }
 
 rawdat<-moddat<-stats::model.frame(coxfit)
 mod.mat<-stats::model.matrix(coxfit)

 newdat<-rawdat[1,-1]
 newdat[1,]<-bsx

 var.ord<-c(grep(trt,colnames(newdat)),grep(mediator,colnames(newdat)))
 if( length(covars)>0 ){ for(g in covars){
  var.ord <- c(var.ord, grep(g,colnames(newdat)))
 }}

 colnames(newdat)[var.ord]<-c(trt,mediator,covars)
 colnames(rawdat)[-1][var.ord]<-c(trt,mediator,covars)

 bs.surv<-survival::survfit(coxfit,newdata=newdat)
 Ht<-stats::stepfun(bs.surv$time,c(0,-log(bs.surv$surv)))(tseq)

 St.pred<-exp(-Ht%o%as.numeric(exp(mod.mat%*%coef(coxfit))))
  
 St.a0m0<-apply(St.pred[,rawdat[,trt]==levels(x[,trt])[[1]]],1,mean)
 St.a1m1<-apply(St.pred[,rawdat[,trt]==levels(x[,trt])[[2]]],1,mean)
  
 moddat[,grep(trt,colnames(moddat))]<-1
 mod.mat1<-stats::model.matrix(stats::formula(coxfit),moddat)[,-1]
 St.med<-exp(-Ht%o%as.numeric(exp(mod.mat1%*%stats::coef(coxfit))))
 St.a1m0<-apply(St.med[,rawdat[,trt]==levels(x[,trt])[[1]]],1,mean)

 #### Alternative calculation to adjust for effects of prognostic covariates ####
 if( FALSE ){
  datr<-rawdat
  datr[,1] <- NULL

  St.a1m0.LE<-numeric()
  formu <- paste0("St ~ 1 ")
  if(length(covars) > 0){
   for(g in covars){ formu <- paste(formu, " + ", g, sep = "") }
  }

  for(j in 1:length(tseq)){
    datr$St<-St.med[j,]
    w <- datr[,trt]==levels(x[,trt])[[1]]
    datrs <- base::subset(datr,w)
    fit.logit<-stats::glm(stats::as.formula(formu),data=datrs,family=quasibinomial(link = "logit"))
    St.a1m0.LE[j]<-mean(predict(fit.logit,datr,type="response"))
  }
 }
 ################################################################################

 # total treatment effect
 RR_tot<-St.a1m1/(St.a0m0*as.numeric(St.a0m0>0)+(1-as.numeric(St.a0m0>0)))
 DIF_tot<-log(St.a1m1)-log(St.a0m0) 

 # direct effect of Trt on risk ratio
 RR_d<-St.a1m0/(St.a0m0*as.numeric(St.a0m0>0)+(1-as.numeric(St.a0m0>0)))
 #RR_d.LE<-St.a1m0.LE/(St.a0m0*as.numeric(St.a0m0>0)+(1-as.numeric(St.a0m0>0)))
 DIF_d<-log(St.a1m0)-log(St.a0m0)

 # indirect effect of Trt on risk ratio
 RR_m<-St.a1m1/(St.a1m0*as.numeric(St.a1m0>0)+(1-as.numeric(St.a1m0>0)))
 #RR_m.LE<-St.a1m1/(St.a1m0.LE*as.numeric(St.a1m0.LE>0)+(1-as.numeric(St.a1m0.LE>0)))
 DIF_m<-log(St.a1m1)-log(St.a1m0)

 # mediation proportion
 PROP_m<-(St.a1m1-St.a1m0)/(St.a1m1-St.a0m0)
 #PROP_m.LE<-(St.a1m1-St.a1m0.LE)/(St.a1m1-St.a0m0)

 #out <- as.list(data.frame(tseq,St.a0m0,St.a1m1,St.a1m0,St.a1m0.LE,RR_d,RR_d.LE,RR_m,RR_m.LE,RR_tot,DIF_d,DIF_m,DIF_tot,PROP_m,PROP_m.LE))
 out <- as.list(data.frame(tseq,St.a0m0,St.a1m1,St.a1m0,RR_d,RR_m,RR_tot,DIF_d,DIF_m,DIF_tot,PROP_m))

 attributes(out)[["scenario"]] <- attributes(x)$scenario

 lctr <- x[,trt]==levels(x[,trt])[[1]]
 x_ctr <- base::subset(x,lctr)
 x_trt <- base::subset(x,!lctr)
 x_ctr <- base::subset(x_ctr,!is.na(x_ctr[,mediator]))
 x_trt <- base::subset(x_trt,!is.na(x_trt[,mediator])) 

 attributes(out)[["beta.ctr"]] <- c(1,1)
 attributes(out)[["beta.trt"]] <- c(1,1)

 if( nrow(x_ctr)>0 ){
  attributes(out)[["beta.ctr"]] <- c(1+sum(x_ctr[,mediator]==1), 1+nrow(x_ctr)-sum(x_ctr[,mediator]==1))
 }
 if( nrow(x_trt)>0 ){
  attributes(out)[["beta.trt"]] <- c(1+sum(x_trt[,mediator]==1), 1+nrow(x_trt)-sum(x_trt[,mediator]==1))
 }

 attributes(out)[["summary"]] <- c( "median.ctr"=round(tseq[which.min(abs(St.a0m0-0.5))],2),
                                    "FUP.ctr"=summary( x_ctr[,y] )[[3]], 
                                    "FUP.ctr_lb"=summary( x_ctr[,y] )[[2]],
                                    "FUP.ctr_ub"=summary( x_ctr[,y] )[[5]],
                                    "median.trt"=round(tseq[which.min(abs(St.a1m1-0.5))],2),
                                    "FUP.trt"=summary( x_trt[,y] )[[3]], 
                                    "FUP.trt_lb"=summary( x_trt[,y] )[[2]],
                                    "FUP.trt_ub"=summary( x_trt[,y] )[[5]] )

 class(out) <- c("riskratio_est",class(out))
 
 return(out)
}


#' Calculate true survival risk ratios from MPTI Weibull multi-state mediation model
#'
#' @param x MPTI Weibull model parameters
#' @param maxt time duration upper bound
#' @param tseq time grid 
#' @param npart length out of time partition
#' @param nsims number of iterations for simulating response model
#' @importFrom stats rbinom
#' @return list of risk ratio calculations over time grid
#' @export

riskratio.mpti_weib_parm <-function(x,maxt=30,tseq=NULL,npart=300,nsims=1000){

 if(is.null(tseq)){ tseq <- seq(0,maxt,length.out=npart) }

 for(ii in 1:length(x$beta)){ x$beta[[ii]] <- x$beta[[ii]][1:2] }

 ### simulate with gamma ####
 pr.ctr <- exp(x$gamma[1])/(1+exp(x$gamma[1]))
 pr.trt <- exp(x$gamma[1]+x$gamma[2])/(1+exp(x$gamma[1]+x$gamma[2]))

 St.00 <- sapply(tseq, FUN=s_mpti_weib, z=c(0,0), weib_nv=x$weib_nv, weib_lamb=x$weib_lamb, beta=x$beta)
 St.01 <- sapply(tseq, FUN=s_mpti_weib, z=c(0,1), weib_nv=x$weib_nv, weib_lamb=x$weib_lamb, beta=x$beta)
 St.11 <- sapply(tseq, FUN=s_mpti_weib, z=c(1,1), weib_nv=x$weib_nv, weib_lamb=x$weib_lamb, beta=x$beta)
 St.10 <- sapply(tseq, FUN=s_mpti_weib, z=c(1,0), weib_nv=x$weib_nv, weib_lamb=x$weib_lamb, beta=x$beta)

 St.a0m0 <- St.a1m0 <-  St.a1m1 <- St.00

 for(ii in 1:nsims){
  if( stats::rbinom(1,1,pr.ctr)==1 ){
   St.a0m0 <- rbind(St.a0m0, St.01)
   St.a1m0 <- rbind(St.a1m0, St.11)
  } else{
   St.a0m0 <- rbind(St.a0m0, St.00)
   St.a1m0 <- rbind(St.a1m0, St.10)
  }
  if( stats::rbinom(1,1,pr.trt)==1 ){
   St.a1m1 <- rbind(St.a1m1, St.11)
  } else{
   St.a1m1 <- rbind(St.a1m1, St.10)
  }
 }

 St.a0m0 <- colMeans(St.a0m0[-1,])
 St.a1m0 <- colMeans(St.a1m0[-1,])
 St.a1m1 <- colMeans(St.a1m1[-1,])
 
 # total treatment effect
 RR_tot<-St.a1m1/(St.a0m0*as.numeric(St.a0m0>0)+(1-as.numeric(St.a0m0>0)))
 DIF_tot<-log(St.a1m1)-log(St.a0m0) 

 # direct effect of Trt on risk ratio
 RR_d<-St.a1m0/(St.a0m0*as.numeric(St.a0m0>0)+(1-as.numeric(St.a0m0>0)))
 DIF_d<-log(St.a1m0)-log(St.a0m0)

 # indirect effect of Trt on risk ratio
 RR_m<-St.a1m1/(St.a1m0*as.numeric(St.a1m0>0)+(1-as.numeric(St.a1m0>0)))
 DIF_m<-log(St.a1m1)-log(St.a1m0)

 # mediation proportion
 PROP_m<-(St.a1m1-St.a1m0)/(St.a1m1-St.a0m0)

 out <- as.list(data.frame(tseq,St.a0m0,St.a1m1,St.a1m0,RR_d,RR_m,RR_tot,DIF_d,DIF_m,DIF_tot,PROP_m))

 attributes(out)[["scenario"]] <- x$scenario
 attributes(out)[["pr.ctr"]] <- pr.ctr
 attributes(out)[["pr.trt"]] <- pr.trt

 attributes(out)[["summary"]] <- c( "median.ctr"=round(tseq[which.min(abs(St.a0m0-0.5))],2),
                                    "median.trt"=round(tseq[which.min(abs(St.a1m1-0.5))],2) )

 class(out) <- c("riskratio_tru",class(out))

 return(out)
}
