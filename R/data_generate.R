#' Randomly generate treatment assignment for n pts
#'
#' @param n total sample size
#' @param p allocation probability to treatment
#' @return vector of assignments

gen_a <- function(n,p=NULL,seed=1){
  set.seed(seed)
 if(length(p)<1){ p <- 0.5 }
 out <- rep(0,n)
 out[ sample(n, ceiling(n*p), replace = FALSE) ] <- 1
 return(out)
}

#' Randomly generate independent patient baseline/prognostic covariate from binomial and normal distributions
#'
#' @param n total sample size
#' @param pb number of binary covariates
#' @param pc number of continuous covariates
#' @param prob vector of probabilities for binary covariates, default 0.5
#' @param mu vector of mean for continuous covariates, default 0
#' @param sigma vector of sigma for continuous covariates, default 1
#' @importFrom stats rbinom rnorm
#' @return model matrix

 gen_x <- function(n, pb=1, pc=1, prob=NULL, mu=NULL, sigma=NULL,seed=1){
  set.seed(seed)
  p <- pb+pc # number of covariates
  if(is.null(prob)) prob<-rep(0.5,pb)
  if(is.null(mu)) mu<-rep(0,pc)
  if(is.null(sigma)) sigma<-rep(1,pc)

  X.mat<-matrix(nrow=n,ncol=p)
  if(pb>0){
   for(j in 1:pb){
    X.mat[,j]<-stats::rbinom(n,1,prob[j])
   }
  }
  if(pc>0){
   for(k in 1:pc){
    X.mat[,pb+k]<-stats::rnorm(n,mu[k],sigma[k])
   }
  }

  return(X.mat)
 }

#' Randomly generate independent patient mediators from binomial and normal distributions
#'
#' @param n total sample size
#' @param qb number of binary mediators
#' @param qc number of continuous normal mediator
#' @param a vector of treatment arm indicators
#' @param x matrix of baseline covariates
#' @param int vector of qb+qc intercepts in linear predictor for mediator models, default 1
#' @param coef_a vector of qb+qc coefficients for a in linear predictor of surrogate models, default 0
#' @param coef_x matrix of coefficients for x in linear predictor of surrogate models, one row per mediator, default 0
#' @param err_sigm sd of error terms in normal mediators, length of qc, default 0.5
#' @importFrom stats rbinom rnorm
#' @return model matrix

gen_m <- function(n, qb=1, qc=1, a=NULL, x=NULL, int=NULL, coef_a=NULL, coef_x=NULL, err_sigm=NULL, seed=1){
  set.seed(seed)
  q_n <- qb+qc # number of mediators
  if(is.null(int)) int<-rep(1,q_n)
  if(is.null(coef_a)) coef_a<-rep(0,q_n)
  if(is.null(coef_x)) coef_x<-matrix(0,nrow=q_n,ncol=ncol(x))
  if(is.null(err_sigm)) err_sigm<-rep(0.5,qc)
  modmat<-cbind(1,a,x)

  M_mat<-matrix(nrow=n,ncol=q_n)
  if(qb>0){
   for(j in 1:qb){
    coef<-c(int[j],coef_a[j],coef_x[j,])
    lp_tmp<-modmat%*%coef
    M_mat[,j]<-stats::rbinom(n,1,exp(lp_tmp)/(1+exp(lp_tmp)))
   }
  }
  if(qc>0){
   for(k in 1:qc){
    coef<-c(int[qb+k],coef_a[qb+k],coef_x[qb+k,])
    lp_tmp<-modmat%*%coef
    M_mat[,qb+k]<-lp_tmp+stats::rnorm(n,0,err_sigm[k])
   }
  }

 return(M_mat)
}


#' Randomly generate observation durations for mediators
#'
#' @param n total sample size
#' @param q number of mediator that will not be observed at baseline
#' @param eos_dur end of study duration, used if the mediation time is larger than all visit time
#' @param weib_shape vector of q weibull shape pars of incidence time for mediation variables
#' @param weib_scale vector of q weibull scale pars of incidence time for mediation variables
#' @param by_visit whether to use grid of study visits for the incidence time
#' @param mvis mean number of visits for each patient
#' @param fix_gapt whether fix the gap time at mgapt, if not will be generated from exp dist with mean mgapt
#' @param mgapt mean gap time between two visits
#' @importFrom stats rbinom rnorm rpois rweibull rexp
#' @return list of generated and observed observation durations

gen_m_obs <- function(n, q=1, eos_dur=50, weib_shape=NULL, weib_scale=NULL, by_visit=TRUE, mvis=NULL, fix_gapt=TRUE, mgapt=NULL,seed=1){
  set.seed(seed)
  if(is.null(weib_shape)) weib_shape<-rep(1,q)
  if(is.null(weib_scale)) weib_scale<-rep(3,q)
  if(is.null(mvis)) mvis<-5
  if(is.null(mgapt)) mgapt<-1.5

  Mt_mat<-Mt_obs<-matrix(nrow=n,ncol=q)

  for(k in 1:q){
    Mt_mat[,k]<-stats::rweibull(n,shape=weib_shape[k],scale=weib_scale[k])
  }
  if(by_visit){
    nvis<-1+stats::rpois(n,mvis)
    for(i in 1:n){
      if(!fix_gapt) gap_t<-c(0,stats::rexp(nvis[i],1/mgapt),eos_dur)
      if(fix_gapt) gap_t<-c(0,rep(mgapt,nvis[i]))
      vist<-c(cumsum(gap_t),eos_dur)
      for(j in 1:q) Mt_obs[i,j]<-min(vist[Mt_mat[i,j]<=vist])
    }
  } else Mt_obs<-Mt_mat

 return(list(trueMt=Mt_mat,obsMt=Mt_obs))
}


#' Randomly generate survival endpoint duration from Weibull PH model
#'
#' @param n total sample size
#' @param a treatment arm
#' @param x baseline covariates
#' @param m mediators
#' @param coef_a coef of treatment a
#' @param coef_x coefs of x
#' @param coef_m coefs of m
#' @param weibpars weib baseline distribution
#' @importFrom stats runif
#' @return vector of endpoint durations

gen_dur_weibph <- function(n, a=matrix(nrow=n,ncol=0), x=matrix(nrow=n,ncol=0), m=matrix(nrow=n,ncol=0),
                           coef_a=0, coef_x=rep(0,ncol(x)), coef_m=rep(0,ncol(m)), weibpars=list(nv=1,lamb=0.05)){

  mat_s<-cbind(a,m,x)
  lp_s<-mat_s%*%c(coef_a,coef_m,coef_x)
  U<-stats::runif(n)
  T<-(-log(U)/exp(lp_s)/weibpars$lamb)^(1/weibpars$nv)
  return(as.numeric(T))
}


# Multi-state model optimization function (weibull - MPTI model)
# cumulative intensity functions for transitions 12 and 13

t_search <- function(t,u,z,weib_nv,weib_lamb,beta){
  Hazd12<-weib_lamb$lamb12*t^(weib_nv$nv12)*exp(z%*%beta$beta12)
  Hazd13<-weib_lamb$lamb13*t^(weib_nv$nv13)*exp(z%*%beta$beta13)
  return((exp(-Hazd12-Hazd13)-u)^2)
}


#' Generate peicewise accural
#'
#' @param n total sample size
#' @param accTime time parition for simulating accrual
#' @param accExN expected enrollment for each bin of accTime
#' @importFrom stats runif
#' @return sorted enrollment times for n patients
#' @export

gen_accural <- function(n,accTime=40,accExN=NULL){

 out <- rep(NA,n)

 if(length(accTime)<2){

  out <- stats::runif(n,0,accTime)

 } else{

  dE <- c(accExN[[1]],base::diff(accExN))
  grp <- factor(sort(sample(1:length(accExN), n, replace=TRUE, prob=dE/sum(dE))))

  w <- which(grp==levels(grp)[[1]])
  out[w] <- stats::runif(length(w),0,accTime[[1]])

  for(ii in 2:length(levels(grp))){
   w <- which(grp==levels(grp)[[ii]])
   out[w] <- stats::runif(length(w),accTime[[ii-1]],accTime[[ii]])
  }

 }

 return(sort(out))
}


#' Generate dataset from Weibull multi-state mediation model with binary surrogate (mediator)
#'
#' @param parm parameter list of class mpti_weib_parm
#' @param n total sample size
#' @param nb number of binary prognostic covariates
#' @param nc number of continous (numeric) prognostic covariates
#' @param prtrt prob of treatment
#' @param accTime time parition for simulating accrual
#' @param accExN expected enrollment for each bin of accTime
#' @param addFUt added follow-up duration in months before final analysis
#' @param CenUpLim upper bound of right censoring for final analysis (longest duration possible to observe)
#' @param tunits label for time units
#' @param eos_stats compute end of study stats
#' @importFrom stats runif optimize rbinom
#' @return simulated dataset data frame
#' @export

gen_mstate <- function(parm,n,nb=1,nc=1,prtrt=0.5,accTime=40,accExN=NULL,addFUt=0,CenUpLim=NULL,tunits="months",eos_stats=TRUE){

  if( !inherits(parm,"mpti_weib_parm")  ){
   stop("parm must be of class mpti_weib_parm")
  }

  if( length(parm$gamma) != (2+nb+nc) ){ stop("length of parm$gamma must equal 2 + nb + nc") }

  eos_dur <- max(accTime)+addFUt

  accrt <- gen_accural(n=n,accTime=accTime,accExN=accExN)

  if(length(CenUpLim)<1){ CenUpLim <- eos_dur }

  a <- gen_a(n=n,p=prtrt)

  x <- NULL

  if( (nb+nc)>0 ){
   x <- gen_x(n=n, pb=nb, pc=nc, prob=parm$prob, mu=parm$mu, sigma=parm$sigma)
  }

  coef_x <- NULL

  if( length(x)>0 ){
   coef_x <- matrix(parm$gamma[3:(2+nb+nc)],1,nb+nc)
  }

  resp <- gen_m(n=n,qb=nb,qc=0,a=a,x=x,int=parm$gamma[1],coef_a=parm$gamma[2],coef_x=coef_x)

  if( length(parm$weib_shape)>0 && length(parm$weib_scale)>0 ){
   resp_t <- gen_m_obs(n=n,q=1,weib_shape=parm$weib_shape,weib_scale=parm$weib_scale)$trueMt
  } else{
   resp_t <- gen_m_obs(n=n,q=1,eos_dur=eos_dur)$trueMt
  }

  resp_t[resp==0]<-NA

  modmat<-cbind(a,resp,x)

  U<-stats::runif(n)
  if(parm$weib_nv$nv12==parm$weib_nv$nv13){
    T<-(-log(U)/(parm$weib_lamb$lamb12*exp(modmat%*%parm$beta$beta12)+parm$weib_lamb$lamb13*exp(modmat%*%parm$beta$beta13)))^(1/parm$weib_nv$nv12)
  }else{
    T<-Obj<-numeric()
    for(i in 1:n){
      dif<-1; d<-0
      while(dif>0.001 & d<100){
       tfit0<-stats::optimize(t_search,interval=c(0,1+d),u=U[i],z=modmat[i,],weib_nv=parm$weib_nv,weib_lamb=parm$weib_lamb,beta=parm$beta)
       T[i]<-tfit0$minimum
       Obj[i]<-dif<-tfit0$objective
       d<-d+1
      }
    }
  }
  hazd12<-parm$weib_nv$nv12*parm$weib_lamb$lamb12*T^(parm$weib_nv$nv12-1)*exp(modmat%*%parm$beta$beta12)
  hazd13<-parm$weib_nv$nv13*parm$weib_lamb$lamb13*T^(parm$weib_nv$nv13-1)*exp(modmat%*%parm$beta$beta13)
  pd1<-hazd12/(hazd12+hazd13)
  delta1<-stats::rbinom(n,1,pd1)
  T1<-T2<-T
  T2.1<-(T1^(parm$weib_nv$nv23)-log(stats::runif(n))/parm$weib_lamb$lamb23/exp(modmat%*%parm$beta$beta23))^(1/parm$weib_nv$nv23)
  T2[delta1==1]<-T2.1[delta1==1]

  C<-apply(cbind(accrt+stats::runif(n,0,CenUpLim),eos_dur),1,min)
  delta2<-as.numeric(T2<=C)
  delta1<-as.numeric(T1<=C)*delta1
  T1.obs<-apply(cbind(T1,C),1,min)
  T2.obs<-apply(cbind(T2,C),1,min)

  mydat<-data.frame(cbind(accrt,resp_t,T1.obs,T2.obs,delta1,delta2,a,resp,resp_t,resp,x))
  names(mydat)<-c("accT","resp_t","PFS","OS","PFSevent","OSevent","trt","resp","latent_resp_t","latent_resp","x1","x2")

  dthf <- which( mydat$OSevent==1 )
  if( length(dthf)>0 ){
   mydat$PFSevent[dthf] <- 1
  }

  ww <- which( mydat$resp==1 )

  if( length(ww)>0 ){ for(ii in ww){

   if( mydat$resp_t[[ii]] > mydat$PFS[[ii]] ){
    mydat$resp_t[[ii]] <- NA
    mydat$resp[[ii]] <- 0
   }

  }}

  mydat$DOR <- NA
  mydat$DORevent <- NA
  uresp <- which(mydat$resp==1)
  if( length(uresp)>0 ){
   mydat[uresp,"DOR"] <- mydat[uresp,"PFS"]
   mydat[uresp,"DORevent"] <- mydat[uresp,"PFSevent"]
  }

  mydat$TTP <- mydat$PFS
  mydat$TTPevent <- mydat$PFSevent
  udeath <- which(mydat$OSevent==1)
  if( length(udeath)>0 ){ for(ii in udeath){
   if( mydat$PFS[[ii]]==mydat$OS[[ii]] ){
    mydat$TTPevent[[ii]] <- 0
   }
  }}

  mydat[,"USUBJID"] <- as.character(1:nrow(mydat))
  mydat[,"enrolled"] <- TRUE

  attributes(mydat)[["type"]] <- "bimed"
  attributes(mydat)[["scenario"]] <- parm$scenario
  attributes(mydat)[["eos"]] <- eos_dur
  attributes(mydat)[["tunits"]] <- tunits

  class(mydat) <- c("mpti_weib_dat",class(mydat))

  #### End of study analysis ####
  if( eos_stats ){
   attributes(mydat)[["eos_stats"]] <- eos_analysis(dat=mydat)
  }

  #### Add multi-state pathway info ####
  out <- add_pathlos_data(dat=mydat)

  attributes(out)[["args"]] <- list(n=n,nb=nb,nc=nc,prtrt=prtrt,accTime=accTime,accExN=accExN,addFUt=addFUt,CenUpLim=CenUpLim,tunits=tunits)

  dr <- NULL; fup <- NULL; dths <- NULL
  for(ii in as.character(unique(out$trt))){
   temp <- base::subset(out, out$trt==ii & out$enrolled)
   if( nrow(temp)>0 ){
    dths[[ii]] <- length(grep(",DEATH",temp$path,fixed=TRUE))
    dr[[ii]] <- length(grep(",DEATH",temp$path,fixed=TRUE))/sum(temp$OS)
    fup[[ii]] <- sum(temp$OS)
   } else{
    dths[[ii]] <- 0; dr[[ii]] <- 0; fup[[ii]] <- 0
   }
   temp <- NULL
  }
  attributes(out)[["deaths"]] <- dths
  attributes(out)[["death_rate"]] <- dr
  attributes(out)[["os_fup"]] <- fup

 return(out)
}


#' Cut multi-state data at interim analysis time tp
#'
#' @param dat either data.frame of class mpti_weib_dat or object of class mstate_data
#' @param ... additional arguments to class specific function calls
#' @export

interim_mstate <- function(dat,...){
 UseMethod("interim_mstate")
}

#' Cut multi-state data at interim analysis time tp mstate_data
#'
#' @param dat data.frame of class mpti_weib_dat
#' @param tp interim analysis trial time
#' @param delt time lag of minimal follow-up required to be "enrolled"
#' @return data.frame of class impti_weib_dat
#' @export

interim_mstate.mstate_data <- function(dat,tp,delt=0){

 fdat <- dat@data
 idat <- interim_mstate(dat=fdat, tp=tp, delt=delt)

 new("imstate_data",
      data =  idat,
      cntrl = dat@cntrl,
      exper = dat@exper,
      resp_levels = dat@resp_levels,
      resp_t = dat@resp_t,
      covar = dat@covar,
      ttp = dat@ttp,
      dor = dat@dor)
}

#' Cut multi-state data at interim analysis time tp mpti_weib_dat
#'
#' @param dat data.frame of class mpti_weib_dat
#' @param tp interim analysis trial time
#' @param delt time lag of minimal follow-up required to be "enrolled"
#' @return data.frame of class impti_weib_dat
#' @export

interim_mstate.mpti_weib_dat <- function(dat,tp,delt=0){

  dat$path <- NULL
  dat$SD <- NULL
  dat$RS <- NULL
  dat$PD <- NULL

  dat$enrolled <- dat$accT <= (tp-delt)

  if( sum(dat$enrolled)<1  ){

   out_1 <- base::subset(dat, rep(FALSE, nrow(dat)))

  } else{

   out_1 <- base::subset(dat,dat$enrolled)

  }

  if( sum(dat$enrolled)==nrow(dat)  ){

   out_2 <- base::subset(dat, rep(FALSE, nrow(dat)))

  } else{

   out_2 <- base::subset(dat,!dat$enrolled)

  }

  if( nrow(out_1)>0 ){

   resp_t1_0  <- out_1$accT + out_1$resp_t
   PFS1_0     <- out_1$accT + out_1$PFS
   OS1_0      <- out_1$accT + out_1$OS

   resp1 <- rep(0,nrow(out_1))
   DOR1 <- DORevent1 <- resp_t1 <- rep(NA,nrow(out_1))

   PFS1 <- apply(cbind(PFS1_0,tp),1,min)-out_1$accT
   PFSevent1 <- as.numeric(PFS1_0<=tp & out_1$PFSevent==1)

   OS1 <- apply(cbind(OS1_0,tp),1,min)-out_1$accT
   OSevent1 <- as.numeric(OS1_0<=tp & out_1$OSevent==1)

   w1 <- which( resp_t1_0<=tp & out_1$resp==1 )
   w0 <- base::setdiff( 1:nrow(out_1), w1 )

   if( length(w1)>0 ){
    resp1[w1] <- 1
    resp_t1[w1] <- resp_t1_0[w1] - out_1$accT[w1]
    DOR1[w1] <- PFS1[w1]
    DORevent1[w1] <- PFSevent1[w1]
   }

   out_1$resp <- resp1
   out_1$resp_t <- resp_t1
   out_1$PFS <- PFS1
   out_1$PFSevent <- PFSevent1
   out_1$OS <- OS1
   out_1$OSevent <- OSevent1
   out_1$DOR <- DOR1
   out_1$DORevent <- DORevent1

   n <- nrow(out_1)
   n_resp <- sum(out_1$resp)
   n_PFSevent <- sum(out_1$PFSevent)
   n_OSevent <- sum(out_1$OSevent)

  } else{

   n <- 0
   n_resp <- 0
   n_PFSevent <- 0
   n_OSevent <- 0

  }

  if( nrow(out_2)>0 ){

   out_2$resp <- NA
   out_2$resp_t <- NA
   out_2$PFS <- NA
   out_2$PFSevent <- NA
   out_2$OS <- NA
   out_2$OSevent <- NA
   out_2$DOR <- NA
   out_2$DORevent <- NA

  }

  out <- rbind(out_1,out_2)

  out$TTP <- out$PFS
  out$TTPevent <- out$PFSevent
  udeath <- which(out$OSevent==1)
  if( length(udeath)>0 ){ for(ii in udeath){
   if( !is.na(out$PFS[[ii]]) && !is.na(out$OS[[ii]]) ){
    if( out$PFS[[ii]]==out$OS[[ii]] ){
     out$TTPevent[[ii]] <- 0
    }
   }
  }}

  attributes(out)[["type"]] <- attributes(dat)[["type"]]
  attributes(out)[["scenario"]] <- attributes(dat)[["scenario"]]
  attributes(out)[["eos"]] <- attributes(dat)[["eos"]]
  attributes(out)[["tunits"]] <- attributes(dat)[["tunits"]]
  attributes(out)[["args"]] <- attributes(dat)[["args"]]

  class(out) <- c("impti_weib_dat",class(out))

  attributes(out)[["summary"]] <- c(trial_time=tp, delt=delt, n=n, n.resp=n_resp, n.PFSevent=n_PFSevent, n.OSevent=n_OSevent)

  if( length(attributes(dat)[["eos_stats"]])>0 ){
   attributes(out)[["eos_stats"]] <- attributes(dat)[["eos_stats"]]
  }

  #### Add multi-state pathway info ####
  if( any(out$enrolled) ){
   out <- add_pathlos_data(dat=out)
  }

  dr <- NULL; fup <- NULL; dths <- NULL
  for(ii in as.character(unique(out$trt))){
   temp <- base::subset(out, out$trt==ii & out$enrolled)
   if( nrow(temp)>0 ){
    dths[[ii]] <- length(grep(",DEATH",temp$path,fixed=TRUE))
    dr[[ii]] <- length(grep(",DEATH",temp$path,fixed=TRUE))/sum(temp$OS)
    fup[[ii]] <- sum(temp$OS)
   } else{
    dths[[ii]] <- 0; dr[[ii]] <- 0; fup[[ii]] <- 0
   }
   temp <- NULL
  }
  attributes(out)[["deaths"]] <- dths
  attributes(out)[["death_rate"]] <- dr
  attributes(out)[["os_fup"]] <- fup

  return(out)
}


#' Add multi-state path and length of stay data at interim analysis time tp
#'
#' @param dat data.frame of class mpti_weib_dat or impti_weib_dat
#' @param tp interim analysis trial time
#' @param delt time lag of minimal follow-up required to be "enrolled"
#' @importFrom stats na.omit
#' @return data.frame of class impti_weib_dat
#' @export

add_pathlos_data <- function(dat,tp=NULL,delt=NULL){

 #### org length of stay data ####
 comp_los <- function(dat,lab){

  los <- list( "SD"= as.numeric(stats::na.omit(dat$SD)),
               "RS" = as.numeric(stats::na.omit(dat$RS)),
               "PD" = as.numeric(stats::na.omit(dat$PD)) )

  ldat <- as.data.frame(list( x=NA, y=NA, col=NA ))
  ldat <- base::subset(ldat, FALSE)
  for(ii in names(los)){
   temp <- as.data.frame(list(x=los[[ii]], y=rep(ii,length(los[[ii]]))))
   if( nrow(temp)>0 ){
    if(ii=="SD"){
     temp$col <- "#808080"
    } else if(ii=="PD"){
     temp$col <- "#3D4D8AFF"
    } else if(ii=="RS"){
     temp$col <- "#FDE725FF"
    }
    ldat <- rbind(ldat,temp)
   }
  }

  ldat$trt <- lab

  return(ldat)
 }

 ########################################################################

 if( !inherits(dat,"mpti_weib_dat") && !inherits(dat,"impti_weib_dat") ){
  stop("dat must be of class mpti_weib_dat or impti_weib_dat")
 }

 if( length(tp)<1 ){

  if( inherits(dat,"impti_weib_dat") ){
   tp <- attributes(dat)[["summary"]][["trial_time"]]
  } else{
   tp <- attributes(dat)[["eos"]]
  }

 }

 if( length(delt)<1 ){

  if( inherits(dat,"impti_weib_dat") ){
   delt <- attributes(dat)[["summary"]][["delt"]]
  } else{
   delt <- 0
  }

 }

 if( !inherits(dat,"impti_weib_dat") && length(tp)>0 ){
  if( tp < attributes(dat)[["eos"]] ){
   dat <- interim_mstate(dat=dat, tp=tp, delt=delt)
  }
 }

 final <- TRUE

 if( tp < attributes(dat)[["eos"]] ){
  final <- FALSE
 }

 eval_i <- 1:nrow(dat)

 if( inherits(dat,"impti_weib_dat") ){
  eval_i <- which(dat$enrolled)
 }

 ### add multi-state transition data ###
 dat$path <- NA
 dat$SD <- NA
 dat$RS <- NA
 dat$PD <- NA

 if( length(eval_i)<1 ){
  #stop("dat has no enrollment")
  ldat <- base::subset(as.data.frame(list( x=NA, y=NA, col=NA, trt=NA )),FALSE)
  attributes(dat)[["los_data"]] <- ldat
  return( dat )
 }

 ### iterate for each patient ###

 for(ii in eval_i){

  t1 <- "SD"
  t2 <- "."
  t3 <- "."
  t4 <- "OBS"

  dat$SD[[ii]] <- dat$OS[[ii]]

  if( !dat$enrolled[[ii]] ){ ### not yet enrolled
   dat$path[[ii]] <- NA
   next
  }

  if( dat$resp[[ii]]==1 ){

   dat$SD[[ii]] <- dat$resp_t[[ii]]
   t2 <- "RS"

  } else if( dat$PFSevent[[ii]]==1 && dat$OS[[ii]]!=dat$PFS[[ii]] ){

   dat$SD[[ii]] <- dat$PFS[[ii]]
   dat$PD[[ii]] <- dat$OS[[ii]] - dat$PFS[[ii]]
   t2 <- "PD"

  }

  if( t2=="RS" ){

   dat$RS[[ii]] <- dat$OS[[ii]] - dat$resp_t[[ii]]

   if( dat$PFSevent[[ii]]==1 && dat$OS[[ii]]!=dat$PFS[[ii]] ){

    dat$RS[[ii]] <- dat$PFS[[ii]] - dat$resp_t[[ii]]
    dat$PD[[ii]] <- dat$OS[[ii]] - dat$PFS[[ii]]
    t3 <- "PD"

   }

  }

  if( dat$OSevent[[ii]]==1 ){

   t4 <- "DEATH"

  } else if(final || (dat$OS[[ii]] < (tp-dat$accT[[ii]])) ){

   t4 <- "EOS"

  }

  dat$path[[ii]] <- paste(c(t1,t2,t3,t4),collapse=",")

 }

 ldat <- as.data.frame(list( x=NA, y=NA, col=NA, trt=NA ))
 for(jj in unique(dat$trt)){
  temp <- base::subset(dat, dat$trt==jj)
  if( any(temp$enrolled) ){
   ldat_0 <- comp_los(dat=temp, lab=jj)
   ldat <- rbind(ldat, ldat_0)
  }
 }

 ldat <- base::subset(ldat, !is.na(ldat$trt))

 attributes(dat)[["los_data"]] <- ldat

 return( dat )
}



#' Randomly generate treatment assignment for n pts
#'
#' @param n total sample size
#' @param a treatment indicator
#' @param x covariate variables
#' @param nb number of binary variables in weibull model
#' @param nc number of continuous variables in weibull model
#' @return vector of assignments
#'
gen_heterogeneous_data <- function(parm, n, a, x, nb,nc,prtrt=0.5,accTime=40,accExN=NULL,addFUt=0,CenUpLim=NULL,tunits="months",eos_stats=TRUE,seed=1){
  set.seed(seed)
  eos_dur <- max(accTime)+addFUt

  accrt <- gen_accural(n=n,accTime=accTime,accExN=accExN)

  if(length(CenUpLim)<1){ CenUpLim <- eos_dur }

  #a <- gen_a(n=n,p=prtrt)

  #x <- NULL

  #if( (nb+nc)>0 ){
  #  x <- gen_x(n=n, pb=nb, pc=nc, prob=parm$prob, mu=parm$mu, sigma=parm$sigma)
  #}

  coef_x <- NULL

  if( ncol(x)>0 ){
    coef_x <- matrix(parm$gamma[3:(2+nb+nc)],1,nb+nc)
  }

  resp <- gen_m(n=n,qb=1,qc=0,a=a,x=x,int=parm$gamma[1],coef_a=parm$gamma[2],coef_x=coef_x)

  if( length(parm$weib_shape)>0 && length(parm$weib_scale)>0 ){
    resp_t <- gen_m_obs(n=n,q=1,weib_shape=parm$weib_shape,weib_scale=parm$weib_scale)$trueMt
  } else{
    resp_t <- gen_m_obs(n=n,q=1,eos_dur=eos_dur)$trueMt
  }

  resp_t[resp==0]<-NA

  modmat<-cbind(a,resp,x)

  U<-stats::runif(n)
  if(parm$weib_nv$nv12==parm$weib_nv$nv13){
    T<-(-log(U)/(parm$weib_lamb$lamb12*exp(modmat%*%parm$beta$beta12)+parm$weib_lamb$lamb13*exp(modmat%*%parm$beta$beta13)))^(1/parm$weib_nv$nv12)
  }else{
    T<-Obj<-numeric()
    for(i in 1:n){
      dif<-1; d<-0
      while(dif>0.001 & d<100){
        tfit0<-stats::optimize(t_search,interval=c(0,1+d),u=U[i],z=modmat[i,],weib_nv=parm$weib_nv,weib_lamb=parm$weib_lamb,beta=parm$beta)
        T[i]<-tfit0$minimum
        Obj[i]<-dif<-tfit0$objective
        d<-d+1
      }
    }
  }
  hazd12<-parm$weib_nv$nv12*parm$weib_lamb$lamb12*T^(parm$weib_nv$nv12-1)*exp(modmat%*%parm$beta$beta12)
  hazd13<-parm$weib_nv$nv13*parm$weib_lamb$lamb13*T^(parm$weib_nv$nv13-1)*exp(modmat%*%parm$beta$beta13)
  pd1<-hazd12/(hazd12+hazd13)
  delta1<-stats::rbinom(n,1,pd1)
  T1<-T2<-T
  T2.1<-(T1^(parm$weib_nv$nv23)-log(stats::runif(n))/parm$weib_lamb$lamb23/exp(modmat%*%parm$beta$beta23))^(1/parm$weib_nv$nv23)
  T2[delta1==1]<-T2.1[delta1==1]

  C<-apply(cbind(accrt+stats::runif(n,0,CenUpLim),eos_dur),1,min)
  delta2<-as.numeric(T2<=C)
  delta1<-as.numeric(T1<=C)*delta1
  T1.obs<-apply(cbind(T1,C),1,min)
  T2.obs<-apply(cbind(T2,C),1,min)

  mydat<-data.frame(cbind(accrt,resp_t,T1.obs,T2.obs,delta1,delta2,a,resp,resp_t,resp,x))
  names(mydat)<-c("accT","resp_t","PFS","OS","PFSevent","OSevent","trt","resp","latent_resp_t","latent_resp",colnames(x))

  dthf <- which( mydat$OSevent==1 )
  if( length(dthf)>0 ){
    mydat$PFSevent[dthf] <- 1
  }

  ww <- which( mydat$resp==1 )

  if( length(ww)>0 ){ for(ii in ww){

    if( mydat$resp_t[[ii]] > mydat$PFS[[ii]] ){
      mydat$resp_t[[ii]] <- NA
      mydat$resp[[ii]] <- 0
    }

  }}

  mydat$DOR <- NA
  mydat$DORevent <- NA
  uresp <- which(mydat$resp==1)
  if( length(uresp)>0 ){
    mydat[uresp,"DOR"] <- mydat[uresp,"PFS"]
    mydat[uresp,"DORevent"] <- mydat[uresp,"PFSevent"]
  }

  mydat$TTP <- mydat$PFS
  mydat$TTPevent <- mydat$PFSevent
  udeath <- which(mydat$OSevent==1)
  if( length(udeath)>0 ){ for(ii in udeath){
    if( mydat$PFS[[ii]]==mydat$OS[[ii]] ){
      mydat$TTPevent[[ii]] <- 0
    }
  }}

  mydat[,"USUBJID"] <- as.character(1:nrow(mydat))
  mydat[,"enrolled"] <- TRUE

  attributes(mydat)[["type"]] <- "bimed"
  attributes(mydat)[["scenario"]] <- parm$scenario
  attributes(mydat)[["eos"]] <- eos_dur
  attributes(mydat)[["tunits"]] <- tunits

  class(mydat) <- c("mpti_weib_dat",class(mydat))

  #### End of study analysis ####
  if( eos_stats ){
    attributes(mydat)[["eos_stats"]] <- eos_analysis(dat=mydat)
  }

  #### Add multi-state pathway info ####
  out <- add_pathlos_data(dat=mydat)

  attributes(out)[["args"]] <- list(n=n,nb=nb,nc=nc,prtrt=prtrt,accTime=accTime,accExN=accExN,addFUt=addFUt,CenUpLim=CenUpLim,tunits=tunits)

  dr <- NULL; fup <- NULL; dths <- NULL
  for(ii in as.character(unique(out$trt))){
    temp <- base::subset(out, out$trt==ii & out$enrolled)
    if( nrow(temp)>0 ){
      dths[[ii]] <- length(grep(",DEATH",temp$path,fixed=TRUE))
      dr[[ii]] <- length(grep(",DEATH",temp$path,fixed=TRUE))/sum(temp$OS)
      fup[[ii]] <- sum(temp$OS)
    } else{
      dths[[ii]] <- 0; dr[[ii]] <- 0; fup[[ii]] <- 0
    }
    temp <- NULL
  }
  attributes(out)[["deaths"]] <- dths
  attributes(out)[["death_rate"]] <- dr
  attributes(out)[["os_fup"]] <- fup

  return(out)
}
