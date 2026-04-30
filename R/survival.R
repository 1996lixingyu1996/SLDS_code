#' Specify Weibull MPTI model parameters for binary surrogate

#' @param beta coefficients for treatment, surrogate response, and prognostic factors (x) in each multistate model 
               ## beta12 -> treatment (direct effect), surrogate response, x for hazard of PFS
               ## beta13 -> treatment (direct effect), surrogate response, x for hazard of OS
               ## beta23 -> treatment (direct effect), surrogate response, x for hazard of OS after progression
#' @param gamma vector of intercepts and regression coefficients for treatment and baseline covariates x in surrogate response model
#' @param weib_nv baseline hazard Weibull shape parameters for PFS (12), OS (13), and OS (23) after progression
               ## nv12 -> shape for PFS
               ## nv13 -> shape for OS
               ## nv23 -> shape for OS after progression
#' @param weib_lamb baseline hazard Weibull rate parameters for PFS (12), OS (13), and OS (23) after progression
               ## lamb12 -> baseline rate for PFS
               ## lamb13 -> baseline rate for OS
               ## lamb23 -> baseline rate for OS after progression
#' @param prob binomial probability vector for binary prognostic covariates
#' @param mu mean vector for numeric prognostic covariates
#' @param sigma stand dev vector for numeric prognostic covariates
#' @param weib_shape Weibull shape paramter for time-to-observing surrogate response
#' @param weib_scale Weibull scale paramter for time-to-observing surrogate response
#' @param scenario text string to label the scenario
#' @return list of model parameters
#' @export

parm_mpti_weib <- function(beta,gamma,weib_nv,weib_lamb,prob=NULL,mu=NULL,sigma=NULL,weib_shape=NULL,weib_scale=NULL,scenario=NULL){
 
 if( length(beta)<1 ){ stop("beta required") }
 if( length(gamma)<1 ){ stop("gamma required") }
 if( length(weib_nv)<1 ){ stop("weib_nv required") }
 if( length(weib_lamb)<1 ){ stop("weib_lamb required") }

 if( length(prob)>0 ){ if( prob<0 || prob>1 ){
  stop("prob must be on unit interval (0,1)")
 }}

 if( length(sigma)>0 ){ if( sigma<0 ){
  stop("sigma must be > 0")
 }}

 if( length(weib_shape)>0 ){ if( weib_shape<0 ){
  stop("weib_shape must be > 0")
 }} 

 if( length(weib_scale)>0 ){ if( weib_scale<0 ){
  stop("weib_scale must be > 0")
 }} 

 out <- list(beta=beta,gamma=gamma,weib_nv=weib_nv,weib_lamb=weib_lamb,prob=prob,mu=mu,sigma=sigma,weib_shape=weib_shape,weib_scale=weib_scale,scenario=scenario)

 class(out) <- c("mpti_weib_parm",class(out))

 return(out)
}


#' Survival function for Weibull MPTI model over interval 0 to t (probability of state 1 or 2)

#' @param parm MPTI Weibull model parameters
#' @param t time interval upper bound (0, t)
#' @param z vector of model matrix row (treatment, mediatior, prognostic covariates)
#' @param npart length out time grid
#' @return vector of survival probabilities over time 0, t
#' @export

surv_mpti_weib <- function(parm,t,z,npart=300){
 
 tseq <- seq(0,t,length.out=npart)

 if( !inherits(parm,"mpti_weib_parm")  ){
  stop("parm must be of class mpti_weib_parm")
 }
 
 out <- sapply(tseq, FUN=s_mpti_weib, z=z, weib_nv=parm$weib_nv, weib_lamb=parm$weib_lamb, beta=parm$beta)

 return(out)
}


# Survival function for Weibull MPTI model at time t (probability of state 1 or 2)

s_mpti_weib <- function(t,z,weib_nv,weib_lamb,beta){

 p_11 <- function(u,t,z,weib_nv,weib_lamb,beta){
  Hazd12<-weib_lamb$lamb12*((t^weib_nv$nv12)-(u^weib_nv$nv12))*exp(z%*%beta$beta12)
  Hazd13<-weib_lamb$lamb13*((t^weib_nv$nv13)-(u^weib_nv$nv13))*exp(z%*%beta$beta13)
  return(exp(-Hazd12-Hazd13))
 }

 dp_12 <- function(u,t,z,weib_nv,weib_lamb,beta){
  p_22 <- exp(-weib_lamb$lamb23*((t^weib_nv$nv23)-(u^weib_nv$nv23))*exp(z%*%beta$beta23))
  lam_12 <- weib_nv$nv12*weib_lamb$lamb12*(u^(weib_nv$nv12-1))*exp(z%*%beta$beta12)
  p_11 <- p_11(u=0,t=u,z=z,weib_nv=weib_nv,weib_lamb=weib_lamb,beta=beta)
  return(p_22*lam_12*p_11)
 }

 rsums <- function(ub,z,weib_nv,weib_lamb,beta,l=1000){
  g <- seq(0,ub,length.out=l)
  delt <- (g[2]-g[1])
  grid <- g[-length(g)]+(delt/2)
  return( sum( sapply(grid,FUN=dp_12,t=ub,z=z,weib_nv=weib_nv,weib_lamb=weib_lamb,beta=beta)*delt ) )
 }

 out <- p_11(u=0,t=t,z=z,weib_nv=weib_nv,weib_lamb=weib_lamb,beta=beta)
 out <- out + rsums(ub=t,z=z,weib_nv=weib_nv,weib_lamb=weib_lamb,beta=beta)

 return(out)
}

