  
  ilogit <- function(x){ 
   if(x==Inf){ return(1) }
   if(x==-Inf){ return(0) }
   return( exp(x)/(1+exp(x)) ) 
  }
  
  stat_summaries <- function(t,accT,b,trt,alternative){
   time <- apply(cbind(accT+t,b),1,min)-accT
   status <- as.numeric(accT+t<=b)
   out <- c(p=logrank_z(time,status,trt,alternative), hr=cox_hr(time,status,trt))
   out2 <- medians(time,status,trt)
   out <- c(out,out2)
   return(out)
  }

#' Organize posterior predictive distr as data.frame

#' @param x object of class pp_eos
#' @param tp trial time point
#' @return data.frame
#' @export

  pp_to_df <- function(x,tp){

   if( !inherits(x,"pp_eos") ){
    stop("x must be of class pp_eos")
   }

   out <- as.data.frame( list(type=c(rep("p",ncol(x)),rep("hr",ncol(x)),rep("med",2*ncol(x))),
                             x=rep(rep(tp,ncol(x)),4),
                             y=c(x[1,],x[2,],x[3,],x[4,]),
                             arm=c(rep(NA,2*ncol(x)),rep(gsub("trt=","",rownames(x)[3]),ncol(x)),
                                                      rep(gsub("trt=","",rownames(x)[4]),ncol(x)) )) )
   class(out) <- c(class(out),"pp_eos_df")                                      
   return(out)
  }

  eos_analysis <- function(dat){

   out <- list()
   out$OS <- list()
   out$OS$two.sided <- logrank_z(time=dat$OS, status=dat$OSevent, trt=dat$trt, alternative="two.sided")
   out$OS$less <- logrank_z(time=dat$OS, status=dat$OSevent, trt=dat$trt, alternative="less")
   out$OS$greater <- logrank_z(time=dat$OS, status=dat$OSevent, trt=dat$trt, alternative="greater")
   out$OS$hr <- cox_hr(time=dat$OS, status=dat$OSevent, trt=dat$trt)
   out$OS$medians <- medians(time=dat$OS, status=dat$OSevent, trt=dat$trt)

   out$PFS <- list()
   out$PFS$two.sided <- logrank_z(time=dat$PFS, status=dat$PFSevent, trt=dat$trt, alternative="two.sided")
   out$PFS$less <- logrank_z(time=dat$PFS, status=dat$PFSevent, trt=dat$trt, alternative="less")
   out$PFS$greater <- logrank_z(time=dat$PFS, status=dat$PFSevent, trt=dat$trt, alternative="greater")
   out$PFS$hr <- cox_hr(time=dat$PFS, status=dat$PFSevent, trt=dat$trt)
   out$PFS$medians <- medians(time=dat$PFS, status=dat$PFSevent, trt=dat$trt)

   lctr <- dat[,"trt"]==0
   x_ctr <- base::subset(dat,lctr)
   x_ctr <- base::subset(x_ctr,!is.na(x_ctr[,"resp"]))
   x_trt <- base::subset(dat,!lctr)
   x_trt <- base::subset(x_trt,!is.na(x_trt[,"resp"])) 

   out$resp <- list()
   out$resp$beta.ctr <- c(1,1)
   out$resp$beta.trt <- c(1,1)

   if( nrow(x_ctr)>0 ){
    out$resp$beta.ctr <- c(1+sum(x_ctr[,"resp"]==1), 1+nrow(x_ctr)-sum(x_ctr[,"resp"]==1))
   }

   if( nrow(x_trt)>0 ){
    out$resp$beta.trt <- c(1+sum(x_trt[,"resp"]==1), 1+nrow(x_trt)-sum(x_trt[,"resp"]==1))
   }

   out$DOR <- list()
   #out$DOR$two.sided <- logrank_z(time=dat$DOR, status=dat$DORevent, trt=dat$trt, alternative="two.sided")
   #out$DOR$less <- logrank_z(time=dat$DOR, status=dat$DORevent, trt=dat$trt, alternative="less")
   #out$DOR$greater <- logrank_z(time=dat$DOR, status=dat$DORevent, trt=dat$trt, alternative="greater")
   #out$DOR$hr <- cox_hr(time=dat$DOR, status=dat$DORevent, trt=dat$trt)
   out$DOR$medians <- medians(time=dat$DOR, status=dat$DORevent, trt=dat$trt)

   return(out)
  }

  #' @importFrom stats pnorm

  logrank_z <- function(time,status,trt,alternative){
   
   out <- 1
   all_t <- unique(sort(time[status==1]))

   if( length(all_t)>1 ){
    weight_t <- rep(1,length(all_t))
    mat_allt <- matrix(all_t,nrow=length(all_t),ncol=length(time))
    mat_time <- matrix(time,nrow=length(all_t),ncol=length(time),byrow=TRUE)
    mat_status <- matrix(status,nrow=length(all_t),ncol=length(time),byrow=TRUE)

    grps <- levels(factor(trt))
    K <- length(grps)
    d <- Y <- matrix(nrow=length(all_t),ncol=K)
  
    dmat <- mat_time==mat_allt & mat_status==1 
    Ymat <- mat_time>=mat_allt
    for(k in 1:K){
     d[,k] <- rowSums(dmat[,trt==grps[k]])
     Y[,k] <- rowSums(Ymat[,trt==grps[k]])
    }
  
    numr <- sum(weight_t*(d[,1]-Y[,1]*rowSums(d)/rowSums(Y)))
    denom <- sqrt(sum(weight_t^2*Y[,1]/rowSums(Y)*(1-Y[,1]/rowSums(Y))*(rowSums(Y)-rowSums(d))/(rowSums(Y)-1)*rowSums(d),na.rm=TRUE))
    Z <- numr/denom

    if(alternative=="greater") out <- 1-stats::pnorm(Z)
    if(alternative=="less") out <- stats::pnorm(Z)
    if(alternative=="two.sided") out <- 2*stats::pnorm(-abs(Z))
   }

   return(out)
  }

  #' @importFrom survival coxph Surv

  cox_hr <- function(time,status,trt){
   out <- NA
   fit <- NULL
   suppressWarnings( try( fit <- survival::coxph( survival::Surv(time,status) ~ trt ), silent=TRUE ) )
   if( length(fit)>0 ){ 
    out <- summary(fit)$coefficients[[2]]
   }
   return(out)
  }

  #' @importFrom survival survfit Surv

  medians <- function(time,status,trt){
   out <- c("trt=0"=NA,"trt=1"=NA)
   fit <- NULL
   suppressWarnings( try( fit <- survival::survfit( survival::Surv(time,status) ~ trt ), silent=TRUE ) )
   if( length(fit)>0 ){ if( length(nrow(summary(fit)$table))==1 ){
    out <- summary(fit)$table[,"median"]
   }}
   return(out)    
  }


#' Compute posterior predictive trial summaries

#' @param x object of class pp_eos
#' @param threshold list of thresholds
#' @return list of summaries
#' @export

comp_pp_eos <- function(x, threshold=list(p=0.05, hr=0.8)){

 if( !inherits(x,"pp_eos") ){
  stop("x must be of class pp_eos")
 }

 eos <- attributes(x)[["eos"]]
 eos_units <- attributes(x)[["tunits"]]

 mp <- summary(x[1,])
 mhr <- summary(x[2,])
 m1 <- summary(x[3,])
 m2 <- summary(x[4,])
 md <- summary(x[4,]-x[3,])

 tlabs <- gsub("trt\\=","",rownames(x)[3:4])

 out <- list( end_of_study = paste0(eos," ",eos_units),
              thresholds = threshold,
              trial_time = attributes(x)[["summary"]][[1]] )

 out[[ "prob of stat signif" ]] <- sum( x["p",] < threshold$p )/ncol(x) # prob of stat signif at threshold
 out[[ "prob of clinical signif" ]] <- sum( x["hr",] < threshold$hr )/ncol(x) # prob of clinical signif at threshold

 out[[ paste0("median survival difference ",tlabs[[2]]," - ",tlabs[[1]]) ]] <- m2[[3]] - m1[[3]]

 out[["pval distr"]] <- mp
 out[["hr distr"]] <- mhr
 out[[ paste0("median survival ",tlabs[[1]]) ]] <- m1
 out[[ paste0("median survival ",tlabs[[2]]) ]] <- m2
 out[[ paste0("median difference distr ",tlabs[[2]]," - ",tlabs[[1]]) ]] <- md

 x_0 <- stats::na.omit(x[3:4,])
 
 out[[ paste0("prob median survival ",tlabs[[2]]," > ",tlabs[[1]]) ]] <- sum(x_0[2,]>x_0[1,],na.rm=TRUE)/ncol(x_0)

 out[["eos_stats"]] <- attributes(x)[["eos_stats"]]

 out[["eos_concord"]] <- attributes(x)[["eos_concord"]]

 out[["type"]] <- attributes(x)[["type"]]

 if( length(attributes(x)[["mod_formula_resp"]])>0 ){
  out[["mod_formula_resp"]] <- attributes(x)[["mod_formula_resp"]]
 }
 if( length(attributes(x)[["mod_formula_surv"]])>0 ){
  out[["mod_formula_surv"]] <- attributes(x)[["mod_formula_surv"]]
 }
 if( length(attributes(x)[["mod_formula"]])>0 ){
  out[["mod_formula"]] <- attributes(x)[["mod_formula"]]
 }
 
 out[["alternative"]] <- attributes(x)[["alternative"]]

 return(out)
}


#' Summarize posterior predictive trial summaries

#' @param x object of class pp_eos
#' @param threshold list of thresholds
#' @return list of summaries
#' @export

summary_pp_eos <- function(x, threshold=list(p=0.05, hr=0.8)){

 if( !inherits(x,"pp_eos") ){
  stop("x must be of class pp_eos")
 }

 print( attributes(x)[["posterior summaries"]] )

 print( plot_pp_eos( x=x, threshold=threshold ) )

 return(comp_pp_eos( x=x, threshold=threshold ))
}


#' Summarize concordance with eos stats 

concord_eos <- function(x,y,a){

 out <- c(  OS.med.p.dist_l=0,
            OS.med.p.dist_u=0, 
            OS.med.p.pch=NA, 
            OS.IQR.p.dist=NA, 
            OS.IQR.p.pch_l=NA, 
            OS.IQR.p.pch_u=NA, 

            OS.med.hr.dist_l=0,
            OS.med.hr.dist_u=0, 
            OS.med.hr.pch=NA, 
            OS.IQR.hr.dist=NA, 
            OS.IQR.hr.pch_l=NA, 
            OS.IQR.hr.pch_u=NA,

            OS.med.med.d.dist_l=0,
            OS.med.med.d.dist_u=0, 
            OS.med.med.d.pch=NA, 
            OS.IQR.med.d.dist=NA, 
            OS.IQR.med.d.pch_l=NA, 
            OS.IQR.med.d.pch_u=NA, 

            PFS.med.p.dist_l=0,
            PFS.med.p.dist_u=0, 
            PFS.med.p.pch=NA, 
            PFS.IQR.p.dist=NA, 
            PFS.IQR.p.pch_l=NA, 
            PFS.IQR.p.pch_u=NA, 

            PFS.med.hr.dist_l=0,
            PFS.med.hr.dist_u=0, 
            PFS.med.hr.pch=NA, 
            PFS.IQR.hr.dist=NA, 
            PFS.IQR.hr.pch_l=NA, 
            PFS.IQR.hr.pch_u=NA,

            PFS.med.med.d.dist_l=0,
            PFS.med.med.d.dist_u=0, 
            PFS.med.med.d.pch=NA, 
            PFS.IQR.med.d.dist=NA, 
            PFS.IQR.med.d.pch_l=NA, 
            PFS.IQR.med.d.pch_u=NA 

            )

 if( length(y)>0 ){

  if( y$OS[[a]] < median(x[1,],na.rm=TRUE) ){

   out[["OS.med.p.dist_u"]] <- median(x[1,],na.rm=TRUE) - y$OS[[a]]

  } else if( y$OS[[a]] > median(x[1,],na.rm=TRUE) ){

   out[["OS.med.p.dist_l"]] <- y$OS[[a]] - median(x[1,],na.rm=TRUE)   

  }

  out[["OS.med.p.pch"]] <- 100*( (median(x[1,],na.rm=TRUE)/y$OS[[a]]) - 1)

  out[["OS.IQR.p.dist"]] <- max( abs(y$OS[[a]]-quantile(x[1,],c(0.25,0.75),na.rm=TRUE)), na.rm=TRUE  )

  out[["OS.IQR.p.pch_l"]] <- 100*( (quantile(x[1,],c(0.25),na.rm=TRUE)/y$OS[[a]]) - 1)

  out[["OS.IQR.p.pch_u"]] <- 100*( (quantile(x[1,],c(0.75),na.rm=TRUE)/y$OS[[a]]) - 1)


  if( y$OS[["hr"]] < median(x[2,],na.rm=TRUE) ){

   out[["OS.med.hr.dist_u"]] <- median(x[2,],na.rm=TRUE) - y$OS[["hr"]]

  } else if( y$OS[["hr"]] > median(x[2,],na.rm=TRUE) ){

   out[["OS.med.hr.dist_l"]] <- y$OS[["hr"]] - median(x[2,],na.rm=TRUE)   

  }

  out[["OS.med.hr.pch"]] <- 100*( (median(x[2,],na.rm=TRUE)/y$OS[["hr"]]) - 1)

  out[["OS.IQR.hr.dist"]] <- max( abs(y$OS[["hr"]]-quantile(x[2,],c(0.25,0.75),na.rm=TRUE)), na.rm=TRUE  )

  out[["OS.IQR.hr.pch_l"]] <- 100*( (quantile(x[2,],c(0.25),na.rm=TRUE)/y$OS[["hr"]]) - 1)

  out[["OS.IQR.hr.pch_u"]] <- 100*( (quantile(x[2,],c(0.75),na.rm=TRUE)/y$OS[["hr"]]) - 1)


  #if( abs(y$OS$medians[[2]]-y$OS$medians[[1]]) < abs(median(x[4,],na.rm=TRUE)-median(x[3,],na.rm=TRUE)) ){
  #
  # out[["OS.med.med.d.dist_u"]] <- abs(median(x[4,],na.rm=TRUE)-median(x[3,],na.rm=TRUE)) - abs(y$OS$medians[[2]]-y$OS$medians[[1]])
  #
  #} else if( abs(y$OS$medians[[2]]-y$OS$medians[[1]]) > abs(median(x[4,],na.rm=TRUE)-median(x[3,],na.rm=TRUE)) ){
  #
  # out[["OS.med.med.d.dist_l"]] <- abs(y$OS$medians[[2]]-y$OS$medians[[1]]) - abs(median(x[4,],na.rm=TRUE)-median(x[3,],na.rm=TRUE)) 
  #
  #}
  #
  #out[["OS.med.med.d.pch"]] <- 100*( (abs((median(x[4,],na.rm=TRUE)-median(x[3,],na.rm=TRUE)))/abs((y$OS$medians[[2]]-y$OS$medians[[1]]))) - 1)
  #
  #out[["OS.IQR.med.d.dist"]] <- max( abs(abs((y$OS$medians[[2]]-y$OS$medians[[1]]))-quantile(abs(x[4,]-x[3,]),c(0.25,0.75),na.rm=TRUE)), na.rm=TRUE  )
  #
  #out[["OS.IQR.med.d.pch_l"]] <- 100*( (quantile(abs(x[4,]-x[3,]),c(0.25),na.rm=TRUE)/abs((y$OS$medians[[2]]-y$OS$medians[[1]]))) - 1)
  #
  #out[["OS.IQR.med.d.pch_u"]] <- 100*( (quantile(abs(x[4,]-x[3,]),c(0.75),na.rm=TRUE)/abs((y$OS$medians[[2]]-y$OS$medians[[1]]))) - 1)

 }

 return(out)
}


