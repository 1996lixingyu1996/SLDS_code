#' General class for Weibull Generated multistate data
#'
#' @export
setClass("mpti_weib_dat")

#' General class for Weibull Generated multistate interim data
#'
#' @export
setClass("impti_weib_dat")

#' General Multistate data class
#'
#' A class representing a multistate data objects not simulated from a multistate model
#'
#' @slot data data.frame of class mpti_weib_dat containing patient-level outcomes and treatment
#' @slot cntrl name of the unique subject identifier
#' @slot exper name of the duration since study start to accrual (var must be positive numeric value)
#' @slot resp_levels vector of stratum of resp variable that indicate a "response" (e.g. RECIST PR or CR)
#' @slot resp_t logical if variable name of duration from start of study to "response" is included
#' @slot covar name of prognostic covariates
#' @slot ttp logical if time to progression duration variable is included
#' @slot dor logical if duration of response variable is included
#' @import methods
#'
#' @export
setClass("mstate_data",
         slots = list(
           data = "mpti_weib_dat",
           cntrl = "character",
           exper = "character",
           resp_levels = "character",
           resp_t = "logical",
           covar = "character",
           ttp = "logical",
           dor = "logical"
         )
)

#' Create a Multistate data class Object
#'
#' Create a new instance of the \code{mstate_dat} class.
#'
#' @param data data.frame containing patient-level outcomes and treatment
#' @param resp name of the surrogate response variable (var is factor)
#' @param resp_levels vector of stratum of resp variable that indicate a "response" (e.g. RECIST PR or CR)
#' @param surv name of overall survival duration variable (var is numeric)
#' @param surv_evnt name of overall survival event indicator (integer 0-1)
#' @param surv_csnr name of overall survival right-censoring indicator (integer 0-1)
#' @param pfs name of progression-free survival duration variable (var is numeric)
#' @param pfs_evnt name of progression-free survival event indicator (integer 0-1)
#' @param pfs_csnr name of progression-free survival right-censoring indicator (integer 0-1)
#' @param accT name of the duration since study start to accrual (var must be positive numeric value)
#' @param id name of the unique subject identifier
#' @param trt name of the treatment indicator variable (var can be integer or factor)
#' @param cntrl name of the control arm group (in treatment variable)
#' @param exper name of the experimental arm group (in treatment variable)
#' @param resp_t name of the variable describing the time of observed response
#' @param covar name of prognostic covariates
#' @param ttp name of time-to-progression duration variable (var is numeric; needed only if missing PFS & Survival)
#' @param ttp_evnt name of time-to-progression event indicator (integer 0-1; needed only if missing PFS & Survival)
#' @param ttp_csnr name of time-to-progression right-censoring indicator (integer 0-1; needed only if missing PFS & Survival)
#' @param cycle_dur_weeks duration of first cycle in weeks
#' @param tunits label for time units
#' @param accTime time parition for accrual
#' @param accExN expected enrollment for each bin of accTime
#' @param addFUt added follow-up duration in months before final analysis
#' @return An object of class \code{mstate_dat}.
#' @importFrom tibble is_tibble
#' @importFrom stats relevel runif
#' @export
create_mstate <- function( data, resp, resp_levels, surv, surv_evnt=NA, surv_csnr=NA, pfs, pfs_evnt=NA, pfs_csnr=NA,
                           accT, id=NA, trt=NA, cntrl=NA, exper=NA, resp_t=NA, covar=NA,
                           ttp=NA, ttp_evnt=NA, ttp_csnr=NA, cycle_dur_weeks=4,
                           tunits="months", accTime=NULL, accExN=NULL, addFUt=NULL) {

  if( tibble::is_tibble(data) ){
   data <- as.data.frame(data)
  }

  if( !is.data.frame(data) ){
   stop("data must be class data.frame or tibble")
  }

  if( !(resp %in% names(data) ) ){
   stop("resp not found in data")
  }

  if( !is.factor(data[,resp]) ){
   stop("resp must be factor class")
  }

  if( !is.character(resp_levels) ){
   resp_levels <- as.character(resp_levels)
  }

  if( sum(!(resp_levels %in% levels(data[,resp])))>0 ){
   stop("resp_level not found in levels of resp variable")
  }

  if( !(surv %in% names(data)) ){
   stop("surv not found in data")
  }

  if( !(pfs %in% names(data)) ){
   stop("pfs not found in data")
  }  

  if( !is.numeric(data[,surv]) ){
   stop("survival must be numeric class")
  }

  if( !is.numeric(data[,pfs]) ){
   stop("pfs must be numeric class")
  }  

  if( min(data[,surv])<0 ){
   stop("survival must be positive")
  }

  if( min(data[,pfs])<0 ){
   stop("pfs must be positive")
  }

  if( is.na(surv_evnt) && is.na(surv_csnr) ){
   
   stop("both survival event & censoring indicators are NA")

  } else if( is.na(surv_csnr) ){
   
   if( !(surv_evnt %in% names(data)) ){
    stop("survival event indicator not found in data")
   }

   if( is.numeric(data[,surv_evnt]) ){

    if( sum(!(unique(data[,surv_evnt]) %in% 0:1))>0 ){
     stop("numeric survival event indicator must be integer valued with 0 & 1")
    }

   } else{

    if( sum(!(unique(as.character(data[,surv_evnt])) %in% as.character(0:1)))>0 ){
     stop("non-numeric survival event indicator must contain character 0 & 1")
    }

   }

  } else{

   if( !(surv_csnr %in% names(data)) ){
    stop("survival censor indicator not found in data")
   }

   if( is.numeric(data[,surv_csnr]) ){

    if( sum(!(unique(data[,surv_csnr]) %in% 0:1))>0 ){
     stop("numeric survival censoring indicator must be integer valued with 0 & 1")
    }

   } else{

    if( sum(!(unique(as.character(data[,surv_csnr])) %in% as.character(0:1)))>0 ){
     stop("non-numeric survival censoring indicator must contain character 0 & 1")
    }

   }

  }

  if( is.na(pfs_evnt) && is.na(pfs_csnr) ){
   
   stop("both survival event & censoring indicators are NA")

  } else if( is.na(pfs_csnr) ){
   
   if( !(pfs_evnt %in% names(data)) ){
    stop("survival event indicator not found in data")
   }

   if( is.numeric(data[,pfs_evnt]) ){

    if( sum(!(unique(data[,pfs_evnt]) %in% 0:1))>0 ){
     stop("numeric survival event indicator must be integer valued with 0 & 1")
    }

   } else{

    if( sum(!(unique(as.character(data[,pfs_evnt])) %in% as.character(0:1)))>0 ){
     stop("non-numeric survival event indicator must contain character 0 & 1")
    }

   }

  } else{

   if( !(pfs_csnr %in% names(data)) ){
    stop("survival censor indicator not found in data")
   }

   if( is.numeric(data[,pfs_csnr]) ){

    if( sum(!(unique(data[,pfs_csnr]) %in% 0:1))>0 ){
     stop("numeric survival censoring indicator must be integer valued with 0 & 1")
    }

   } else{

    if( sum(!(unique(as.character(data[,pfs_csnr])) %in% as.character(0:1)))>0 ){
     stop("non-numeric survival censoring indicator must contain character 0 & 1")
    }

   }

  }

  if( !(accT %in% names(data)) ){

   stop("accrual duration variable not found")

  } else if( !is.numeric(data[,accT]) ){

   data[,accT] <- as.numeric(data[,accT])

  }

  if( !is.numeric(data[,accT]) ){ stop("accrual duration variable must be numeric") }

  if( min(data[,accT])<0 ){ stop("accrual duration variable must be positive") }

  if(!is.na(trt)){ 

   if( trt %in% names(data) ){

    if( length(unique(data[,trt]))>2 ){

     if( is.na(cntrl) && is.na(exper) ){
      stop("cntrl and exper must be specified if treatment has >2 levels")
     }

     if( !(cntrl %in% unique(data[,trt])) || !(exper %in% unique(data[,trt])) ){
      stop("cntrl or exper not found in treatment")
     }

     data <- base::subset( data, data[,trt] %in% c(cntrl,exper) )

    } 
   
    if( !is.factor(data[,trt]) ){ 
     data[,trt] <- as.factor(as.character(data[,trt]))
     if( !is.na(cntrl) ){
      if( cntrl %in% levels(data[,trt]) ){ 
       data[,trt] <- stats::relevel(data[,trt], ref=cntrl) 
      } else{
       stop("cntrl level not found in treatment variable")
      }
     }
    }

    if( !is.na(exper) ){
     if( !(exper %in% levels(data[,trt])) ){
      stop("exper level not found in treatment variable")
     }
    }

    if( length(levels(data[,trt]))==2 ){
   
     cntrl <- levels(data[,trt])[[1]]
     exper <- levels(data[,trt])[[2]]
   
    } else{

     exper <- levels(data[,trt])[[1]]
     cntrl <- NA

    }

   } else{
    stop("trt not found in data")
   }

  }
  ###########

  out <- data

  keep_v <- NULL

  if( !is.na(id) ){ 
   if( id %in% names(data) ){ 
    if( !is.character(data[,id]) ){
     data[,id] <- as.character(data[,id])
    }
   } else{
    stop("id variable not found in data")
   }
   out[,"USUBJID"] <- data[,id]
  } else{
   out[,"USUBJID"] <- as.character(1:nrow(out))
  }
  keep_v <- c(keep_v,"USUBJID")

  out[,"accT"] <- data[,accT]
  keep_v <- c(keep_v,"accT")

  out[,"resp"] <- NA
  w <- which( as.character(data[,resp]) %in% resp_levels )
  if( length(w)>0 ){
   out[w,"resp"] <- 1
  }
  w <- which( !(as.character(data[,resp]) %in% resp_levels) )
  if( length(w)>0 ){
   out[w,"resp"] <- 0
  }  
  keep_v <- c(keep_v,"resp")

  out[,"OS"] <- data[,surv]
  if( is.na(surv_evnt) ){
   out[,"OSevent"] <-  1-as.numeric(as.character(data[,surv_csnr]))
  } else{
   out[,"OSevent"] <-  as.numeric(as.character(data[,surv_evnt]))
  }
  keep_v <- c(keep_v,"OS","OSevent")

  out[,"trt"] <- NA

  if(!is.na(trt)){

   if( length(unique(data[,trt]))==2 ){
    w <- which( as.character(data[,trt]) == levels(data[,trt])[[1]] )
    if( length(w)>0 ){
     out[w,"trt"] <- 0
    }
    w <- which( as.character(data[,trt]) == levels(data[,trt])[[2]] )
    if( length(w)>0 ){
     out[w,"trt"] <- 1
    }  
   } else{
    out[,"trt"] <- 1
   }

  } else{

   out[,"trt"] <- 1

   if( is.na(exper) && !is.na(cntrl) ){
    exper <- cntrl
    cntrl <- NA
   } else if( is.na(exper) && is.na(cntrl) ){
    exper <- "treatment"
   }  

  }
  keep_v <- c(keep_v,"trt")

  dur_pfs <- NULL; evnt_pfs <- NULL
  if( !is.na(pfs) ){ if( pfs %in% names(data) ){ if( is.numeric(data[,pfs]) ){ if( min(data[,pfs])>=0 ){
   if( !(is.na(pfs_evnt) && is.na(pfs_csnr)) ){
    if( is.na(pfs_csnr) ){
     if( pfs_evnt %in% names(data) ){
      if( sum(!(unique(as.character(data[,pfs_evnt])) %in% as.character(0:1)))==0 ){
       evnt_pfs <- as.numeric(as.character(data[,pfs_evnt]))
       dur_pfs <- data[,pfs] 
      }
     }
    } else{
     if( pfs_csnr %in% names(data) ){
      if( sum(!(unique(as.character(data[,pfs_csnr])) %in% as.character(0:1)))==0 ){
       evnt_pfs <- 1-as.numeric(as.character(data[,pfs_csnr]))
       dur_pfs <- data[,pfs] 
      }
     }      
    }
   }
  }}}}

  if( length(dur_pfs)>0 && length(evnt_pfs)>0 ){
   out[,"PFS"] <- dur_pfs
   out[,"PFSevent"] <- evnt_pfs
   keep_v <- c(keep_v,"PFS","PFSevent")
  }

  if( "PFS" %in% keep_v ){

   out[,"TTP"] <- out[,"PFS"]
   out[,"TTPevent"] <- out[,"PFSevent"]
   for( ii in 1:nrow(out) ){
    if( out[ii,"PFSevent"]==1 ){
     if( ( out[ii,"PFS"]==out[ii,"OS"] ) && ( out[ii,"OSevent"]==1 ) ){
      out[ii,"TTPevent"] <- 0
     }
    }
   }  
   keep_v <- c(keep_v,"TTP","TTPevent")

  } else{

   dur_ttp <- NULL; evnt_ttp <- NULL
   if( !is.na(ttp) ){ if( ttp %in% names(data) ){ if( is.numeric(data[,ttp]) ){ if( min(data[,ttp])>=0 ){
    if( !(is.na(ttp_evnt) && is.na(ttp_csnr)) ){
     if( is.na(ttp_csnr) ){
      if( ttp_evnt %in% names(data) ){
       if( sum(!(unique(as.character(data[,ttp_evnt])) %in% as.character(0:1)))==0 ){
        evnt_ttp <- as.numeric(as.character(data[,ttp_evnt]))
        dur_ttp <- data[,ttp] 
       }
      }
     } else{
      if( ttp_csnr %in% names(data) ){
       if( sum(!(unique(as.character(data[,ttp_csnr])) %in% as.character(0:1)))==0 ){
        evnt_ttp <- 1-as.numeric(as.character(data[,ttp_csnr]))
        dur_ttp <- data[,ttp] 
       }
      }      
     }
    }
   }}}}

   if( length(dur_ttp)>0 && length(evnt_ttp)>0 ){
    out[,"TTP"] <- dur_ttp
    out[,"TTPevent"] <- evnt_ttp
    keep_v <- c(keep_v,"TTP","TTPevent")
   }

   out[,"PFS"] <- out[,"TTP"]
   out[,"PFSevent"] <- out[,"TTPevent"]
   for( ii in 1:nrow(out) ){
    if( out[ii,"TTPevent"]==0 ){
     if( ( out[ii,"TTP"]==out[ii,"OS"] ) && ( out[ii,"OSevent"]==1 ) ){
      out[ii,"PFSevent"] <- 1
     }
    }
   }
   keep_v <- c(keep_v,"PFS","PFSevent")

  }

  if( "PFS" %in% keep_v ){

   out[,"DOR"] <- out[,"PFS"]
   out[,"DORevent"] <- out[,"PFSevent"]
   w <- which(out[,"resp"]==0)
   if( length(w)>0 ){
    out[w,"DOR"] <- NA
    out[w,"DORevent"] <- NA
   }
   keep_v <- c(keep_v,"DOR","DORevent")   

  } 

  out[,"resp_t"] <- NA
  if( !is.na(resp_t) ){ 
   if( resp_t %in% names(data) ){ if( is.numeric(data[,resp_t]) ){
    out[,"resp_t"] <- data[,resp_t]
   }}
  } else if( sum(out$resp)>0 && tunits=="months"){  ### implement random response time for months ###
   out[which(out$resp==1), "resp_t"] <- stats::runif(sum(out$resp==1),(cycle_dur_weeks*7)/(365.25/12),out[which(out$resp==1),"PFS"]*0.67)
  }
  keep_v <- c(keep_v,"resp_t")

  #### covariates ####
  covar <- NA

  if( !(length(covar)==1 && is.na(covar)[[1]]) ){
   covar <- base::subset(covar, covar %in% names(data) )
   if( length(covar)>0 ){ for(ii in covar){

    if( is.numeric(data[,ii]) ){
     out[,ii] <- data[,ii]
     keep_v <- c(keep_v,ii)
    } else if( is.factor(data[,ii]) ){
     if( length(levels(data[,ii]))<=5 ){
      out[,ii] <- data[,ii]
      keep_v <- c(keep_v,ii)
     }
    } else if(is.character(data[,ii])){
     if( length(levels(as.factor(data[,ii])))<=5 ){
      out[,ii] <- as.factor(data[,ii])
      keep_v <- c(keep_v,ii)
     }
    }

   }}
  
  }
  
  out <- out[,keep_v]

  attributes(out)[["type"]] <- "user"
  attributes(out)[["scenario"]] <- NA
  attributes(out)[["eos"]] <- round(max(out$OS,na.rm=TRUE))
  attributes(out)[["tunits"]] <- tunits

  #### covariates not implemented ####
  nb <- nc <- 0

  if(FALSE){
   if( !(length(covar)==1 && !is.na(covar)[[1]]) ){
    for(ii in covar){

     if( is.numeric(out[,ii]) ){
      
       if( length(unique(out[,ii]))==2 ){
        nb <- nb + 1
       } else{
        nc <- nc + 1
       }
    
     } else{
      nb <- nb + as.numeric( is.numeric(out[,ii]) )
     }
   
    }
   }
  }

  if( length(accTime)==1 ){

   if( is.numeric(accTime) ){ if( accTime < 0 ){
     accTime <- NULL
   }}

   if( length(accTime)<1 || !is.numeric(accTime) ){
    accTime <- round(2*max(out$OS,na.rm=TRUE)/3)
   }

   if( accTime > attributes(out)[["eos"]] - addFUt ){
    accTime <- attributes(out)[["eos"]] - addFUt
   }

  } else{

   if( length(accTime)<1 ){
    
    accTime <- round(2*max(out$OS,na.rm=TRUE)/3)
   
   } else if( length(accTime)!=length(accExN) ){

     accTime <- round(2*max(out$OS,na.rm=TRUE)/3)
     accExN <- NULL

   }

  }

  if( length(addFUt)==0){
   addFUt <- 0
  }

  if( length(accTime)==1 ){

   if( accTime > attributes(out)[["eos"]] - addFUt ){

    accTime <- attributes(out)[["eos"]] - addFUt

   } else{

    addFUt <- attributes(out)[["eos"]] - accTime

   }

  }

  #### Add multi-state pathway info ####
  class(out) <- c("mpti_weib_dat",class(out))

  out[,"enrolled"] <- TRUE

  out <- add_pathlos_data(dat=out)

  attributes(out)[["args"]] <- list(
                       n=nrow(out), 
                       nb=nb, 
                       nc=nc, 
                       prtrt=round(sum(out$trt==1)/nrow(out),3), 
                       accTime=accTime, 
                       accExN=accExN, 
                       addFUt=addFUt, 
                       CenUpLim=attributes(out)[["eos"]], 
                       tunits=tunits )


  new("mstate_data",
      data =  out,
      cntrl = as.character(cntrl),
      exper = as.character(exper),
      resp_levels = as.character(resp_levels),
      resp_t = sum(!is.na(out[,"resp_t"]))>0,
      covar = as.character(covar),
      ttp = ifelse("TTP" %in% names(out), TRUE, FALSE),
      dor = ifelse("DOR" %in% names(out), TRUE, FALSE))
}

#' General Multistate interim data class
#'
#' A class representing a multistate data objects not simulated from a multistate model
#'
#' @slot data data.frame of class impti_weib_dat containing patient-level outcomes and treatment
#' @slot cntrl name of the unique subject identifier
#' @slot exper name of the duration since study start to accrual (var must be positive numeric value)
#' @slot resp_levels vector of stratum of resp variable that indicate a "response" (e.g. RECIST PR or CR)
#' @slot resp_t logical if variable name of duration from start of study to "response" is included
#' @slot covar name of prognostic covariates
#' @slot pfs logical if progression-free survival duration variable is included
#' @slot ttp logical if time to progression duration variable is included
#' @slot dor logical if duration of response variable is included
#' @import methods
#'
#' @export
setClass("imstate_data",
         slots = list(
           data = "impti_weib_dat",
           cntrl = "character",
           exper = "character",
           resp_levels = "character",
           resp_t = "logical",
           covar = "character",
           ttp = "logical",
           dor = "logical"
         )
)
