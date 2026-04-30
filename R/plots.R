
#' Plotting function for survival mediation modeling
#'
#' @param x list of class riskratio
#' @param type text string indicating type of plot "riskratio", "survival", or "medprop"
#' @return plot log risk ratios
#' @export

plot_rr <- function(x,type="riskratio"){
 if( !inherits(x,"riskratio_tru") && !inherits(x,"riskratio_est") ){
  stop("x must be obj of class riskratio")
 }
 if( type=="riskratio" ){ return( plot_rr_riskratio(x) ) 
 } else if( type=="survival" ){ return( plot_rr_survival(x) ) 
 } else if( type=="medprop" ){ return( plot_rr_medprop(x) ) 
 }
}

# Plot log risk ratio from mediation analysis

plot_rr_riskratio <- function(x){
 if( inherits(x,"riskratio_tru") ){
  return( plot_rr_riskratio_tru(x) )
 } else if( inherits(x,"riskratio_est") ){
  return( plot_rr_riskratio_est(x) )
 } else{
  return(NULL)
 }
}

# Plot survival functions from mediation analysis

plot_rr_survival <- function(x){
 if( inherits(x,"riskratio_tru") ){
  return( plot_rr_survival_tru(x) )
 } else if( inherits(x,"riskratio_est") ){
  return( plot_rr_survival_est(x) )
 } else{
  return(NULL)
 }
}

# Plot survival functions from mediation analysis

plot_rr_medprop <- function(x){
 if( inherits(x,"riskratio_tru") ){
  return( plot_rr_medprop_tru(x) )
 } else if( inherits(x,"riskratio_est") ){
  return( plot_rr_medprop_est(x) )
 } else{
  return(NULL)
 }
}


# Plot log risk ratio from mediation analysis given parms

#' @import ggplot2
#' @import patchwork

plot_rr_riskratio_tru <- function(x){

dat0 <- as.data.frame( list(time=x$tseq, y=NA, effect=NA) )
datd <- dat0; datd$y <- x$DIF_d; datd$effect <- "Direct"
datm <- dat0; datm$y <- x$DIF_m; datm$effect <- "Mediated"
datt <- dat0; datt$y <- x$DIF_tot; datt$effect <- "Total"
dat <- rbind( datd, datm, datt )

pp <- ggplot(data=dat, aes(x=time, y=y, group = effect, color = effect)) + 
  geom_hline(yintercept = 0, color="gray") +
  geom_line( size=1.5) + ggplot2::theme_minimal() + 
  theme(
    plot.title = element_text(color = "black", face = "bold", size = 12),
    axis.title.x = element_text(color = "black", face = "bold", size = 12),
    axis.title.y = element_text(color = "black", face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12)) +
  labs(
    y = paste0("log Risk Ratios", "\n"),
    x = paste0("\n","time"),
    title = "True log Risk Ratios",
    subtitle = paste0("Scenario = ",attributes(x)$scenario) )

df_means <- as.data.frame( list(mean=c(attributes(x)$pr.ctr,attributes(x)$pr.trt),arm=c("control","treat")) )
lower_bound <- min(c(attributes(x)$pr.ctr,attributes(x)$pr.trt)); upper_bound <- max(c(attributes(x)$pr.ctr,attributes(x)$pr.trt))
distance <- round(abs(attributes(x)$pr.ctr - attributes(x)$pr.trt), 2)
midpoint <- max(c(attributes(x)$pr.ctr,attributes(x)$pr.trt)) - distance/2
    
    p <- ggplot2::ggplot(df_means, mapping = ggplot2::aes(x = 0, y = mean, group=arm, color=arm)) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      geom_segment( aes(y=lower_bound, x=0, xend=0, yend=upper_bound, linewidth=6), color="darkgray", show.legend=FALSE ) +
      ggplot2::geom_point(alpha =.8, size = 6) +
      ggplot2::guides( fill = "none" ) +
      scale_color_manual(values = c("#E15759", "#4E79A7"), name = "arm") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank(),
                    plot.title = element_text(color = "black", face = "bold", size = 12),
                    axis.title.x = element_text(color = "black", face = "bold", size = 12),
                    axis.title.y = element_text(color = "black", face = "bold", size = 12),
                    axis.text.x = element_text(color = "black", face = "bold", size = 12) ) +
  labs( y = paste0("\n","probability of response"), 
    x = "",    
    title = "Surrogate response rates",
    subtitle = paste0("difference = ",distance) )                     
    
 out <- pp / p  + plot_layout(heights = c(6, 1))

 return(out)
}


# Plot log risk ratio from mediation analysis

#' @import ggplot2
#' @import patchwork

plot_rr_riskratio_est <- function(x){

dat0 <- as.data.frame( list(time=x$tseq, y=NA, effect=NA) )
datd <- dat0; datd$y <- x$DIF_d; datd$effect <- "Direct"
datm <- dat0; datm$y <- x$DIF_m; datm$effect <- "Mediated"
datt <- dat0; datt$y <- x$DIF_tot; datt$effect <- "Total"
dat <- rbind( datd, datm, datt )

pp <- ggplot(data=dat, aes(x=time, y=y, group = effect, color = effect)) + 
  geom_hline(yintercept = 0, color="gray") +
  geom_line( size=1.5) + ggplot2::theme_minimal() + 
  theme(
    plot.title = element_text(color = "black", face = "bold", size = 12),
    axis.title.x = element_text(color = "black", face = "bold", size = 12),
    axis.title.y = element_text(color = "black", face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12))
  if( length(attributes(x)$scenario)>0 ){ 
   pp <- pp + labs(
    y = paste0("log Risk Ratios", "\n"),
    x = paste0("\n","time"),
    title = "Estimated log Risk Ratios",
    subtitle = paste0("Generated from true scenario = ",attributes(x)$scenario) )
  } else{
   pp <- pp + labs(
    y = paste0("log Risk Ratios", "\n"),
    x = paste0("\n","time"),
    title = "Estimated log Risk Ratios")
  }

 x1 <- stats::rbeta(10000,attributes(x)$beta.ctr[[1]],attributes(x)$beta.ctr[[2]])
 x2 <- stats::rbeta(10000,attributes(x)$beta.trt[[1]],attributes(x)$beta.trt[[2]])
 
 df <- as.data.frame( list(x=c(x1, x2), arm=factor(c(rep("control",length(x1)),rep("treated",length(x2)))) ))

 p <- ggplot2::ggplot(df, aes(color=arm, fill=arm) ) +  
      ggplot2::geom_density(mapping = ggplot2::aes(x = x, group = arm), linetype=0, alpha = 0.2) +
                             geom_vline(xintercept = mean(x1), color="darkgray") +
                             geom_vline(xintercept = mean(x2), color="darkgray") +
                             ggplot2::theme_minimal() +
                             ggplot2::theme(
                               panel.grid.major.y = ggplot2::element_blank(),
                               panel.grid.minor.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank(),
                               axis.text.y=element_blank(),
                               plot.title = ggplot2::element_text( color = "black", face = "bold", size = 12),
                               axis.title.x = ggplot2::element_text( color = "black", face = "bold", size = 12),
                               axis.text.x = ggplot2::element_text( color = "black", face = "bold", size = 12)
                             ) +
                             scale_fill_manual(values = c("#E15759","#4E79A7") ) +
                             scale_color_manual(values = c("#E15759","#4E79A7") ) +
                             ggplot2::labs(
                               x = paste0("\n", "response rate"),
                               y = "",
                               title = paste0("Surrogate response rate posterior distributions")
                              )

 out <- pp / p  + plot_layout(heights = c(6, 1))

 return(out)
}



# Plot mediation proportion from mediation analysis given parms

#' @import ggplot2
#' @import patchwork

plot_rr_medprop_tru <- function(x){

dat <- as.data.frame( list(time=x$tseq, y=x$PROP_m) )

pp <- ggplot(data=dat, aes(x=time, y=y)) + 
  geom_hline(yintercept = 0, color="darkgray") +
  geom_hline(yintercept = 1, color="darkgray") +
  geom_line( size=1.5, color="navy" ) + ggplot2::theme_minimal() + 
  theme(
    plot.title = element_text(color = "black", face = "bold", size = 12),
    axis.title.x = element_text(color = "black", face = "bold", size = 12),
    axis.title.y = element_text(color = "black", face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12)) +
  labs(
    y = paste0("Mediation proportion", "\n"),
    x = paste0("\n","time"),
    title = "True Mediation proportion",
    subtitle = paste0("Scenario = ",attributes(x)$scenario) )

df_means <- as.data.frame( list(mean=c(attributes(x)$pr.ctr,attributes(x)$pr.trt),arm=c("control","treat")) )
lower_bound <- min(c(attributes(x)$pr.ctr,attributes(x)$pr.trt)); upper_bound <- max(c(attributes(x)$pr.ctr,attributes(x)$pr.trt))
distance <- round(abs(attributes(x)$pr.ctr - attributes(x)$pr.trt), 2)
midpoint <- max(c(attributes(x)$pr.ctr,attributes(x)$pr.trt)) - distance/2
    
    p <- ggplot2::ggplot(df_means, mapping = ggplot2::aes(x = 0, y = mean, group=arm, color=arm)) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      geom_segment( aes(y=lower_bound, x=0, xend=0, yend=upper_bound, linewidth=6), color="darkgray", show.legend=FALSE ) +
      ggplot2::geom_point(alpha =.8, size = 6) +
      ggplot2::guides( fill = "none" ) +
      scale_color_manual(values = c("#E15759", "#4E79A7"), name = "arm") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank(),
                    plot.title = element_text(color = "black", face = "bold", size = 12),
                    axis.title.x = element_text(color = "black", face = "bold", size = 12),
                    axis.title.y = element_text(color = "black", face = "bold", size = 12),
                    axis.text.x = element_text(color = "black", face = "bold", size = 12) ) +
  labs( y = paste0("\n","probability of response"), 
    x = "",    
    title = "Surrogate response rates",
    subtitle = paste0("difference = ",distance) )                     
    
 out <- pp / p  + plot_layout(heights = c(6, 1))

 return(out)
}


# Plot log risk ratio from mediation analysis

#' @import ggplot2
#' @import patchwork

plot_rr_medprop_est <- function(x){

dat <- as.data.frame( list(time=x$tseq, y=x$PROP_m) )

pp <- ggplot(data=dat, aes(x=time, y=y)) + 
  geom_hline(yintercept = 0, color="darkgray") +
  geom_hline(yintercept = 1, color="darkgray") +
  geom_line( size=1.5, color="navy" ) + ggplot2::theme_minimal() + 
  theme(
    plot.title = element_text(color = "black", face = "bold", size = 12),
    axis.title.x = element_text(color = "black", face = "bold", size = 12),
    axis.title.y = element_text(color = "black", face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12))
  if( length(attributes(x)$scenario)>0 ){
   pp <- pp + labs(
    y = paste0("Mediation proportion", "\n"),
    x = paste0("\n","time"),
    title = "Estimated mediation proportion",
    subtitle = attributes(x)$scenario )
  } else{
   pp <- pp + labs(
    y = paste0("Mediation proportion", "\n"),
    x = paste0("\n","time"),
    title = "Estimated mediation proportion")    
  }

 x1 <- stats::rbeta(10000,attributes(x)$beta.ctr[[1]],attributes(x)$beta.ctr[[2]])
 x2 <- stats::rbeta(10000,attributes(x)$beta.trt[[1]],attributes(x)$beta.trt[[2]])
 
 df <- as.data.frame( list(x=c(x1, x2), arm=factor(c(rep("control",length(x1)),rep("treated",length(x2)))) ))

 p <- ggplot2::ggplot(df, aes(color=arm, fill=arm) ) +  
      ggplot2::geom_density(mapping = ggplot2::aes(x = x, group = arm), linetype=0, alpha = 0.2) +
                             geom_vline(xintercept = mean(x1), color="darkgray") +
                             geom_vline(xintercept = mean(x2), color="darkgray") +
                             ggplot2::theme_minimal() +
                             ggplot2::theme(
                               panel.grid.major.y = ggplot2::element_blank(),
                               panel.grid.minor.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank(),
                               axis.text.y=element_blank(),
                               plot.title = ggplot2::element_text( color = "black", face = "bold", size = 12),
                               axis.title.x = ggplot2::element_text( color = "black", face = "bold", size = 12),
                               axis.text.x = ggplot2::element_text( color = "black", face = "bold", size = 12)
                             ) +
                             scale_fill_manual(values = c("#E15759","#4E79A7") ) +
                             scale_color_manual(values = c("#E15759","#4E79A7") ) +
                             ggplot2::labs(
                               x = paste0("\n", "response rate"),
                               y = "",
                               title = paste0("Surrogate response rate posterior distributions")
                              )

 out <- pp / p  + plot_layout(heights = c(6, 1))

 return(out)
}



# Plot survival from mediation analysis given parms

#' @import ggplot2
#' @import patchwork

plot_rr_survival_tru <- function(x){

dat0 <- as.data.frame( list(time=x$tseq, y=NA, effect=NA) )
datd <- dat0; datd$y <- x$St.a0m0; datd$effect <- "Control"
datm <- dat0; datm$y <- x$St.a1m0; datm$effect <- "Counterfactual"
datt <- dat0; datt$y <- x$St.a1m1; datt$effect <- "Treated"
dat <- rbind( datd, datm, datt )

pp <- ggplot(data=dat, aes(x=time, y=y, group = effect, color = effect)) + 
  geom_hline(yintercept = 0, color="gray") +
  geom_line( size=1.5) + ggplot2::theme_minimal() + 
  geom_segment(aes(x=0, y=0.5, yend=0.5, xend=attributes(x)$summary[["median.ctr"]]), col="gray") + 
  geom_segment(aes(x=0, y=0.5, yend=0.5, xend=attributes(x)$summary[["median.trt"]]), col="gray") +  
  geom_segment(aes(x=attributes(x)$summary[["median.ctr"]], y=0, yend=0.5, xend=attributes(x)$summary[["median.ctr"]]), col="gray") +  
  geom_segment(aes(x=attributes(x)$summary[["median.trt"]], y=0, yend=0.5, xend=attributes(x)$summary[["median.trt"]]), col="gray") +  
  theme(
    plot.title = element_text(color = "black", face = "bold", size = 12),
    axis.title.x = element_text(color = "black", face = "bold", size = 12),
    axis.title.y = element_text(color = "black", face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12)) +
  labs(
    y = paste0("Survival", "\n"),
    x = paste0("\n","time"),
    title = "True Survival functions",
    subtitle = paste0("Scenario = ",attributes(x)$scenario) )

df_means <- as.data.frame( list(mean=c(attributes(x)$pr.ctr,attributes(x)$pr.trt),arm=c("control","treat")) )
lower_bound <- min(c(attributes(x)$pr.ctr,attributes(x)$pr.trt)); upper_bound <- max(c(attributes(x)$pr.ctr,attributes(x)$pr.trt))
distance <- round(abs(attributes(x)$pr.ctr - attributes(x)$pr.trt), 2)
midpoint <- max(c(attributes(x)$pr.ctr,attributes(x)$pr.trt)) - distance/2
    
    p <- ggplot2::ggplot(df_means, mapping = ggplot2::aes(x = 0, y = mean, group=arm, color=arm)) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      geom_segment( aes(y=lower_bound, x=0, xend=0, yend=upper_bound, linewidth=6), color="darkgray", show.legend=FALSE ) +
      ggplot2::geom_point(alpha =.8, size = 6) +
      ggplot2::guides( fill = "none" ) +
      scale_color_manual(values = c("#E15759", "#4E79A7"), name = "arm") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     panel.grid.minor.y = ggplot2::element_blank(),
                    plot.title = element_text(color = "black", face = "bold", size = 12),
                    axis.title.x = element_text(color = "black", face = "bold", size = 12),
                    axis.title.y = element_text(color = "black", face = "bold", size = 12),
                    axis.text.x = element_text(color = "black", face = "bold", size = 12) ) +
  labs( y = paste0("\n","probability of response"), 
    x = "",    
    title = "Surrogate response rates",
    subtitle = paste0("difference = ",distance) )                     
    
 out <- pp / p  + plot_layout(heights = c(6, 1))

 return(out)
}


# Plot survival from mediation analysis

#' @import ggplot2
#' @import patchwork

plot_rr_survival_est <- function(x){

dat0 <- as.data.frame( list(time=x$tseq, y=NA, effect=NA) )
datd <- dat0; datd$y <- x$St.a0m0; datd$effect <- "Control"
datm <- dat0; datm$y <- x$St.a1m0; datm$effect <- "Counterfactual"
datt <- dat0; datt$y <- x$St.a1m1; datt$effect <- "Treated"
dat <- rbind( datd, datm, datt )

pp <- ggplot(data=dat, aes(x=time, y=y, group = effect, color = effect)) + 
  geom_hline(yintercept = 0, color="gray") +
  geom_line( linewidth=1.5 ) + ggplot2::theme_minimal() + 
  geom_segment(aes(x=0, y=0.5, yend=0.5, xend=attributes(x)$summary[["median.ctr"]]), col="gray") + 
  geom_segment(aes(x=0, y=0.5, yend=0.5, xend=attributes(x)$summary[["median.trt"]]), col="gray") + 
  geom_segment(aes(x=attributes(x)$summary[["median.ctr"]], y=0, yend=0.5, xend=attributes(x)$summary[["median.ctr"]]), col="gray") +  
  geom_segment(aes(x=attributes(x)$summary[["median.trt"]], y=0, yend=0.5, xend=attributes(x)$summary[["median.trt"]]), col="gray") +   
  theme(
    plot.title = element_text(color = "black", face = "bold", size = 12),
    axis.title.x = element_text(color = "black", face = "bold", size = 12),
    axis.title.y = element_text(color = "black", face = "bold", size = 12),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(color = "black", face = "bold", size = 12),
    axis.text.y = element_text(color = "black", face = "bold", size = 12))
  if( length(attributes(x)$scenario)>0 ){     
   pp <- pp + labs(
    y = paste0("Survival", "\n"),
    x = paste0("\n","time"),
    title = "Estimated Survival functions",
    subtitle = paste0("Generated from true scenario = ",attributes(x)$scenario) )
  } else{
   pp <- pp + labs(
    y = paste0("Survival", "\n"),
    x = paste0("\n","time"),
    title = "Estimated Survival functions" )
  }    

 x1 <- stats::rbeta(10000,attributes(x)$beta.ctr[[1]],attributes(x)$beta.ctr[[2]])
 x2 <- stats::rbeta(10000,attributes(x)$beta.trt[[1]],attributes(x)$beta.trt[[2]])
 
 df <- as.data.frame( list(x=c(x1, x2), arm=factor(c(rep("control",length(x1)),rep("treated",length(x2)))) ))

 p <- ggplot2::ggplot(df, aes(color=arm, fill=arm) ) +  
      ggplot2::geom_density(mapping = ggplot2::aes(x = x, group = arm), linetype=0, alpha = 0.2) +
                             geom_vline(xintercept = mean(x1), color="darkgray") +
                             geom_vline(xintercept = mean(x2), color="darkgray") +
                             ggplot2::theme_minimal() +
                             ggplot2::theme(
                               panel.grid.major.y = ggplot2::element_blank(),
                               panel.grid.minor.y = ggplot2::element_blank(),
                               axis.ticks.y = ggplot2::element_blank(),
                               axis.text.y=element_blank(),
                               plot.title = ggplot2::element_text( color = "black", face = "bold", size = 12),
                               axis.title.x = ggplot2::element_text( color = "black", face = "bold", size = 12),
                               axis.text.x = ggplot2::element_text( color = "black", face = "bold", size = 12)
                             ) +
                             scale_fill_manual(values = c("#E15759","#4E79A7") ) +
                             scale_color_manual(values = c("#E15759","#4E79A7") ) +
                             ggplot2::labs(
                               x = paste0("\n", "response rate"),
                               y = "",
                               title = paste0("Surrogate response rate posterior distributions")
                              )

 out <- pp / p  + plot_layout(heights = c(6, 1))

 return(out)
}


#' Plot Survival end of study posterior predictive trial summaries

#' @param x object of class pp_eos
#' @param threshold list of thresholds
#' @importFrom viridis inferno turbo viridis
#' @importFrom stats median quantile
#' @import ggplot2
#' @import patchwork
#' @return ggplot2 object
#' @export

plot_pp_eos <- function(x, threshold=list(p=0.05, hr=0.8) ){

 if( !inherits(x,"pp_eos") ){
  stop("x must be of class pp_eos")
 }

 eos <- attributes(x)[["eos"]]
 eos_units <- attributes(x)[["tunits"]]

 col_pal <- viridis::inferno(128)[1:100]

 pp <- sum( x["p",] < threshold$p )/ncol(x)

 col_p <- col_pal[[ max(round( 100*pp ),1) ]]

 df <- as.data.frame( list(x=x["p",]), stringsAsFactors = TRUE )

 if( is.null(eos) || is.null(eos_units) ){
  #sub1 = paste0("prob of stat signif = ",round(pp,3))
  xlab1 <- paste0("\n", "log-rank test p-value")
 } else{
  #sub1 = paste0("prob of stat signif at ",eos," ",eos_units," = ",round(pp,3))
  xlab1 <- paste0("\n", "log-rank test p-value at ",eos," ",eos_units)
 }
 sub1 = paste0("prob of stat signif (",threshold$p,") = ",round(pp,3))

 p1 <- ggplot2::ggplot(df, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram( #ggplot2::aes(y = after_stat(density)),
        fill = col_p, position = "identity", bins = 50) +
      ggplot2::theme_minimal() +
      ggplot2::geom_vline(xintercept = threshold$p, color="maroon") +
      ggplot2::theme( panel.grid.major.y = ggplot2::element_blank(),
                      panel.grid.minor.y = ggplot2::element_blank(),
                      panel.grid.major.x = ggplot2::element_blank(),
                      panel.grid.minor.x = ggplot2::element_blank(),
                      axis.ticks.y = ggplot2::element_blank(),
                      axis.text.y=ggplot2::element_blank(),
                      plot.title = ggplot2::element_text( color = "black", face = "bold", size = 12),
                      axis.title.x = ggplot2::element_text( color = "black", face = "bold", size = 12),
                      axis.text.x = ggplot2::element_text( color = "black", face = "bold", size = 12) ) +
        ggplot2::labs( x = xlab1,
                       y = "",
                       #title = paste0("Predictive distribution: log-rank test"),
                       title = paste0("Log-rank test (after ",attributes(x)[["summary"]][[1]],")"),
                       subtitle = sub1 )


 pp_hr <- sum( x["hr",] < threshold$hr )/ncol(x)

 col_hr <- col_pal[[ max(round( 100*pp_hr ),1) ]]

 df <- as.data.frame( list(x=x["hr",]), stringsAsFactors = TRUE )

 if( is.null(eos) || is.null(eos_units) ){
  #sub2 = paste0("prob of clinical signif = ",round(pp_hr,3))
  xlab2 <- paste0("\n", "hazard ratio")
 } else{
  #sub2 = paste0("prob of clinical signif at ",eos," ",eos_units," = ",round(pp_hr,3))
  xlab2 <- paste0("\n", "hazard ratio at ",eos," ",eos_units)
 }
 sub2 = paste0("prob clinical signif (",threshold$hr,") = ",round(pp_hr,3))

 p2 <- ggplot2::ggplot(df, ggplot2::aes(x = x)) +
      ggplot2::geom_histogram( #ggplot2::aes(y = ggplot2::after_stat(density)),
       fill = col_hr, position = "identity", bins = 50) +
      ggplot2::theme_minimal() +
      ggplot2::geom_vline(xintercept = threshold$hr, color="maroon") +
      ggplot2::theme( panel.grid.major.y = ggplot2::element_blank(),
                      panel.grid.minor.y = ggplot2::element_blank(),
                      panel.grid.major.x = ggplot2::element_blank(),
                      panel.grid.minor.x = ggplot2::element_blank(),                      
                      axis.ticks.y = ggplot2::element_blank(),
                      axis.text.y=ggplot2::element_blank(),
                      plot.title = ggplot2::element_text( color = "black", face = "bold", size = 12),
                      axis.title.x = ggplot2::element_text( color = "black", face = "bold", size = 12),
                      axis.text.x = ggplot2::element_text( color = "black", face = "bold", size = 12) ) +
        ggplot2::labs( x = xlab2,
                       y = "",
                       #title = paste0("Predictive distribution: hazard ratio"),
                       title = paste0("Hazard ratio (after ",attributes(x)[["summary"]][[1]],")"),
                       subtitle = sub2 )


 x1 <- x[3,]
 x2 <- x[4,]
 m1 <- stats::median(x1,na.rm=TRUE)
 m2 <- stats::median(x2,na.rm=TRUE)
 if( m1>m2 ){
  q1 <- stats::quantile(x1,0.75,na.rm=TRUE)
  q2 <- stats::quantile(x2,0.25,na.rm=TRUE)
 } else{
  q1 <- stats::quantile(x1,0.25,na.rm=TRUE)
  q2 <- stats::quantile(x2,0.75,na.rm=TRUE)
 }

 tlabs <- gsub("trt\\=","",rownames(x)[3:4])

 if( is.null(eos) || is.null(eos_units) ){
  xlab3 <- paste0("\n", "Median survival")
 } else{
  xlab3 <- paste0("\n", "Median survival in ",eos_units)
 }

 col_mpal <- viridis::turbo(100)

 df <- as.data.frame( list(x=c(x1, x2), arm=factor(c(rep(tlabs[[1]],length(x1)),rep(tlabs[[2]],length(x2)))) ))

  p3 <- ggplot2::ggplot(data=df, ggplot2::aes(x = x, y = arm)) +
      ggplot2::geom_violin(ggplot2::aes(fill=arm), trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width")
  ggb <- ggplot2::ggplot_build(p3)

  cg <- seq(ggb$layout$panel_params[[1]]$x$continuous_range[[1]],
            ggb$layout$panel_params[[1]]$x$continuous_range[[2]], length.out=100)

  c1 <- col_mpal[100:1][[ which.min(abs(cg-q1)) ]]
  c2 <- col_mpal[100:1][[ which.min(abs(cg-q2)) ]]

  vcol <- c(c1,c2); names(vcol) <- tlabs

  p3 <- ggplot2::ggplot(data=df, ggplot2::aes(x = x, y = arm)) +
      ggplot2::geom_violin(ggplot2::aes(fill=arm), trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width") +
      ggplot2::scale_fill_manual(values = vcol, na.value = "#7F7F7F") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_text(color = "black", face = "bold", size = 12),
      axis.title.y = ggplot2::element_text(color = "black", face = "bold", size = 12),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "black", face = "bold", size = 10),
      axis.text.y = ggplot2::element_text(color = "black", face = "bold", size = 10)) +
      ggplot2::scale_x_continuous(n.breaks = 5) + 
      ggplot2::labs(
       y = paste0("study arm", "\n"),
       x = xlab3 )

 ggb <- ggplot2::ggplot_build(p3)
 
 p3 <- p3 +
      ggplot2::geom_segment( ggplot2::aes(y=ggb$layout$panel_params[[1]]$y$continuous_range[[1]], 
                                 x=m1, xend=m2, yend=ggb$layout$panel_params[[1]]$y$continuous_range[[1]]), size=1, color="maroon", show.legend=FALSE ) #+
      #ggplot2::geom_segment( ggplot2::aes(y=ggb$layout$panel_params[[1]]$y$continuous_range[[1]]-0.1, 
      #                           x=m1, xend=m1, yend=ggb$layout$panel_params[[1]]$y$continuous_range[[1]]), size=1, color=vcol[[1]], show.legend=FALSE ) +
      #ggplot2::geom_segment( ggplot2::aes(y=ggb$layout$panel_params[[1]]$y$continuous_range[[1]]-0.1, 
      #                           x=m2, xend=m2, yend=ggb$layout$panel_params[[1]]$y$continuous_range[[1]]), size=1, color=vcol[[2]], show.legend=FALSE )

 out <- ( p1 + p2 ) / p3 + plot_layout(widths = c(1, 1), heights = c(1.75, 1))

 return(out)
}


#' Plot end of study posterior predictive trial summaries for sequence of analysis times

#' @param dat data.frame of class mpti_weib_dat
#' @param df object of class pp_eos_df
#' @param alternative A data.frame
#' @importFrom viridis inferno viridis
#' @importFrom stats median
#' @importFrom scales pretty_breaks
#' @import ggplot2
#' @import patchwork
#' @return ggplot2 object
#' @export

plot_seq_pp_eos <- function(dat, df, alternative="two.sided"){

 plot_pp_p <- function(df,target,eos,tunits){

  df$type <- NULL
  df$arm <- NULL

  df <- rbind( df, c(eos,target) )
  df$fill <- factor(df$x)

  med <-  NULL
  for(ii in unique(df$x)){
   med <- c(med, stats::median(base::subset(df,df$x==ii)$y,na.rm=TRUE))
  }

  col_pal1 <- viridis::inferno(128)[101:1]
  col_pal2 <- viridis::viridis(128)[101:1]
  col_tar <- c(col_pal2[ round(100*med[1:(length(med)-1)])+1 ],col_pal1[ round(100*med[[length(med)]])+1 ] )
 
  p <- ggplot2::ggplot(data=df, ggplot2::aes(x = x, y = y), show.legend=FALSE) + 
       ggplot2::geom_hline(yintercept = 0.05, color="grey", size=0.5) +
       ggplot2::geom_violin(ggplot2::aes(fill=fill), trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width", show.legend=FALSE) +
       ggplot2::scale_fill_manual(values = col_tar, na.value = "#7F7F7F") +
       ggplot2::geom_point( ggplot2::aes(x=eos, y=med[[length(med)]]), size=8, 
                                                     col=col_tar[[length(col_tar)]],
                                                   shape=18 ) +    
       ggplot2::theme_minimal() +
       ggplot2::theme( 
        panel.grid.major.y = ggplot2::element_blank(),
                      panel.grid.minor.y = ggplot2::element_blank(),
                      panel.grid.major.x = ggplot2::element_blank(),
                      panel.grid.minor.x = ggplot2::element_blank(), 
        plot.title = ggplot2::element_text( color = "black", face = "bold", size = 12),
                      axis.title.x = ggplot2::element_text( color = "black", face = "bold", size = 11),
                      axis.title.y = ggplot2::element_text( color = "black", face = "bold", size = 11),
                      axis.text.x = ggplot2::element_text( color = "black", face = "bold", size = 10),
                      axis.text.y = ggplot2::element_text( color = "black", face = "bold", size = 10) ) +
        ggplot2::labs( x = paste0("\n","trial time (",tunits,")"),
                       y = paste0("post pred distr", "\n"),
                       #title = paste0("Predictive distribution: log-rank test"),
                       title = paste0("Log-rank test prediction") ) +
                       ggplot2::scale_x_continuous(breaks = unique(df$x)) +
                       ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10))        

  return(p)
 }

 plot_pp_hr <- function(df,target,eos,tunits){

  df$type <- NULL
  df$arm <- NULL

  df <- rbind( df, c(eos,target) )
  df$fill <- factor(df$x)

  med <-  NULL
  for(ii in unique(df$x)){
   med <- c(med, min( stats::median(base::subset(df,df$x==ii)$y,na.rm=TRUE), 1/stats::median(base::subset(df,df$x==ii)$y,na.rm=TRUE)) )
  }

  col_pal1 <- viridis::inferno(128)[1:101]
  col_pal2 <- viridis::viridis(128)[101:1]
  col_tar <- c(col_pal2[ round(100*med[1:(length(med)-1)])+1 ],col_pal1[ round(100*med[[length(med)]])+1 ] )
 
  p <- ggplot2::ggplot(data=df, ggplot2::aes(x = x, y = y), show.legend=FALSE) +
       ggplot2::geom_hline(yintercept = 1, color="grey", size=0.5) +
       ggplot2::geom_violin(ggplot2::aes(fill=fill), trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width", show.legend=FALSE) +
       ggplot2::scale_fill_manual(values = col_tar, na.value = "#7F7F7F") +
       ggplot2::geom_point( ggplot2::aes(x=eos, y=med[[length(med)]]), size=8, 
                                                     col=col_tar[[length(col_tar)]],
                                                   shape=18 ) +    
       ggplot2::theme_minimal() +
       ggplot2::theme( 
        panel.grid.major.y = ggplot2::element_blank(),
                      panel.grid.minor.y = ggplot2::element_blank(),
                      panel.grid.major.x = ggplot2::element_blank(),
                      panel.grid.minor.x = ggplot2::element_blank(), 
        plot.title = ggplot2::element_text( color = "black", face = "bold", size = 12),
                      axis.title.x = ggplot2::element_text( color = "black", face = "bold", size = 11),
                      axis.title.y = ggplot2::element_text( color = "black", face = "bold", size = 11),
                      axis.text.x = ggplot2::element_text( color = "black", face = "bold", size = 10),
                      axis.text.y = ggplot2::element_text( color = "black", face = "bold", size = 10) ) +
        ggplot2::labs( x = paste0("\n","trial time (",tunits,")"),
                       y = paste0("post pred distr", "\n"),
                       #title = paste0("Predictive distribution: log-rank test"),
                       title = paste0("Hazard Ratio prediction") ) +
                       ggplot2::scale_x_continuous(breaks = unique(df$x)) +
                       ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10))        
  return(p)
 }

 plot_pp_med <- function(df,target,eos,tunits){

  df$type <- NULL

  df <- rbind( df, c(eos,target[[1]],unique(df$arm)[[1]]) )
  df <- rbind( df, c(eos,target[[2]],unique(df$arm)[[2]]) )

  ord <- NULL
  for(ii in unique(df$x)){
   ord <- c(ord, paste0(ii,"-",unique(df$arm)[[1]]))
   ord <- c(ord, paste0(ii,"-",unique(df$arm)[[2]]))
  }

  df$x2 <- paste0(df$x,"-",df$arm)
  df$y <- as.numeric(df$y)
  df$x2 <- factor(df$x2, levels=ord)

  med1 <- med2 <- NULL
  for(ii in base::setdiff(unique(df$x),eos)){
   temp <- base::subset(df,df$x==ii)
   temp1 <- base::subset(temp, temp$arm==unique(temp$arm)[[1]])$y
   temp2 <- base::subset(temp, temp$arm==unique(temp$arm)[[2]])$y
   med1 <- c(med1, stats::median(temp1,na.rm=TRUE) )
   med2 <- c(med2, stats::median(temp2,na.rm=TRUE) )
  }
  med1 <- c(med1, target[[1]])
  med2 <- c(med2, target[[2]])

  p3 <- ggplot2::ggplot(data=df, ggplot2::aes(x = x2, y = y)) +
        ggplot2::geom_violin(, trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width")
  suppressWarnings( ggb <- ggplot2::ggplot_build(p3) )

  cg <- seq(ggb$layout$panel_params[[1]]$x$continuous_range[[1]],
            ggb$layout$panel_params[[1]]$x$continuous_range[[2]], length.out=100)

  col_pal1 <- viridis::inferno(128)[1:100]
  col_pal2 <- viridis::viridis(128)[100:1]

  c1 <- c2 <- NULL
  for(ii in med1[-length(med1)]){
   c1 <- c(c1, col_pal2[ which.min(abs(cg-ii)) ] )
  }
  for(ii in med1[-length(med2)]){
   c2 <- c(c2, col_pal2[ which.min(abs(cg-ii)) ] )
  }
  c1 <- c(c1, col_pal1[ which.min(abs(cg-med1[[length(med1)]])) ] )
  c2 <- c(c2, col_pal1[ which.min(abs(cg-med2[[length(med2)]])) ] )

  col_tar <- NULL
  for(ii in 1:length(c1)){
   col_tar <- c(col_tar, c1[[ii]])
   col_tar <- c(col_tar, c2[[ii]])
  }

   p <- ggplot2::ggplot(data=df, ggplot2::aes(x = x2, y = y), show.legend=FALSE) + 
        ggplot2::geom_violin(ggplot2::aes(fill=x2), trim = TRUE, draw_quantiles = c(0.25, 0.5, 0.75), scale = "width", show.legend=FALSE) +
        ggplot2::scale_fill_manual(values = col_tar, na.value = "#7F7F7F") +
        ggplot2::geom_point( ggplot2::aes(x=paste0(eos,"-",unique(df$arm)[[1]]), y=target[[1]]), col=c1[[length(c1)]], size=8, shape=18, show.legend=FALSE ) +
        ggplot2::geom_point( ggplot2::aes(x=paste0(eos,"-",unique(df$arm)[[2]]), y=target[[2]]), col=c2[[length(c2)]], size=8, shape=18, show.legend=FALSE ) +                                                   
        ggplot2::theme_minimal() +
        ggplot2::theme( 
         panel.grid.major.y = ggplot2::element_blank(),
                       panel.grid.minor.y = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       panel.grid.minor.x = ggplot2::element_blank(), 
         plot.title = ggplot2::element_text( color = "black", face = "bold", size = 12),
                      axis.title.x = ggplot2::element_text( color = "black", face = "bold", size = 11),
                       axis.title.y = ggplot2::element_text( color = "black", face = "bold", size = 11),
                       axis.text.x = ggplot2::element_text( color = "black", face = "bold", size = 10),
                       axis.text.y = ggplot2::element_text( color = "black", face = "bold", size = 10) ) +
        ggplot2::labs( x = paste0("\n","trial time (",tunits,")"),
                        y = paste0("post pred distr", "\n"),
                        #title = paste0("Predictive distribution: log-rank test"),
                        title = paste0("Median survival prediction") ) +
                        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10))        

  return(p)
 }

 #############################################################################

 if( !inherits(df,"pp_eos_df") ){
  stop("x must be of class pp_eos_df")
 }

 eos <- attributes(dat)[["eos"]]

 targets <- c( attributes(dat)[["eos_stats"]]$OS[[alternative]],
              attributes(dat)[["eos_stats"]]$OS[["hr"]],
              attributes(dat)[["eos_stats"]]$OS[["medians"]][[1]],attributes(dat)[["eos_stats"]]$OS[["medians"]][[2]])

 p1 <- plot_pp_p(df=base::subset(df,df$type=="p"),target=targets[[1]],eos=eos,tunits=attributes(dat)[["tunits"]])
 p2 <- plot_pp_hr(df=base::subset(df,df$type=="hr"),target=targets[[2]],eos=eos,tunits=attributes(dat)[["tunits"]])
 p3 <- plot_pp_med(df=base::subset(df,df$type=="med"),target=targets[3:4],eos=eos,tunits=attributes(dat)[["tunits"]])

 out <- ( p1 + p2 ) / p3 + plot_layout(widths = c(1, 1), heights = c(1.75, 1))

 return(out)
}

