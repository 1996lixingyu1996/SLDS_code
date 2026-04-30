#### Expe 6 Global####
library(WeibullMultiState)
library(tidyverse)
#library(randomForestSRC)
library(xgboost)
library(MLearner)

f_exp1 <- function(seed){
  
  options(rf.cores=1, mc.cores = 1)
  prtrt = 0.5
  
  parm.tmp <- parm_mpti_weib(
    ## intercept and coefs for trt, x1, x2 in response model
    gamma = c(log(0.15/0.85), ### intercept (log-odds of response)
              0, ### coefficient for log-odds of treatment effect
              rep(0,30) ### coefficient for prognostic covariate x1
    ), ### coefficient for prognostic covariate x2
    
    ### Parameters controlling baseline TTF models for PFS and OS ####
    ## Weibull shape (nv) and rate (lamb) parameters for PFS (12), OS (13), and OS (23) after progression
    #weib_nv = list(nv12=3, nv13=2, nv23=4),
    weib_nv = list(nv12=1.88, nv13=1.88, nv23=1.88),
    weib_lamb = list(lamb12=(1/55), lamb13=(1/300), lamb23=(1/300)),
    
    ### Coefficients for treatment, response, and prognostic factors (x) in each multistate model ###
    ## beta12 -> treatment (direct effect), response, x for hazard of PFS
    ## beta13 -> treatment (direct effect), response, x for hazard of OS
    ## beta23 -> treatment (direct effect), response, x for hazard of OS after progression
    
    beta = list( beta12=c(0, 0, rep(0,31)),  # no care
                 beta13=c(0, 0, rep(0,31)),  # first value, -1/1
                 beta23=c(0, 0, rep(0,31)) ),
    
    prob = 0.5,  ## binomial probability vector for binary prognostic covariates
    mu = 0, ## mean vector for numeric prognostic covariates
    sigma = 1, ## stand dev vector for numeric prognostic covariates
    weib_shape = 5, ## Weibull shape paramter for time-to-observing surrogate response
    weib_scale = 1.5, ## Weibull scale paramter for time-to-observing surrogate response
    scenario = "scenario_name" ## text string to label the scenario
  )
  
  parm.tmp$gamma = c(log(0.15/0.85), ### intercept (log-odds of response)
                     0, ### coefficient for log-odds of treatment effect ### coefficient for prognostic covariate x1
                     rep(0,10), 0, 0) # 
  
  parm.tmp$beta$beta12 = c(0, 0, rep(0,10), 0, 0)  # 
  parm.tmp$beta$beta13 = c(0, 0, c(-0.1,-0.1,0.1,0.1), rep(0,6), -0.75, -0)  # 
  parm.tmp$beta$beta23 = c(0, 0, c(-0.1,-0.1,0.1,0.1), rep(0,6), -0.75, -0)  # 
  parm.tmp$prob = rep(0.5,0)
  parm.tmp$mu = rep(0,12)
  parm.tmp$sigma = rep(1,12)
  
  n = 1000
  
  parm = parm.tmp
  
  var <- gen_x(n=n, pb=0, pc=12, prob=parm$prob, mu=parm$mu, sigma=parm$sigma, seed = seed)  # matrix is covariate
  var = as.data.frame(var)
  colnames(var) = c(paste0("Continuous", 1:10),"mediator","tau")
  
  var$trt = gen_a(n=n,p=prtrt,seed = seed)
  
  ## var13 is TRT
  ## var11 is mediator
  ## var12 is tau(x,w)
  for (i in 1:n){
    error1 = rnorm(n=1,mean=0,sd=0.01)
    #error2 = rnorm(n=1,mean=0,sd=0.01)
    if(var[i,13]==1){
      var[i,11] = 0.1*var[i,1] + 0.1*var[i,2] +1+error1
      var[i,12] = 1
    }else{
      var[i,11] = 0.1*var[i,1] + 0.1*var[i,2] +error1
      var[i,12] = 0
    }
  }
  
  ####### heterogeneous variable indicator
  
  a = var$trt
  x = as.matrix(var[,c(1:12)])
  n_continuous_var = 12
  
  nb = 0   #n_binary_var + n_heterogeneous_var
  nc = 12   #n_continuous_var
  
  dat = gen_heterogeneous_data(parm = parm, n=n, a = a, x = x, nb = nb, nc = nc, prtrt=0.5,accTime=40,accExN=NULL,addFUt=0,CenUpLim=NULL,tunits="months",eos_stats=TRUE)
  
  dat <- dat %>% dplyr::rename(TRT=trt)
  
  fit = survival::coxph(Surv(OS,OSevent)~TRT, data = dat[dat$Continuous1>0&dat$Continuous2>0,])
  
  fit.all = survival::coxph(Surv(OS,OSevent)~TRT, data = dat)
  
  HR.all = exp(as.numeric(fit.all$coefficients))
  HR.subtype = exp(as.numeric(fit$coefficients))
  
  dat_trt1 = dat[dat$TRT==1,]
  dat_trt0 = dat[dat$TRT==0,]
  
  #fit.tmp = survival::coxph(Surv(OS,OSevent)~Continuous1+Continuous2+Continuous3+Continuous4+Continuous5+Continuous6+Continuous7+Continuous8+Continuous9+Continuous10, data = dat)
  
  dat_save_path = paste0(".../expe_null/data_seed_",seed,".csv")
  rio::export(dat, dat_save_path)
  
  p_min_r = 1
  result_flag = NULL
  for (r in 1:7){
    
    x_all <- model.matrix(~ Continuous1+Continuous2+Continuous3+Continuous4+Continuous5+Continuous6+Continuous7+Continuous8+Continuous9+Continuous10 - 1, data = dat)
    y_time_all <- dat$OS
    y_event_all <- dat$OSevent
    m_all = dat$mediator
    
    ## treatment data
    x_1 <- model.matrix(~ Continuous1+Continuous2+Continuous3+Continuous4+Continuous5+Continuous6+Continuous7+Continuous8+Continuous9+Continuous10 - 1, data = dat_trt1)
    y_time_1 <- dat_trt1$OS
    y_event_1 <- dat_trt1$OSevent
    m_1 <- dat_trt1$mediator
    
    x_0 <- model.matrix(~ Continuous1+Continuous2+Continuous3+Continuous4+Continuous5+Continuous6+Continuous7+Continuous8+Continuous9+Continuous10 - 1, data = dat_trt0)
    y_time_0 <- dat_trt0$OS
    y_event_0 <- dat_trt0$OSevent
    m_0 <- dat_trt0$mediator
    
    xm_1 <- model.matrix(~ mediator+Continuous1+Continuous2+Continuous3+Continuous4+Continuous5+Continuous6+Continuous7+Continuous8+Continuous9+Continuous10 - 1, data = dat_trt1)
    
    ## Estimate m by arms firstly
    dtrain_m_arm1 <- xgb.DMatrix(data = x_1, label = m_1)
    dtrain_m_arm0 <- xgb.DMatrix(data = x_0, label = m_0)
    
    dtrain_m_all = xgb.DMatrix(data = x_all, label = m_all)
    
    params <- list(
      objective = "reg:squarederror",  # Regression
      eta = 0.08,                        # Learning rate
      max_depth = 3,
      lambda=1,
      subsample = 0.7
    )
    
    xgb_model1 = xgb.train(
      params = params,
      data = dtrain_m_arm1,
      nrounds = 100
    )
    
    xgb_model0 = xgb.train(
      params = params,
      data = dtrain_m_arm0,
      nrounds = 100
    )
    
    m1_pred <- predict(xgb_model1, dtrain_m_all)
    m0_pred <- predict(xgb_model0, dtrain_m_all)
    label = ifelse(dat_trt1$OSevent==1,dat_trt1$OS,-dat_trt1$OS)
    dtrain <- xgb.DMatrix(data = xm_1, label = label)
    
    setinfo(dtrain, "label_lower_bound", ifelse(dat_trt1$OSevent == 1, dat_trt1$OS, NA))
    setinfo(dtrain, "label_upper_bound", ifelse(dat_trt1$OSevent == 1, dat_trt1$OS, NA))
    param <- list(
      objective = "survival:cox",
      eval_metric = "cox-nloglik",
      eta = 0.06,
      max_depth = 3
    )
    
    xgb_cox <- xgb.train(
      params = param,
      data = dtrain,
      nrounds = 80,
      subsample=0.67
    )
    
    data_arm1_ct = dat 
    data_arm1_ct$mediator = m1_pred
    xm_tmp <- model.matrix(~ mediator+Continuous1+Continuous2+Continuous3+Continuous4+Continuous5+Continuous6+Continuous7+Continuous8+Continuous9+Continuous10 - 1, data = data_arm1_ct)
    dtrain_ct1 = xgb.DMatrix(data = xm_tmp)
    risk_ct1 = predict(xgb_cox, dtrain_ct1)
    
    data_arm0_ct = dat 
    data_arm0_ct$m = m0_pred
    xm_tmp <- model.matrix(~ mediator+Continuous1+Continuous2+Continuous3+Continuous4+Continuous5+Continuous6+Continuous7+Continuous8+Continuous9+Continuous10 - 1, data = data_arm0_ct)
    dtrain_ct0 = xgb.DMatrix(data = xm_tmp)
    risk_ct0 = predict(xgb_cox, dtrain_ct0)
    
    indirect = (log(risk_ct1) - log(risk_ct0))
    
    #### risk  score
    dat$indirect_effect = indirect  #indir_trt_effect
    dis3 = as.data.frame(matrix(rep(0,n*n),nrow = n))
    for (i in 1:n){
      dis3[i,] = rep(as.numeric(dat[i,"indirect_effect"]),n) - as.numeric(dat$indirect_effect)
    }
    dis3 = (dis3)**2
    tsne_result = Rtsne::Rtsne(dis3,dims=2,is_distance=TRUE,verbose=FALSE,max_iter = 3000, theta = 0)
    
    kmeans_result_2 = kmeans(tsne_result$Y, centers = 2, iter.max = 50, nstart = 10)
    kmeans_result_3 = kmeans(tsne_result$Y, centers = 3, iter.max = 50, nstart = 10)
    kmeans_result_4 = kmeans(tsne_result$Y, centers = 4, iter.max = 50, nstart = 10)
    kmeans_result_5 = kmeans(tsne_result$Y, centers = 5, iter.max = 50, nstart = 10)
    
    dat$kmeans2 = as.factor(kmeans_result_2$cluster)
    dat$kmeans3 = as.factor(kmeans_result_3$cluster)
    dat$kmeans4 = as.factor(kmeans_result_4$cluster)
    dat$kmeans5 = as.factor(kmeans_result_5$cluster)
    
    pval_flag_rs = 1
    data_flag = NULL
    result = list()
    for (j in 2:5){
      
      x_name = paste0("Continuous",1:10)
      c_name = paste0("kmeans",j)  
      profiles <- tree_fit(Y=dat[,c_name], X=dat[,x_name], seed=1234,min_leaf = 200)
      
      
      if(length(profiles$trees)>0){
        for (k in 1:length(profiles$trees)){
          
          if(has_repeated_vars(profiles$trees[[k]])){
            next
          }
          
          dat = predict_path(profiles$trees[[k]], newdata = dat)
          dat$leaf_rs = dat$leaf
          profiles$trees[[k]]
          data_tmp = dat
          fit_lm1 = lm(mediator~leaf+TRT+leaf*TRT, data = data_tmp)
          fit_lm0 = lm(mediator~TRT, data = data_tmp)
          pval_M_rs <- as.numeric(na.omit(stats::anova(fit_lm0, fit_lm1)[[6]]))
          
          fit_surv1 = coxph(Surv(OS,OSevent)~mediator+leaf, data = data_tmp)
          fit_surv0 = coxph(Surv(OS,OSevent)~mediator, data = data_tmp)
          pval_Y_rs = as.numeric(na.omit(stats::anova(fit_surv0, fit_surv1)[[4]]))
          
          fit_surv3 = coxph(Surv(OS,OSevent)~leaf+TRT+leaf*TRT, data = data_tmp)
          fit_surv2 = coxph(Surv(OS,OSevent)~TRT, data = data_tmp)
          pval_Y1_rs = as.numeric(na.omit(stats::anova(fit_surv2, fit_surv3)[[4]]))
          
          if(pval_M_rs < pval_flag_rs){
            pval_flag_rs = pval_M_rs
            result[[1]] = profiles
            result[[2]] = pval_flag_rs
            result[[3]] = pval_Y_rs
            result[[4]] = pval_Y1_rs
            result[[5]] = k
            result[[8]] = 1
          }
        }
      }
      
      
    }
    
    result[[6]] = HR.all
    result[[7]] = HR.subtype
    result[[8]] = r
    
    if(is.null(result[[2]])){
      next
    }
    if (result[[2]]<=p_min_r){
      p_min_r = result[[2]]
      result_flag = result
    }
  }
  
  save_path = paste0(".../expe_null/seed_",seed,".rds")
  
  rio::export(result_flag, save_path)
}


library(Rtsne)
library(rio)
library(snowfall)
library(xgboost)
library(survival)
library(parallel)

sfInit(parallel = TRUE, cpus = (detectCores() - 2))
sfLibrary(dplyr)
sfLibrary(rio)
sfLibrary(rpart)
sfLibrary(xgboost)
sfLibrary(base)
sfLibrary(WeibullMultiState)
sfLibrary(Rtsne)
sfLibrary(MLearner)
sfLibrary(survival)

my_new_folder = ".../expe_null"
if (!dir.exists(my_new_folder)) {
  dir.create(my_new_folder)
}

sfSource(".../decision_tree.R")

result_vector = sfLapply(1:100, f_exp1)

sfStop()




extract_rpart_thresholds <- function(tree_model) {
  # 
  if (!inherits(tree_model, "rpart")) {
    stop("must rpart!")
  }
  
  # 
  frame <- tree_model$frame
  splits <- tree_model$splits
  var_names <- as.character(frame$var)
  
  # 
  node_ids <- as.numeric(rownames(frame))
  internal_nodes <- which(var_names != "<leaf>")
  
  # 
  thresholds <- data.frame(
    node = node_ids[internal_nodes],
    variable = var_names[internal_nodes],
    threshold = NA_real_
  )
  
  # ）
  split_vars <- rownames(splits)
  split_vals <- splits[, "index"]
  
  # 
  for (i in seq_along(thresholds$variable)) {
    var <- thresholds$variable[i]
    
    # ）
    match_idx <- which(split_vars == var)[1]
    if (!is.na(match_idx)) {
      thresholds$threshold[i] <- split_vals[match_idx]
      split_vars[match_idx] <- NA  # 
    }
  }
  
  return(thresholds)
}


density_weighted_error <- function(a, b, mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1) {
  f_a <- dnorm(a, mean = mu1, sd = sigma1)
  f_b <- dnorm(b, mean = mu2, sd = sigma2)
  
  error1 <- abs(a - mu1) / (1 + f_a)
  error2 <- abs(b - mu2) / (1 + f_b)
  
  total_error <- error1^2 + error2^2
  return(total_error)
}

folder_path = ".../expe_null"

true_n_vec = NULL
sim_n_vec = NULL
p_med_vec = NULL
prop_med_vec = NULL
var1_list = NULL
var2_list = NULL
var12_list = NULL
n_var_list = NULL
weighted_error_list = NULL
trt_coef_list = NULL
a_list = NULL
b_list = NULL
true_trt_coef_list = NULL
p_1_list = NULL
p_2_list = NULL
p_3_list = NULL
p_min_list = NULL
p_profile_method = NULL
p_1_Y_list = NULL
p_2_Y_list = NULL
p_3_Y_list = NULL
index_list = NULL
p_Y_list = NULL
for (i in 1:100){
  print(i)
  data = rio::import(paste0(folder_path,"/data_seed_", i,".csv"))
  
  result = rio::import(paste0(folder_path,"/seed_", i,".rds"))
  
  p_1 = result[[2]] #+ result[[4]]
  p_2 = result[[3]] #+ result[[8]]
  p_3 = result[[4]] #+ result[[12]]
  #min_p = min(p_1,p_2,p_3)
  p_1_list = c(p_1_list, p_1)
  #p_1_Y_list = c(p_1_Y_list, result[[4]])
  p_2_list = c(p_2_list, p_2)
  #p_2_Y_list = c(p_2_Y_list, result[[9]])
  p_3_list = c(p_3_list, p_3)
  #p_3_Y_list = c(p_3_Y_list, result[[14]])
  
  #p_min_list = c(p_min_list, min_p)
  
  k = result[[5]]
  profile = result[[1]]$trees[[k]]
  var = names(result[[1]]$trees[[k]]$variable.importance)
  index_list = c(index_list, k)
  p_Y_list = c(p_Y_list, result[[4]])
  #p_profile_method = c(p_profile_method, 1)
  
  
  data_new = predict_path(profile, newdata = data)
  
  if("Continuous1"%in%var){
    var1_list = c(var1_list, i)
  }else{
    var1_list = c(var1_list, -1)
  }
  if("Continuous2"%in%var){
    var2_list = c(var2_list, i)
  }else{
    var2_list = c(var2_list, -1)
  }
  if("Continuous1"%in%var & "Continuous2"%in%var){
    var12_list = c(var12_list, i)
  }else{
    var12_list = c(var12_list, -1)
  }
  n_var_list = c(n_var_list, length(var))
  
  tmp = data[data$Continuous1>0&data$Continuous2>0,]
  sim_n = dim(tmp)[1]
  true_n_vec = c(true_n_vec, sim_n)
  
  fit0 = coxph(Surv(OS,OSevent)~TRT, data = tmp)
  true_coef = as.numeric(fit0$coefficients[1])
  true_trt_coef_list = c(true_trt_coef_list,true_coef)
  #dim(data_new[data_new$subtype==3,])
  data_new$leaf = as.factor(data_new$leaf)
  
  n_level_leaf = length(levels(data_new$leaf))
  level_leaf = levels(data_new$leaf)
  
  leaf_flag = 0
  p_flag = 1
  trt_flag = 100
  for (j in 1:n_level_leaf){
    data_leaf = data_new[data_new$leaf==level_leaf[j],]
    
    fit0 = coxph(Surv(OS,OSevent)~TRT, data = data_leaf)
    sum_fit=summary(fit0)
    trt_coef = sum_fit$coefficients[1,2]
    trt_coef
    dim(data_leaf)[1]
    p_med_prop = -1 #NULL#sum_out$n1.p
    
    if (trt_flag > trt_coef){
      trt_flag = trt_coef
      p_flag = p_med_prop
      n_subtype = dim(data_leaf)[1]
      med_prop = -1   #NULL#sum_out$n0
      leaf_flag = level_leaf[j]
    }
    
  }
  
  data_new$leaf_new = ifelse(data_new$leaf==leaf_flag, yes = 1, no = 0)
  data_new$leaf_new = as.factor(data_new$leaf_new)
  x_name = paste0("Continuous",1:10) #colnames(data_tmp)[c(1,3:6,9:20)] #paste0("Cov",1:100)
  c_name = "leaf_new"#paste0("kmeans",j)
  
  profiles <- tree_fit(Y=data_new[,c_name], X=data_new[,x_name], seed=1234)
  if(is.null(profiles)){
    #print(i)
    next
  }
  sim_n_vec = c(sim_n_vec, n_subtype)
  p_med_vec = c(p_med_vec, p_flag)
  trt_coef_list = c(trt_coef_list, trt_flag)
  prop_med_vec = c(prop_med_vec, med_prop)
  thre_df = extract_rpart_thresholds(profiles$trees[[1]])
  
  if(!dim(thre_df[thre_df$variable=="Continuous1",])[1]==0){
    flag = min(abs(thre_df[thre_df$variable=="Continuous1","threshold"]))
    tmp_df = thre_df[thre_df$variable=="Continuous1",]
    if(dim(tmp_df)[1] == 1){
      a = tmp_df[1,"threshold"]
    }else{
      a = tmp_df[abs(tmp_df$threshold) == flag, "threshold"]
    }
  }else{
    a = 1
  }
  #print(a)
  a_list = c(a_list,a)
  if(!dim(thre_df[thre_df$variable=="Continuous2",])[1]==0){
    flag = min(abs(thre_df[thre_df$variable=="Continuous2","threshold"]))
    tmp_df = thre_df[thre_df$variable=="Continuous2",]
    if(dim(tmp_df)[1] == 1){
      b = tmp_df[1,"threshold"]
    }else{
      b = tmp_df[abs(tmp_df$threshold) == flag, "threshold"]
    } 
  }else{
    b = 1
  }
  #print(b)
  b_list = c(b_list,b)
  
  error_flag = density_weighted_error(a=a, b=b, mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1)
  weighted_error_list = c(weighted_error_list, error_flag)
  #plot_profile(profiles$best.tree)
  
}

mean(weighted_error_list)

df = data.frame(error_score = weighted_error_list, X1=var1_list,X2=var2_list,
                X12=var12_list, n_var = n_var_list,sim_n = sim_n_vec,p_med = p_med_vec,
                true_n=true_n_vec,trt_coef = log(trt_coef_list),prop = prop_med_vec,
                thre_x1 = a_list,thre_x2 = b_list, true_trt_coef = true_trt_coef_list,
                p_1_list = p_1_list,p_2_list = p_2_list,p_3_list=p_3_list,
                index_list=index_list)

df_save_path = ".../expe_null.csv"
rio::export(df,df_save_path)





















