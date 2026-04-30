## instructions
## please ensure the xgboost package version is 1.7.9.1
## -----------------------------------------------------------------------------
library(MLearner)
## If MLearner package is not installed, please run code (assume current is code folder) source("./decision_tree.R")


library(survival)
library(BART)
data(ACTG175)
library(survival)
library(rpart)
library(rpart.plot)
library(xgboost)
library(survival)
library(survminer)
library(ggplot2)
set.seed(1)





## -----------------------------------------------------------------------------
seed = 1
set.seed(seed)
data = ACTG175
data$preanti_binary = ifelse(data$preanti>0,yes = "yes",no = "no")
data = data[data$arms%in%c(0,1),]
data$hemo = as.factor(data$hemo)
data$homo = as.factor(data$homo)
data$drugs = as.factor(data$drugs)
data$gender = as.factor(data$gender)
data$oprior = as.factor(data$oprior)
data$z30 = as.factor(data$z30)
data$symptom = as.factor(data$symptom)
data$m = data$cd420
data$M = as.factor(ifelse(data$m>=200, yes = 1, no = 0))


data0 = data[data$arms==0,]
data1 = data[data$arms==1,]


## -----------------------------------------------------------------------------

p_min_r = 1
result_flag = NULL
data_flag = NULL

for (r in 1:5){
  #print("r")
  cat(r)
  set.seed(r)
  x_all <- model.matrix(~ gender+age+wtkg+hemo+homo+drugs+karnof+race+preanti+symptom+cd40 - 1, data = data)
  y_time_all <- data$days
  y_event_all <- data$cens
  m_all = data$m

  ## treatment data
  x_1 <- model.matrix(~ gender+age+wtkg+hemo+homo+drugs+karnof+race+preanti+symptom+cd40 - 1, data = data1)
  y_time_1 <- data1$days
  y_event_1 <- data1$cens
  m_1 <- data1$m

  x_0 <- model.matrix(~ gender+age+wtkg+hemo+homo+drugs+karnof+race+preanti+symptom+cd40 - 1, data = data0)
  y_time_0 <- data0$days
  y_event_0 <- data0$cens
  m_0 <- data0$m

  xm_1 <- model.matrix(~ m+gender+age+wtkg+hemo+homo+drugs+karnof+race+preanti+symptom+cd40 - 1, data = data1)


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


  dtrain <- xgb.DMatrix(data = xm_1, label = data1$days)


  setinfo(dtrain, "label_lower_bound", ifelse(data1$cens == 1, data1$days, NA))
  setinfo(dtrain, "label_upper_bound", ifelse(data1$cens == 1, data1$days, NA))
  param <- list(
    objective = "survival:cox",
    eval_metric = "cox-nloglik",
    eta = 0.05,
    max_depth = 3
  )

  xgb_cox <- xgb.train(
    params = param,
    data = dtrain,
    nrounds = 50,
    subsample=0.67
  )

  data_arm1_ct = data
  data_arm1_ct$m = m1_pred
  xm_tmp <- model.matrix(~ m+gender+age+wtkg+hemo+homo+drugs+karnof+race+preanti+symptom+cd40 - 1, data = data_arm1_ct)
  dtrain_ct1 = xgb.DMatrix(data = xm_tmp)
  risk_ct1 = predict(xgb_cox, dtrain_ct1)

  data_arm0_ct = data
  data_arm0_ct$m = m0_pred
  xm_tmp <- model.matrix(~ m+gender+age+wtkg+hemo+homo+drugs+karnof+race+preanti+symptom+cd40 - 1, data = data_arm0_ct)
  dtrain_ct0 = xgb.DMatrix(data = xm_tmp)
  risk_ct0 = predict(xgb_cox, dtrain_ct0)

  indirect = (log(risk_ct1) - log(risk_ct0))

  data$indir_arms_effect = indirect
  n=dim(data)[1]
  dis3 = as.data.frame(matrix(rep(0,n*n),nrow = n))
  for (i in 1:n){
    dis3[i,] = rep(as.numeric(data[i,"indir_arms_effect"]),n) - as.numeric(data$indir_arms_effect)
  }
  dis3 = (dis3)**2
  tsne_result = Rtsne::Rtsne(dis3,dims=2,is_distance=TRUE,verbdayse=FALSE,max_iter = 2000, theta = 0)

  kmeans_result_2 = kmeans(tsne_result$Y, centers = 2, iter.max = 50, nstart = 10)
  kmeans_result_3 = kmeans(tsne_result$Y, centers = 3, iter.max = 50, nstart = 10)
  kmeans_result_4 = kmeans(tsne_result$Y, centers = 4, iter.max = 50, nstart = 10)
  kmeans_result_5 = kmeans(tsne_result$Y, centers = 5, iter.max = 50, nstart = 10)
  kmeans_result_6 = kmeans(tsne_result$Y, centers = 6, iter.max = 50, nstart = 10)

  data$kmeans2 = as.factor(kmeans_result_2$cluster)
  data$kmeans3 = as.factor(kmeans_result_3$cluster)
  data$kmeans4 = as.factor(kmeans_result_4$cluster)
  data$kmeans5 = as.factor(kmeans_result_5$cluster)
  data$kmeans6 = as.factor(kmeans_result_6$cluster)

  result = list()
  pval_flag_prob = 1


  for (j in 2:6){
    #print("j")
    #cat(j)
    x_name = c("gender","age","wtkg","hemo","homo","drugs","karnof", "preanti","race", "str2","symptom","cd40")
    c_name = paste0("kmeans",j)

    profiles <- tree_fit(Y=data[,c_name], X=data[,x_name], seed=1234,maxdepth = 2,min_leaf = 200)

    for (l in 1:length(profiles$trees)){

      if(has_repeated_vars(profiles$trees[[l]])){
        next
      }
      dat = predict_path(profiles$trees[[l]], newdata = data)

      if(!is.null(dat$leaf)){
        dat$leaf_rs = dat$leaf
        data_tmp = dat
        fit_lm1 = lm(cd420~leaf+arms+leaf*arms, data = data_tmp)
        fit_lm0 = lm(cd420~arms, data = data_tmp)
        pval_M_prob <- as.numeric(na.omit(stats::anova(fit_lm0, fit_lm1)[[6]]))

        fit_surv1 = coxph(Surv(days,cens)~cd420+leaf, data = data_tmp)
        fit_surv0 = coxph(Surv(days,cens)~cd420, data = data_tmp)
        pval_Y_prob = as.numeric(na.omit(stats::anova(fit_surv0, fit_surv1)[[4]]))

        fit_surv3 = coxph(Surv(days,cens)~leaf+arms+leaf*arms, data = data_tmp)
        fit_surv2 = coxph(Surv(days,cens)~arms, data = data_tmp)
        pval_Y1_prob = as.numeric(na.omit(stats::anova(fit_surv2, fit_surv3)[[4]]))

        if(pval_M_prob < pval_flag_prob){
          pval_flag_prob = pval_M_prob
          result[[1]] = profiles
          result[[2]] = pval_flag_prob
          result[[3]] = pval_Y_prob
          result[[4]] = pval_Y1_prob
          result[[5]] = l
          result[[6]] = r
        }
      }

    }

  }

  if (result[[2]]<=p_min_r){
    p_min_r = result[[2]]
    result_flag = result
    data_flag = data
  }

}

plot_profile(result_flag[[1]]$trees[[result_flag[[5]]]])
data = data_flag
summary(data$indir_arms_effect)
quantile(data$indir_arms_effect,0.9)
quantile(data$indir_arms_effect,0.1)

data$indir_arms_effect = pmin(data$indir_arms_effect ,0.008)
data$indir_arms_effect = pmax(data$indir_arms_effect ,-1.13)



data[,"subgroup"] = rep(NA, 1054)
for (i in 1:1054){
  if(data[i,"cd40"]<298){
    data[i,"subgroup"] = "subgroup1"
  }else if(data[i,"preanti"]<8&data[i,"cd40"]>=298){
    data[i,"subgroup"] = "subgroup2"
  }else if(data[i,"preanti"]>=8&data[i,"cd40"]>=298){
    #print(i)
    data[i,"subgroup"] = "subgroup3"
  }
}


data %>%
  group_by(subgroup) %>%
  summarise(mean_te = mean(indir_arms_effect, na.rm = TRUE))


data$subtype = as.factor(data$subtype)
data$subgroup = as.factor(data$subgroup)


cols <- c(
  "subgroup1" = "#1b9e77",
  "subgroup2" = "#d95f02",
  "subgroup3" = "#7570b3"
)

summary(data$subtype)
boxplot = ggplot(data, aes(x = subgroup, y = indir_arms_effect)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = subgroup),width = 0.15, alpha = 0.4, size = 1) +
  theme_bw() +
  labs(
    x = "",
    y = "Indirect Treatment Effect Contional on Covariates"
  )+ stat_summary(fun = mean, geom = "point",
                  shape = 18, size = 3, color = "red")+
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),      #
    axis.title = element_text(size = 16),     #
    text = element_text(size = 16)            #
  )


