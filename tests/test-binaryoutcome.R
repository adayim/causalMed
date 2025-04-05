#data preparation
devtools::install_github("CausalInference/gfoRmula")
library(gfoRmula)
load("C:/Users/16074/AppData/Local/Temp/MicrosoftEdgeDownloads/1b793faa-5323-43a2-9da1-1eb8ced81ace/binary_eofdata.rda")
binary_eofdata$time_f <- ifelse(binary_eofdata$time <= 1, 0,
                                ifelse(binary_eofdata$time <= 3, 1,
                                       ifelse(binary_eofdata$time <= 5, 2, 3)))
data2<-binary_eofdata
data2$treat <- ifelse(is.na(data2$treat), NA, ifelse(data2$treat > 0, 1, 0))
data2$time_f<-as.factor(data2$time_f)

remotes::install_github("adayim/causalMed", ref = "dev")
library(causalMed)
library(dplyr)
library(data.table)
library(progressr)
library(progress)
library(data.table)
data1<-binary_eofdata
data1$treat <- ifelse(is.na(data1$treat), NA, ifelse(data1$treat > 0, 1, 0))
data1 <- data1 %>%
  group_by(id_num) %>%  
  mutate(
    lag1_treat = lag(treat, 1),  
    lag1_cov1 = lag(cov1, 1),    
    cumavg_treat = cummean(ifelse(is.na(treat), 0, treat)),  
    cumavg_cov1 = cummean(ifelse(is.na(cov1), 0, cov1))     
  ) %>%
  ungroup()
data1$time<-as.numeric(binary_eofdata$time)
data1$treat2<-as.numeric(binary_eofdata$treat)
data1$time_f<-as.factor(data1$time_f)
data1 <- data1 %>%
  mutate(
    lag1_treat = ifelse(time == 1, NA, lag1_treat),
    lag1_cov1 = ifelse(time == 1, NA, lag1_cov1)
  )

#test for binary outcome
results <- data.table(
  seed = integer(),
  bz_natural = numeric(),
  bz_intervention = numeric(),
  cs_natural = numeric(),
  cs_intervention = numeric(),
  q_n=numeric(),
  q_i=numeric()
)
for (seed_num in 500:650) {
  mod1<-spec_model(cov1 ~ lag1_treat + lag1_cov1  +cov3 + time_f,var_type  = "binomial",mod_type = "covariate")
  mod2<-spec_model(treat ~ lag1_treat + cumavg_cov1 +cov3 + time_f,var_type= "binomial",mod_type = "exposure")
  mod3<-spec_model(outcome ~  treat + cov1 + lag1_cov1 + cov3,var_type = "binomial",mod_type = "outcome")
  models1<-list(mod1,mod2,mod3)
  res2<- gformula(data=data1,
                  id_var='id_num',
                  base_vars="cov3",
                  exposure="treat",
                  time_var="time",
                  models=models1,
                  intervention = list(
                    natural = NULL,
                    always = 1
                  ),
                  ref_int = 1,
                  init_recode = c("lag1_treat = 0","lag1_cov1=0","time_f=factor(0,levels=levels(data1$time_f))","cumavg_cov1=0","cumavg_treat=0"),
                  in_recode = c("lag1_treat = treat","lag1_cov1=cov1","time_f=as.factor(ifelse(time <= 1, 0,
                                ifelse(time <= 3, 1,
                                       ifelse(time <= 5, 2, 3))))","cumavg_cov1=cummean(ifelse(is.na(cov1), 0, cov1))","cumavg_treat=cummean(ifelse(is.na(treat), 0, treat))"),
                  out_recode = NULL,
                  return_fitted = T,
                  mc_sample = 100000,
                  return_data = T,
                  R = 5,
                  quiet = T,
                  seed=seed_num)  
  ymodel <- outcome ~  treat + cov1 + lag1_cov1 + cov3
  res<- gfoRmula::gformula(obs_data = data2,
                           outcome_type = 'binary_eof', id = 'id_num',
                           time_name = 'time', covnames = c('cov1', 'treat', 'time_f'),
                           outcome_name = 'outcome', covtypes = c('binary', 'binary', 'categorical time'),
                           covparams = list(covmodels = c(cov1 ~ lag1_treat + lag1_cov1 + cov3 + time_f,
                                                          treat ~ lag1_treat + cumavg_cov1 +cov3 + time_f,
                                                          NA)), 
                           ymodel = ymodel,
                           intervention1.treat = list(threshold, 1, Inf),
                           int_descript = 'Threshold - lower bound 1',
                           histories = c(lagged, cumavg),
                           histvars = list(c('treat', 'cov1'), c('cov1')), basecovs = c("cov3"),
                           seed = seed_num, parallel = TRUE,
                           nsimul = 100000, ncores = 2,sim_data_b = T)
  result_dt1 <- res$result[["g-form mean"]]
  result_dt2<-res2$effect_size[["Est"]]
  bz_natural <- result_dt1[[1]]
  bz_intervention <- result_dt1[[2]]
  cs_natural<-result_dt2[[2]]
  cs_intervention <- result_dt2[[1]]
  
  q_n <- (cs_natural - bz_natural) / bz_natural * 100
  q_i <- (cs_intervention - bz_intervention) / bz_intervention * 100
  results <- rbind(results, data.table(
    seed = seed_num,
    bz_natural = bz_natural,
    bz_intervention = bz_intervention,
    cs_natural = cs_natural,
    cs_intervention = cs_intervention,
    q_n=q_n,
    q_i=q_i
  ))
}
print(results)
mean(results$bz_natural)
mean(results$bz_intervention)
mean(results$cs_natural)
mean(results$cs_intervention)
mean(results$q_n)
max(results$q_n)
min(results$q_n)
mean(results$q_i)
max(results$q_i)
min(results$q_i)
