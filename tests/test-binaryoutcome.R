#data preparation
devtools::install_github("CausalInference/gfoRmula")
library(gfoRmula)
data(binary_eofdata)
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
data1$time_f<-as.numeric(data1$time_f)
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
for (seed_num in 600:650) {
  mod1<-spec_model(cov1 ~ lag1_treat + lag1_cov1  +cov3 + time_f,var_type  = "binomial",mod_type = "covariate")
  mod2<-spec_model(treat ~ lag1_treat + cumavg_cov1 +cov3 + time_f,var_type= "binomial",mod_type = "exposure")
  mod3<-spec_model(outcome ~  lag1_treat+treat + cov1 + lag1_cov1 + cov3,var_type = "binomial",mod_type = "outcome")
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
                  init_recode = c("lag1_treat = 0","lag1_cov1=0","time_f=0","cumavg_cov1=0","cumavg_treat=0"),
                  in_recode = c("lag1_treat = treat","lag1_cov1=cov1","time_f=ifelse(time <= 1, 0, ifelse(time <= 3, 1, ifelse(time <= 5, 2,3)))",
                                "cumavg_cov1=cummean(ifelse(is.na(cov1), 0, cov1))","cumavg_treat=cummean(ifelse(is.na(treat), 0, treat))"),
                  out_recode = NULL,
                  return_fitted = F,
                  mc_sample = 50000,
                  return_data = F,
                  R = 5,
                  quiet = T,
                  seed=seed_num)  
  ymodel <- outcome ~  lag1_treat+treat + cov1 + lag1_cov1 + cov3
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
  cs_natural<-result_dt2[[1]]
  cs_intervention <- result_dt2[[2]]
  
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
sum <- data.frame(
  mean = sapply(results, mean, na.rm = TRUE),
  sd   = sapply(results, sd,   na.rm = TRUE),
  min  = sapply(results, min,  na.rm = TRUE),
  max  = sapply(results, max,  na.rm = TRUE)
)  
sum
##############################################survival
data(basicdata_nocomp)

datanoco<-subset(basicdata_nocomp,t0<5)
 

  datanocomp<-basicdata_nocomp %>%
    arrange(id, t0) %>%  
    mutate(
      lag1_A = lag(A, 1),  
      lag1_L1 = lag(L1, 1),
      lag1_L2 = lag(L2, 1)
    ) %>%
    ungroup()
  datanocom<-subset(datanocomp,t0<5)
    
   mod1<-spec_model(L1 ~  lag1_L2 +L3 + t0,var_type  = "binomial",mod_type = "covariate")
   mod2<-spec_model(L2 ~ lag1_L1 +lag1_L2 + L3 + t0,var_type  = "normal",mod_type = "covariate")
   
    mod3<-spec_model(A ~  L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0,var_type= "binomial",mod_type = "exposure")
    mod4<-spec_model(Y ~ A + L1 + L2 + L3 +lag1_L1 + lag1_L2 + t0,var_type = "binomial",mod_type = "survival")
    models1<-list(mod1,mod2,mod3,mod4)
    
      ymodel <- Y ~ A + L1 + L2 + L3 + lag1_L1 + lag1_L2 + t0

    
    results <- data.table(
      seed = integer(),
      bz_natural = numeric(),
      bz_intervention1 = numeric(),
      bz_intervention2 = numeric(),
      cs_natural = numeric(),
      cs_intervention1 = numeric(),
      cs_intervention2 = numeric(),
      q_n=numeric(),
      q_i1=numeric(),
      q_i2=numeric()
    )  
  for (seed_num in 500:510) {
    res2<- gformula(data=datanocom,
                    id_var='id',
                    base_vars="L3",
                    exposure="A",
                    time_var="t0",
                    models=models1,
                    intervention = list(
                      natural = NULL,
                      never = 0,
                      always = 1
                    ),
                    #ref_int = 1,
                    init_recode = c("lag1_A = 0","lag1_L1=0","lag1_L2=0"),
                    in_recode =  c("lag1_A = A","lag1_L1=L1","lag1_L2=L2"),
                    out_recode =NULL,
                    return_fitted = T,
                    mc_sample = 10305,
                    return_data = T,
                    R = 5,
                    quiet = T,
                    seed=seed_num) 

    
    res <- gfoRmula::gformula(obs_data = datanoco, id = 'id',
                            time_points = 5,
                            time_name = 't0', covnames = c('L1', 'L2', 'A'),
                            outcome_name = 'Y',
                            outcome_type = 'survival', covtypes = c('binary', 'normal', 'binary'),
                            covparams = list(covmodels = c(L1 ~   lag1_L2 +L3 + t0,
                                                           L2 ~  lag1_L1 +lag1_L2 + L3 + t0,
                                                           A ~  L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0)), 
                            ymodel = ymodel,
                            intervention1.A = list(static, rep(0, 5)),
                            intervention2.A = list(static, rep(1, 5)),
                            int_descript = c('Never treat', 'Always treat'),
                            histories =  c(lagged), histvars = list(c('A', 'L1', 'L2')),
                            basecovs = c('L3'), intcomp = c(1,2),
                            nsimul = 10305,sim_data_b = T,
                            seed = seed_num, nsamples = 0)
                            
    result_dt1 <- res$result[k==max(k),][["g-form risk"]]
    result_dt2<-res2$effect_size[["Est"]]
    bz_natural <- result_dt1[[1]]
    bz_intervention1 <- result_dt1[[2]]
    bz_intervention2 <- result_dt1[[3]]
    cs_natural<-result_dt2[[1]]
    cs_intervention1 <- result_dt2[[2]]
    cs_intervention2 <- result_dt2[[3]]
    q_n <- (cs_natural - bz_natural) / bz_natural * 100
    q_i1 <- (cs_intervention1 - bz_intervention1) / bz_intervention1 * 100
    q_i2 <- (cs_intervention2 - bz_intervention2) / bz_intervention2 * 100
    results <- rbind(results, data.table(
      seed = seed_num,
      bz_natural = bz_natural,
      bz_intervention1 = bz_intervention1,
      bz_intervention2 = bz_intervention2,
      cs_natural = cs_natural,
      cs_intervention1 = cs_intervention1,
      cs_intervention2 = cs_intervention2,
      q_n=q_n,
      q_i1=q_i1,
      q_i2=q_i2
    ))
  }
    #results
    summ <- data.frame(
      mean = sapply(results, mean, na.rm = TRUE),
      sd   = sapply(results, sd,   na.rm = TRUE),
      min  = sapply(results, min,  na.rm = TRUE),
      max  = sapply(results, max,  na.rm = TRUE)
    )  
    summ
    
    table(res$sim_data$`Natural course`[t0==4,A])
    table(res2$sim_data[Intervention=="natural",A])
    
    table(res$sim_data$`Natural course`[t0==4,lag1_A])
    table(res2$sim_data[Intervention=="natural",lag1_A])
    
    res$sim_data$`Natural course`[t0==4,]
    res2$sim_data[Intervention=="natural",]
    
    res$sim_data$`Never treat`[t0==4,]
    res2$sim_data[Intervention=="never",]
    
    
    ########################################################3
    
    
    datanoco<-subset(basicdata_nocomp,t0<3)
        datanocom<-subset(datanocomp,t0<3)

    
    ymodel <- Y ~ A + L1 + L2 + L3 +  lag1_L1 + lag1_L2 + t0
    
    mod1<-spec_model(L1 ~   lag1_L2 +L3 + t0,var_type  = "binomial",mod_type = "covariate")
    mod2<-spec_model(L2 ~  lag1_L1 +lag1_L2 + L3 + t0,var_type  = "normal",mod_type = "covariate")
    
    mod3<-spec_model(A ~ L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0,var_type= "binomial",mod_type = "exposure")
    mod4<-spec_model(Y ~ A + L1 + L2 + L3 +  lag1_L1 + lag1_L2 + t0,var_type = "binomial",mod_type = "survival")
    models1<-list(mod1,mod2,mod3,mod4)
    
    results <- data.table(
      seed = integer(),
      bz_natural = numeric(),
      bz_intervention1 = numeric(),
      bz_intervention2 = numeric(),
      cs_natural = numeric(),
      cs_intervention1 = numeric(),
      cs_intervention2 = numeric(),
      q_n=numeric(),
      q_i1=numeric(),
      q_i2=numeric()
    )  
    for (seed_num in 500:520) {
      res2<- gformula(data=datanocom,
                      id_var='id',
                      base_vars="L3",
                      exposure="A",
                      time_var="t0",
                      models=models1,
                      intervention = list(
                        natural = NULL,
                        never = 0,
                        always = 1
                      ),
                      #ref_int = 1,
                      init_recode = c("lag1_L1=0","lag1_L2=0"),
                      in_recode = NULL,
                      out_recode = c("lag1_L1=L1","lag1_L2=L2"),
                      return_fitted = T,
                      mc_sample = 20000,
                      return_data = T,
                      R = 5,
                      quiet = T,
                      seed=seed_num) 
      
      
      res <- gfoRmula::gformula(obs_data = datanoco, id = 'id',
                                time_points = 5,
                                time_name = 't0', covnames = c('L1', 'L2', 'A'),
                                outcome_name = 'Y',
                                outcome_type = 'survival', covtypes = c('binary', 'normal', 'binary'),
                                covparams = list(covmodels = c(L1 ~   lag1_L2 +L3 + t0,
                                                               L2 ~ lag1_L1 +lag1_L2 + L3 + t0,
                                                               A ~  L1 + L2 + lag1_L1 + lag1_L2 + L3 + t0)), 
                                ymodel = ymodel,
                                intervention1.A = list(static, rep(0, 5)),
                                intervention2.A = list(static, rep(1, 5)),
                                int_descript = c('Never treat', 'Always treat'),
                                histories =  c(lagged), histvars = list(c( 'L1', 'L2')),
                                basecovs = c('L3'), intcomp = c(1,2),
                                nsimul = 10000,
                                seed = seed_num, nsamples = 2)
      
     result_dt1 <- res$result[k==max(k),][["g-form risk"]]
    result_dt2<-res2$effect_size[["Est"]]
    bz_natural <- result_dt1[[1]]
    bz_intervention1 <- result_dt1[[2]]
    bz_intervention2 <- result_dt1[[3]]
    cs_natural<-result_dt2[[1]]
    cs_intervention1 <- result_dt2[[2]]
    cs_intervention2 <- result_dt2[[3]]
    q_n <- (cs_natural - bz_natural) / bz_natural * 100
    q_i1 <- (cs_intervention1 - bz_intervention1) / bz_intervention1 * 100
    q_i2 <- (cs_intervention2 - bz_intervention2) / bz_intervention2 * 100
    results <- rbind(results, data.table(
      seed = seed_num,
      bz_natural = bz_natural,
      bz_intervention1 = bz_intervention1,
      bz_intervention2 = bz_intervention2,
      cs_natural = cs_natural,
      cs_intervention1 = cs_intervention1,
      cs_intervention2 = cs_intervention2,
      q_n=q_n,
      q_i1=q_i1,
      q_i2=q_i2
    ))
    }
    summ <- data.frame(
      mean = sapply(results, mean, na.rm = TRUE),
      sd   = sapply(results, sd,   na.rm = TRUE),
      min  = sapply(results, min,  na.rm = TRUE),
      max  = sapply(results, max,  na.rm = TRUE)
    )  
    summ 
    
    ###############去掉L2还是不行
    datanoco<-subset(basicdata_nocomp,t0<5)
        datanocom<-subset(datanocomp,t0<5)

    
    ymodel <- Y ~ A + L1  + L3 +  lag1_L1  + t0
    
    mod1<-spec_model(L1 ~   lag1_L1  +L3 + t0,var_type  = "binomial",mod_type = "covariate")
    #mod2<-spec_model(L2 ~  lag1_L1  + L3 + t0,var_type  = "normal",mod_type = "covariate")
    
    mod3<-spec_model(A ~ L1  + lag1_L1  + L3 + t0,var_type= "binomial",mod_type = "exposure")
    mod4<-spec_model(Y ~ A + L1  + L3 +  lag1_L1  + t0,var_type = "binomial",mod_type = "survival")
    models1<-list(mod1,mod3,mod4)#mod2,
    
    results <- data.table(
      seed = integer(),
      bz_natural = numeric(),
      bz_intervention1 = numeric(),
      bz_intervention2 = numeric(),
      cs_natural = numeric(),
      cs_intervention1 = numeric(),
      cs_intervention2 = numeric(),
      q_n=numeric(),
      q_i1=numeric(),
      q_i2=numeric()
    )  
    for (seed_num in 500:510) {
      res2<- gformula(data=datanocom,
                      id_var='id',
                      base_vars="L3",
                      exposure="A",
                      time_var="t0",
                      models=models1,
                      intervention = list(
                        natural = NULL,
                        never = 0,
                        always = 1
                      ),
                      #ref_int = 1,
                      init_recode = c("lag1_A = 0","lag1_L1=0"),
                      in_recode = c("lag1_A = A","lag1_L1=L1"),
                      out_recode = NULL,
                      return_fitted = T,
                      mc_sample = 200000,
                      return_data = T,
                      R = 5,
                      quiet = T,
                      seed=seed_num) 
      
      
      res <- gfoRmula::gformula(obs_data = datanoco, id = 'id',
                                time_points = 5,
                                time_name = 't0', covnames = c('L1',  'A'),
                                outcome_name = 'Y',
                                outcome_type = 'survival', covtypes = c('binary',  'binary'),
                                covparams = list(covmodels = c(L1 ~  lag1_L1    +L3 + t0,
                                                               A ~  L1 + lag1_L1  + L3 + t0)), 
                                ymodel = ymodel,
                                intervention1.A = list(static, rep(0, 5)),
                                intervention2.A = list(static, rep(1, 5)),
                                int_descript = c('Never treat', 'Always treat'),
                                histories =  c(lagged), histvars = list(c( 'L1')),
                                basecovs = c('L3'), #intcomp = c(1,2),
                                nsimul = 200000,sim_data_b = T,
                                seed = seed_num, nsamples = 0)
      
     result_dt1 <- res$result[k==max(k),][["g-form risk"]]
    result_dt2<-res2$effect_size[["Est"]]
    bz_natural <- result_dt1[[1]]
    bz_intervention1 <- result_dt1[[2]]
    bz_intervention2 <- result_dt1[[3]]
    cs_natural<-result_dt2[[1]]
    cs_intervention1 <- result_dt2[[2]]
    cs_intervention2 <- result_dt2[[3]]
    q_n <- (cs_natural - bz_natural) / bz_natural * 100
    q_i1 <- (cs_intervention1 - bz_intervention1) / bz_intervention1 * 100
    q_i2 <- (cs_intervention2 - bz_intervention2) / bz_intervention2 * 100
    results <- rbind(results, data.table(
      seed = seed_num,
      bz_natural = bz_natural,
      bz_intervention1 = bz_intervention1,
      bz_intervention2 = bz_intervention2,
      cs_natural = cs_natural,
      cs_intervention1 = cs_intervention1,
      cs_intervention2 = cs_intervention2,
      q_n=q_n,
      q_i1=q_i1,
      q_i2=q_i2
    ))
    }
    summ <- data.frame(
      mean = sapply(results, mean, na.rm = TRUE),
      sd   = sapply(results, sd,   na.rm = TRUE),
      min  = sapply(results, min,  na.rm = TRUE),
      max  = sapply(results, max,  na.rm = TRUE)
    )  
    summ 
    