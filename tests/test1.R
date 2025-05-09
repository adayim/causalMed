#数据准备
devtools::install_github("CausalInference/gfoRmula")
library(gfoRmula)

data(binary_eofdata, package = "gfoRmula")

binary_eofdata$time_f <- ifelse(binary_eofdata$time <= 1,
                                0,
                                ifelse(
                                  binary_eofdata$time <= 3,
                                  1,
                                  ifelse(binary_eofdata$time <= 5, 2, 3)
                                ))
data2 <- binary_eofdata
data2$treat <- ifelse(is.na(data2$treat), NA, ifelse(data2$treat > 0, 1, 0))
data2$time_f <- as.factor(data2$time_f)

# remotes::install_github("adayim/causalMed", ref = "dev")
library(causalMed)
library(dplyr)
library(data.table)
library(progressr)
library(progress)
library(data.table)
data1 <- binary_eofdata
data1 <- data1 %>%
  group_by(id_num) %>%
  mutate(
    lag1_treat = lag(treat, 1),
    # treat 的滞后 1 期
    lag1_cov1 = lag(cov1, 1),
    # cov1 的滞后 1 期
    cumavg_treat = cummean(ifelse(is.na(treat), 0, treat)),
    # treat 的累积平均值
    cumavg_cov1 = cummean(ifelse(is.na(cov1), 0, cov1))     # cov1 的累积平均值
  ) %>%
  ungroup()
data1$treat <- ifelse(is.na(data1$treat), NA, ifelse(data1$treat > 0, 1, 0))
data1$time <- as.numeric(binary_eofdata$time)
data1$treat2 <- as.numeric(binary_eofdata$treat)
data1 <- data1 %>%
  mutate(
    lag1_treat = ifelse(time == 1, NA, lag1_treat),
    lag1_cov1 = ifelse(time == 1, NA, lag1_cov1)
  )

# Parallel back end for windows
plan(multisession, workers = availableCores() - 2)

# Parallel back end for non-windows
# plan(multicore, workers = availableCores() - 2)

#测试1：初始测试
results <- data.table(
  seed = integer(),
  bz_natural = numeric(),
  bz_intervention = numeric(),
  cs_natural = numeric(),
  cs_intervention = numeric(),
  q_n = numeric(),
  q_i = numeric()
)
for (seed_num in 1:100) {
  mod1 <- spec_model(
    cov1 ~ lag1_treat + lag1_cov1  + cov3 + time_f,
    var_type  = "binomial",
    mod_type = "covariate"
  )
  mod2 <- spec_model(
    treat ~ lag1_treat + cumavg_cov1 + cov3 + time_f,
    var_type = "binomial",
    mod_type = "exposure"
  )
  #             spec_model(,family = "binomial",order = ,type = ""),
  #             spec_model(,family = "binomial",order = ,type = ""),
  mod3 <- spec_model(outcome ~  treat + cov1 + lag1_cov1 + cov3,
                     var_type = "binomial",
                     mod_type = "outcome")
  models1 <- list(mod1, mod2, mod3)
  res2 <- gformula(
    data = data1,
    id_var = 'id_num',
    base_vars = "cov3",
    exposure = "treat",
    time_var = "time",
    models = models1,
    intervention = list(natural = NULL, # 自然暴露（不干预）
                        always = 1,    # 始终干预（A=1）
                        #never = 0     # 从不干预（A=0）
                        ref_int = 1),
                        init_recode = c(
                          "lag1_treat = 0",
                          "lag1_cov1=0",
                          "time_f=0",
                          "cumavg_cov1=0",
                          "cumavg_treat=0"
                        ),
                        in_recode = NULL,
                        out_recode = c(
                          "cumavg_treat=cumavg_treat+treat",
                          "cumavg_cov1=cumavg_cov1+cov1",
                          "cov1=cov1",
                          "time_f=time_f"
                        ),
                        return_fitted = T,
                        mc_sample = 250000,
                        return_data = T,
                        R = 50,
                        quiet = T,
                        seed = seed_num
    )
    ymodel <- outcome ~  treat + cov1 + lag1_cov1 + cov3
    res <- gfoRmula::gformula(
      obs_data = data2,
      outcome_type = 'binary_eof',
      id = 'id_num',
      time_name = 'time',
      covnames = c('cov1', 'treat', 'time_f'),
      outcome_name = 'outcome',
      covtypes = c('binary', 'normal', 'categorical time'),
      covparams = list(
        covmodels = c(
          cov1 ~ lag1_treat + lag1_cov1 + cov3 + time_f,
          treat ~ lag1_treat + cumavg_cov1 +
            cov3 + time_f,
          NA
        )
      ),
      ymodel = ymodel,
      intervention1.treat = list(threshold, 1, Inf),
      int_descript = 'Threshold - lower bound 1',
      histories = c(lagged, cumavg),
      histvars = list(c('treat', 'cov1'), c('cov1')),
      basecovs = c("cov3"),
      seed = seed_num,
      parallel = TRUE,
      nsimul = 250000,
      ncores = 2
    )
    result_dt1 <- res$result[["g-form mean"]]
    result_dt2 <- res2$effect_size[["Est"]]
    bz_natural <- result_dt1[[1]]
    bz_intervention <- result_dt1[[2]]
    cs_natural <- result_dt2[[2]]
    cs_intervention <- result_dt2[[1]]

    q_n <- (cs_natural - bz_natural) / bz_natural * 100
    q_i <- (cs_intervention - bz_intervention) / bz_intervention * 100
    results <- rbind(
      results,
      data.table(
        seed = seed_num,
        bz_natural = bz_natural,
        bz_intervention = bz_intervention,
        cs_natural = cs_natural,
        cs_intervention = cs_intervention,
        q_n = q_n,
        q_i = q_i
      )
    )
    }
print(results)


#测试2：删除lag和cum
results <- data.table(
  seed = integer(),
  bz_natural = numeric(),
  bz_intervention = numeric(),
  cs_natural = numeric(),
  cs_intervention = numeric(),
  q_n = numeric(),
  q_i = numeric()
)
for (seed_num in 1:100) {
  mod1 <- spec_model(cov1 ~ cov3 + time_f, var_type  = "binomial", mod_type = "covariate")
  mod2 <- spec_model(treat ~ cov3 + time_f,
                     var_type = "binomial",
                     mod_type = "exposure")
  mod3 <- spec_model(outcome ~  treat + cov1 +  cov3,
                     var_type = "binomial",
                     mod_type = "outcome")
  models1 <- list(mod1, mod2, mod3)
  res2 <- gformula(
    data = data1,
    id_var = 'id_num',
    base_vars = "cov3",
    exposure = "treat",
    time_var = "time",
    models = models1,
    intervention = list(natural = NULL, # 自然暴露（不干预）
                        always = 1    # 始终干预（A=1）
                        #never = 0     # 从不干预（A=0）
                        ),
                        ref_int = 1,
                        init_recode = c("time_f=0"),
                        in_recode = NULL,
                        out_recode = NULL,
                        return_fitted = T,
                        mc_sample = 2501,
                        return_data = T,
                        R = 500,
                        quiet = T,
                        seed = seed_num
    )
    ymodel <- outcome ~  treat + cov1 +  cov3
    res <- gfoRmula::gformula(
      obs_data = data2,
      outcome_type = 'binary_eof',
      id = 'id_num',
      time_name = 'time',
      covnames = c('cov1', 'treat', 'time_f'),
      outcome_name = 'outcome',
      covtypes = c('binary', 'normal', 'categorical time'),
      covparams = list(covmodels = c(
        cov1 ~  cov3 + time_f, treat ~ cov3 + time_f, NA
      )),
      ymodel = ymodel,
      intervention1.treat = list(gfoRmula::static, c(1, 1, 1, 1, 1, 1, 1)),
      #int_descript = c('Always treat (A=1)', 'Never treat (A=0)'),#描述干预，不影响模型的计算或结果
      basecovs = c("cov3"),
      seed = seed_num,
      parallel = TRUE,
      nsimul = 2501,
      ncores = 2
    )
    result_dt1 <- res$result[["g-form mean"]]
    result_dt2 <- res2$effect_size[["Est"]]
    bz_natural <- result_dt1[[1]]
    bz_intervention <- result_dt1[[2]]
    cs_natural <- result_dt2[[2]]
    cs_intervention <- result_dt2[[1]]

    q_n <- (cs_natural - bz_natural) / bz_natural * 100
    q_i <- (cs_intervention - bz_intervention) / bz_intervention * 100
    results <- rbind(
      results,
      data.table(
        seed = seed_num,
        bz_natural = bz_natural,
        bz_intervention = bz_intervention,
        cs_natural = cs_natural,
        cs_intervention = cs_intervention,
        q_n = q_n,
        q_i = q_i
      )
    )
    }
print(results)




#测试3：删除cum
results <- data.table(
  seed = integer(),
  bz_natural = numeric(),
  bz_intervention = numeric(),
  cs_natural = numeric(),
  cs_intervention = numeric(),
  q_n = numeric(),
  q_i = numeric()
)
for (seed_num in 1:100) {
  mod1 <- spec_model(
    cov1 ~ lag1_treat + lag1_cov1  + cov3 + time_f,
    var_type  = "binomial",
    mod_type = "covariate"
  )
  mod2 <- spec_model(treat ~ lag1_treat + cov3 + time_f,
                     var_type = "binomial",
                     mod_type = "exposure")
  mod3 <- spec_model(outcome ~  treat + cov1 + lag1_cov1 + cov3,
                     var_type = "binomial",
                     mod_type = "outcome")
  models1 <- list(mod1, mod2, mod3)
  res2 <- gformula(
    data = data1,
    id_var = 'id_num',
    base_vars = "cov3",
    exposure = "treat",
    time_var = "time",
    models = models1,
    intervention = list(natural = NULL, # 自然暴露（不干预）
                        always = 1    # 始终干预（A=1）
                        #never = 0     # 从不干预（A=0）
                        ),
                        ref_int = 1,
                        init_recode = c("time_f=0", "lag1_treat=0", "lag1_cov1=0"),
                        in_recode = NULL,
                        out_recode = c(
                          "treat=treat",
                          "cov1=cov1",
                          "lag1_treat=treat",
                          "lag1_cov1=cov1"
                        ),
                        return_fitted = T,
                        mc_sample = 2501,
                        return_data = T,
                        R = 50,
                        quiet = T,
                        seed = seed_num
    )
    ymodel <- outcome ~  treat + cov1 + lag1_cov1 + cov3
    res <- gfoRmula::gformula(
      obs_data = data2,
      outcome_type = 'binary_eof',
      id = 'id_num',
      time_name = 'time',
      covnames = c('cov1', 'treat', 'time_f'),
      outcome_name = 'outcome',
      covtypes = c('binary', 'normal', 'categorical time'),
      covparams = list(
        covmodels = c(
          cov1 ~ lag1_treat + lag1_cov1 + cov3 + time_f,
          treat ~ lag1_treat + cov3 + time_f,
          NA
        )
      ),
      ymodel = ymodel,
      intervention1.treat = list(gfoRmula::static, c(1, 1, 1, 1, 1, 1, 1)),
      #int_descript = c('Always treat (A=1)', 'Never treat (A=0)'),#描述干预，不影响模型的计算或结果
      histories = c(gfoRmula::lagged),
      histvars = list(c('treat', 'cov1')),
      basecovs = c("cov3"),
      seed = seed_num,
      parallel = TRUE,
      nsimul = 2501,
      ncores = 2
    )

    result_dt1 <- res$result[["g-form mean"]]
    result_dt2 <- res2$effect_size[["Est"]]
    bz_natural <- result_dt1[[1]]
    bz_intervention <- result_dt1[[2]]
    cs_natural <- result_dt2[[2]]
    cs_intervention <- result_dt2[[1]]

    q_n <- (cs_natural - bz_natural) / bz_natural * 100
    q_i <- (cs_intervention - bz_intervention) / bz_intervention * 100
    results <- rbind(
      results,
      data.table(
        seed = seed_num,
        bz_natural = bz_natural,
        bz_intervention = bz_intervention,
        cs_natural = cs_natural,
        cs_intervention = cs_intervention,
        q_n = q_n,
        q_i = q_i
      )
    )
    }
print(results)



#测试4：删除lag，cum，time_f
results <- data.table(
  seed = integer(),
  bz_natural = numeric(),
  bz_intervention = numeric(),
  cs_natural = numeric(),
  cs_intervention = numeric(),
  q_n = numeric(),
  q_i = numeric()
)
for (seed_num in 1:100) {
  mod1 <- spec_model(cov1 ~ cov3, var_type  = "binomial", mod_type = "covariate")
  mod2 <- spec_model(treat ~ cov3, var_type = "binomial", mod_type = "exposure")
  mod3 <- spec_model(outcome ~  treat + cov1 +  cov3,
                     var_type = "binomial",
                     mod_type = "outcome")
  models1 <- list(mod1, mod2, mod3)
  res2 <- gformula(
    data = data1,
    id_var = 'id_num',
    base_vars = "cov3",
    exposure = "treat",
    time_var = "time",
    models = models1,
    intervention = list(natural = NULL, # 自然暴露（不干预）
                        always = 1    # 始终干预（A=1）
                        #never = 0     # 从不干预（A=0）
                        ),
                        ref_int = 1,
                        init_recode = NULL,
                        in_recode = NULL,
                        out_recode = NULL,
                        return_fitted = T,
                        mc_sample = 2501,
                        return_data = T,
                        R = 100,
                        quiet = T,
                        seed = seed_num
    )
    ymodel <- outcome ~  treat + cov1 +  cov3
    res <- gfoRmula::gformula(
      obs_data = data2,
      outcome_type = 'binary_eof',
      id = 'id_num',
      time_name = 'time',
      covnames = c('cov1', 'treat'),
      outcome_name = 'outcome',
      covtypes = c('binary', 'normal'),
      covparams = list(covmodels = c(cov1 ~  cov3, treat ~ cov3)),
      ymodel = ymodel,
      intervention1.treat = list(gfoRmula::static, c(1, 1, 1, 1, 1, 1, 1)),
      #int_descript = c('Always treat (A=1)', 'Never treat (A=0)'),#描述干预，不影响模型的计算或结果
      basecovs = c("cov3"),
      seed = seed_num,
      parallel = TRUE,
      nsimul = 2501,
      ncores = 2
    )
    result_dt1 <- res$result[["g-form mean"]]
    result_dt2 <- res2$effect_size[["Est"]]
    bz_natural <- result_dt1[[1]]
    bz_intervention <- result_dt1[[2]]
    cs_natural <- result_dt2[[2]]
    cs_intervention <- result_dt2[[1]]

    q_n <- (cs_natural - bz_natural) / bz_natural * 100
    q_i <- (cs_intervention - bz_intervention) / bz_intervention * 100
    results <- rbind(
      results,
      data.table(
        seed = seed_num,
        bz_natural = bz_natural,
        bz_intervention = bz_intervention,
        cs_natural = cs_natural,
        cs_intervention = cs_intervention,
        q_n = q_n,
        q_i = q_i
      )
    )
    }
print(results)
