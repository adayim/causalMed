

################################
# calc_iptw()
################################
#' Calculate total effect, natural direct and indirect effect.
#'
#' @param data Data to be used to ceate matrix.
#' @param id.vars ID variable name.
#' @param time.var Time indexing varaible name.
#' @param baseline.vars A vector of baseline variables and time-fixed variable names.
#' @param exposure.model Exposure model, only binomial outcome supported.
#' @param exposure.vars A vector of exposure variables, including lagged exposure variables.
#' @param mediator.model Mediator model.
#' @param mediator.family Model link function, \code{gussian},\code{binomial} and \code{multinomial}.
#' For multinomial model, \code{\link[nnet]{multinom}} will be sued for model fitting. And
#' \code{\link[stats]{glm}} for  \code{gussian} and \code{binomial} link funciton model.
#' @param outcome Outcome variable name.
#' @param is.expL Is time-time varying exposure, defualt is TRUE.
#' @param censor.model Censor model. Censoring indicator should be coded as remaing uncensored (as 1).
#' @param bound.value Bound weights, between 0 and 0.5.
#'
#'
#' @import data.table
#'
calc_iptw <- function(data,
                      id.vars,
                      time.var,
                      baseline.vars,
                      exposure.model,
                      exposure.vars,
                      mediator.model,
                      mediator.family,
                      outcome,
                      is.expL = TRUE,
                      censor.model = NULL,
                      bound.value = 0){

  tpcall <- match.call()

  if(bound.value <0 | bound.value > 0.5)
    stop("bound.value should be between 0 and 0.5!")

  id.vars <- substitute(id.vars)
  time.var <- substitute(time.var)
  outcome <- substitute(outcome)

  # Select all variables
  all_vars <- unique(c(id.vars, baseline.vars, time.var,
                       all.vars(as.formula(exposure.model)),
                       all.vars(as.formula(mediator.model)),
                       all.vars(as.formula(censor.model)),
                       outcome))

  # Fill data, let all subjects has same time.
  # base_dt <- data.table(data)
  # base_dt <- base_dt[, .(max_time = max(get(time.var), na.rm = TRUE)),
  #                    by = c(id.vars, baseline.vars)]
  # base_dt <- unique(base_dt[, c(id.vars, baseline.vars, "max_time"), with = F])
  #
  # dfm <- expand.grid(id = unique(data[[id.vars]]),
  #                    sort(unique(data[[time.var]])))
  # colnames(dfm)[2] <- time.var
  #
  # dfm <- merge(dfm,
  #              base_dt,
  #              by = id.vars, all = TRUE)
  #
  # dfm <- merge(dfm, data, by = c(id.vars, time.var, baseline.vars), all = TRUE)
  #
  # dfm <- data.table::data.table(dfm[, c(all_vars, "max_time")])
  # data.table::setkeyv(dfm, c(id.vars, time.var))
  #
  # # Censor indication
  # dfm[[censor]] <- as.numeric(dfm[[censor]])
  # dfm[[outcome]] <- as.numeric(dfm[[outcome]])
  #
  # dfm[, c(censor, outcome) := list(ifelse(!is.na(get(outcome)), 0,
  #                                         ifelse(max(get(outcome), na.rm = TRUE) == 0 & is.na(get(outcome)),
  #                                                1, NA)),
  #                                  # If outcome then carry on
  #                                  ifelse(is.na(get(outcome)) & max(get(outcome), na.rm = TRUE) == 1,
  #                                         1, get(outcome))),
  #     by = id.vars]
  #
  # # Keep as missing for those before censor or outcome.
  # dfm[, c(censor, outcome) := list(ifelse(get(time.var) < max_time & is.na(get(censor)),
  #                                         NA, get(censor)),
  #                                  ifelse(get(time.var) < max_time & is.na(get(outcome)),
  #                                         NA, get(outcome)))]
  # # Outcome indication
  # dfm <- data.table(as.data.frame(dfm))
  # dfm[order(get(time.var)), cum_Y := get(outcome) + shift(get(outcome), fill = 0), by = id.vars]
  dfm <- data.table(data)

  fit_mods <- lapply(sort(unique(data[[time.var]])), function(t){
    dt <- dfm[get(time.var) == t, ]
    dt <- dt[order(get(id.vars)), ]

    # Weight for mediators
    ip_m1 <- fit_and_predict(data = dt, a_bar = 1,
                             exposure.var = exposure.vars,
                             family = mediator.family, mod = mediator.model)

    ip_m0 <- fit_and_predict(data = dt, a_bar = 0,
                             exposure.var = exposure.vars,
                             family = mediator.family, mod = mediator.model)

    # Weight for exposures
    ip_a0 <- fit_and_predict(data = dt, a_bar = 0,
                             exposure.var = exposure.vars,
                             family = "binomial", mod = exposure.model,
                             is_aform = TRUE)

    ip_c0 <- fit_and_predict(data = dt, a_bar = 0,
                             exposure.var = exposure.vars,
                             family = "binomial", mod = censor.model,
                             is_aform = TRUE)

    ip_a1 <- fit_and_predict(data = dt, a_bar = 1,
                             exposure.var = exposure.vars,
                             family = "binomial", mod = exposure.model,
                             is_aform = TRUE)
    ip_c1 <- fit_and_predict(data = dt, a_bar = 1,
                             exposure.var = exposure.vars,
                             family = "binomial", mod = censor.model,
                             is_aform = TRUE)

    out <- cbind(ip_c0, ip_a0, ip_c1, ip_a1, ip_m0, ip_m1)
    out[dt$cum_Y > 1, ] <- 1
    out
  })

  final_node <- length(unique(data[[time.var]]))

  cum_a0 <- lapply(fit_mods, function(x)x[, 1:2])
  cum_a0 <- do.call("cbind", cum_a0)
  cum_a0 <- bound(t(apply(cum_a0, 1, cumprod)), bound.value)

  cum_a1 <- lapply(fit_mods, function(x)x[, 3:4])
  cum_a1 <- do.call("cbind", cum_a1)
  cum_a1 <- bound(t(apply(cum_a1, 1, cumprod)), bound.value)

  cum_m0 <-  bound(t(apply(sapply(fit_mods, function(x)x[, 5]), 1, cumprod)), bound.value)
  cum_m1 <-  bound(t(apply(sapply(fit_mods, function(x)x[, 6]), 1, cumprod)), bound.value)

  if(is.expL){
    a_0 <- cum_a0[, 2*final_node]
    a_1 <- cum_a1[, 2*final_node]
  }else{
    a_0 <- cum_a0[, 1]
    a_1 <- cum_a1[, 1]
  }
  m_0 <- cum_m0[, final_node]
  m_1 <- cum_m1[, final_node]

  exposure <- all.vars(as.formula(exposure.model))[1]
  censor <- all.vars(as.formula(censor.model))[1]

  dt <- dfm[get(time.var) == max(unique(data[[time.var]])), ]
  dt <- dt[order(get(id.vars)), ]
  index1 <- dt[[exposure]] == 1 & dt[[censor]] == 1
  index1[is.na(index1)] <- FALSE
  index0 <- dt[[exposure]] == 0 & dt[[censor]] == 1
  index0[is.na(index0)] <- FALSE

  phi11 <- (m_1/m_1)/a_1
  phi10 <- (m_0/m_1)/a_1
  phi00 <- (m_0/m_0)/a_0

  phi11 <- weighted.mean(dt[index1, get(outcome)], phi11[index1], na.rm = TRUE)
  phi10 <- weighted.mean(dt[index1, get(outcome)], phi10[index1], na.rm = TRUE)
  phi00 <- weighted.mean(dt[index0, get(outcome)], phi00[index0], na.rm = TRUE)

  nie <- phi11 - phi10
  nde <- phi10 - phi00

  list(call = tpcall,
       weights = list(cum.m0 = cum_m0,
                      cum.m1 = cum_m1,
                      cum.a0 = cum_a0,
                      cum.a1 = cum_a1),
       indexed = list(index0 = index0,
                      inedx1 = index1),
       iptw.med = c("Total effect"    = nie + nde,
                    "Indirect effect" = nie,
                    "Direct effect"   = nde))

}


################################
# fit_and_predict()
################################
#' Calculate predicted value for weights.
#'
#' @param data Data to be used to fit model.
#' @param a_bar Intervention value on exposure.
#' @param exposure.var A vector of exposure variables, including lagged exposure variables.
#' @param family Model link function, \code{gussian},\code{binomial} and \code{multinomial}.
#' For multinomial model, \code{\link[nnet]{multinom}} will be sued for model fitting. And
#' \code{\link[stats]{glm}} for  \code{gussian} and \code{binomial} link funciton model.
#' @param mod A formula string.
#' @param is_aform Indication of model for exposure or censoring, defualt is for mediator model.
#'
#' @return Returns a predicted numerical vector.
#'
#' @import nnet

fit_and_predict <- function(data, a_bar = NULL, exposure.var, family, mod, is_aform = FALSE){

  newdt <- data <- as.data.frame(data)
  # Create NA mtraix for weights
  ip_weights <- matrix(NA_real_, nrow(data))

  if(!is.null(a_bar))
    newdt[, exposure.var] <- a_bar

  if(family  == "binomial"){
    fit <- glm(as.formula(mod), data = data, family = binomial())

    y_vars <- all.vars(fit$formula)[1]

    if(is_aform){
      if(a_bar == 0){
        ip_weights <- 1 - predict.glm(fit, newdata = newdt, type = "response")
      }else{
        ip_weights <- predict.glm(fit, newdata = newdt, type = "response")
      }

    }else{
      indx0 <- newdt[[y_vars]] == 0 & !is.na(newdt[[y_vars]])
      ip_weights[indx0] <- 1 - predict.glm(fit, newdata = newdt, type = "response")[indx0]

      indx1 <- newdt[[y_vars]] == 1 & !is.na(newdt[[y_vars]])
      ip_weights[indx1] <- predict.glm(fit, newdata = newdt, type = "response")[indx1]
    }

  }

  if(family  == "gaussian"){
    fit <- glm(as.formula(mod), data = data)

    y_vars <- all.vars(fit$formula)[1]
    ip_weights <- dnorm(data[[y_vars]], predict(fit, newdata = newdt),
                        sd(fit$residuals))

  }

  if (family == "multinomial"){

    fit <- nnet::multinom(as.formula(mod), data = data)

    y_vars <- all.vars(formula(fit$terms))[1]
    pred <- as.data.frame(nnet::predict(fit, newdata = newdt, type = "probs"))

    for (i in 1:length(unique(newdt[[y_vars]])))
      ip_weights[newdt[[y_vars]] == sort(unique(na.omit(newdt[[y_vars]])))[i]] <-
      pred[newdt[[y_vars]] == sort(unique(na.omit(newdt[[y_vars]])))[i],i]
  }

  return(ip_weights)

}

bound <- function(x, bounds) {
  bounds <- c(bounds, 1 - bounds)
  stopifnot(length(bounds) == 2 && !anyNA(bounds))
  x[x < min(bounds)] <- min(bounds)
  x[x > max(bounds)] <- max(bounds)
  return(x)
}
