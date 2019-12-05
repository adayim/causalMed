
iptw_med <- function(data,
                     id.vars,
                     time.var,
                     baseline.vars,
                     exposure.model,
                     exposure.vars,
                     mediator.model,
                     mediator.family,
                     outcome,
                     is.expL      = TRUE,
                     censor.model = NULL,
                     bound.value  = 0,
                     R            = 500){

  tpcall <- match.call()

  args <- mget(names(formals()),sys.frame(sys.nframe()))

  if(length(bound.value) > 1)
    stop("bound.value should be length of one!")
  if(bound.value <0 | bound.value > 0.5)
    stop("bound.value should be between 0 and 0.5!")

  id.vars <- substitute(id.vars)
  time.var <- substitute(time.var)
  outcome <- substitute(outcome)

  dfm <- data.table(data)
  dfm[order(get(time.var)), cum_Y := get(outcome) + shift(get(outcome), fill = 0), by = id.vars]

  splited <- split(dfm, f = dfm[[id.vars]])

  sample_id <- unique(dfm[[id.vars]])

  args$data <- substitute(sample_id)

  bootFunc <- function(data, index, ...){
    smp_id <- data[index]
    samp_df <- lapply(seq_along(smp_id), function(x){
      res <- splited[[as.character(smp_id[x])]]
      res$newID <- x
      return(res)
    })
    samp_df <- do.call(rbind, samp_df)

    cl <- list(...)
    cl$data <- substitute(samp_df)
    cl$id.vars <- "newID"
    cl$R <- NULL
    res <- do.call(calc_iptw, cl)
    res$iptw.med
  }
  boot.res <- do.call(boot::boot, c(list(statistic = bootFunc), args))

  cls <- tpcall
  cls[[1]] <- substitute(calc_iptw)
  cls$data <- substitute(dfm)
  cls$R <- NULL

  out <- list(call = tpcall,
              boot.res = list(boot.res),
              confint = extract_boot(boot.res))

  out <- append(out, res)
  class(out) <- "iptw_med"
  return(out)
}


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

  id.vars <- substitute(id.vars)
  time.var <- substitute(time.var)
  outcome <- substitute(outcome)

  dfm <- data.table(data)
  dfm <- dfm[order(get(id.vars)), ]

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

  calc_prod <- function(obj, pos){
    res <- lapply(obj, function(x)x[, pos])
    res <- do.call("cbind", res)
    bound(t(apply(res, 1, cumprod)), bound.value)
  }

  cum_a0 <- calc_prod(fit_mods, 1:2)
  cum_a1 <- calc_prod(fit_mods, 3:4)
  cum_m0 <-  calc_prod(fit_mods, 5)
  cum_m1 <-  calc_prod(fit_mods, 6)

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

  #Returns the intervention match for the closest A node.
  #It will carry forward from the last A point if not available at the next
  interv1 <- calc_match(dfm, 1, "A")
  interv0 <- calc_match(dfm, 0, "A")
  #Check censoring for the closest
  uncens <- calc_match(dfm, 1, names(dfm)[grep('C_',names(dfm))], is_censor = TRUE)
  indx1 <- interv1 &  uncens[, max_time]
  indx0 <- interv0 &  uncens[, max_time]

  phi11 <- 1/a_1
  phi10 <- (m_0/m_1)/a_1
  phi00 <- 1/a_0

  phi11 <- weighted.mean(dt[index1, get(outcome)], phi11[index1], na.rm = TRUE)
  phi10 <- weighted.mean(dt[index1, get(outcome)], phi10[index1], na.rm = TRUE)
  phi00 <- weighted.mean(dt[index0, get(outcome)], phi00[index0], na.rm = TRUE)

  nie <- phi11 - phi10
  nde <- phi10 - phi00

  list(weights = list(cum.m0 = cum_m0,
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
# calc_match()
################################

#' calc_match
#'
#' Determine which patients are uncensored or treatment equals to some value.
#' The none-missing value will be carried forward.
#'
#' @param data Dataset.
#' @param vals Values to be used to determine.
#' @param nodes Varaibles to be used to determine the value.
#' @param is_censor Is the nodes is censoring nodes.
#'
#' @return Returns vector of [numDataRows x 1] I(nodes = vals) from
#' nodes[1] to the nodes just before pernodes.
#'

calc_match <- function(data, vals, nodes, is_censor = FALSE){
  calc <- matrix(nrow=nrow(data), ncol=length(nodes))
  cum_calc <- rep(TRUE, nrow(data))
  
  if(is_censor & length(nodes) == 1){
    return(rep(TRUE, nrow(data)))
  }else{
    if(length(nodes) == 1){
      res <- data[[nodes]] == vals
      res[is.na(res)] <- FALSE
      return(res)
    }else{
      df <- t(apply(data[, nodes], 1, function(x){
        v <- !is.na(x)
        c(NA, x[v])[cumsum(v)+1]
      }))
      return(df == vals)
    }
  }
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

#' @method print print.iptw_med
#'
#' @export
#'
print.iptw_med <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)

  cat("------\n")
  cat("with standard errors based on the non-parametric bootstrap\n---\n")
  cat("Parameter estimates:\n")
  print(round(x$confint, digits = digits))
}
