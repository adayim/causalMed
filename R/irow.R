#' Inverse odds ratio weighting.
#'
#' @description Mediation analysis for time fixed mediator using inverse odds ratio.
#' Output contains total effect, natrual direct effect and natural indirect effect.
#'  The confidence interval will be calculated using 1000 bootstrap with normal approximation.
#'
#' @param data Data set to be sued
#' @param trt Intervention/Exposure variable
#' @param med Name of the mediator(s).
#' @param y Name of the outcome variable.
#' @param cov A vector of the covariates' name.
#' @param time Survival time in survival analysis.
#' @param ref Only for one categorical mediator, set to NULL for numerical mediator or multiple mediator. Reference value of the mediator, where the mediator is evaluated at its refrence value.
#' @param family a description of the error distribution and link function to be used in the outcome model. "cox" and "aalen" can also be added for survival outcome.
#' @param stabilized Stabilized weights, TRUE(default) or FALSE.
#' @param R The number of bootstrap replicates. Default is 1000.
#'
#' @example iorw(data = dt, trt = "a", med = "mt", y = "status", time = "eventtime", family = "cox", cov = c("w1", "w2"))
#'
#' @export
#'

iorw <- function(data,
                 trt,
                 med,
                 y,
                 cov,
                 time       = NULL,
                 family,
                 ref        = NULL,
                 stabilized = TRUE,
                 R          = 1000) {

  cl <- match.call()
  #args <- as.list(match.call(iorw))[-1]
  args <- mget(names(formals()),sys.frame(sys.nframe()))

  if(family %in% c("cox", "aalen") & is.null(time)){
    stop("Time variable name must be provided for survival outcome!")
  }
  if(family %in% c("cox", "aalen")){
    require(survival)
    if(family == "aalen") require(timereg)
  }

  if(length(med) > 1 & !is.null(ref)){
    warning("Reference values were provided along with multiple mediators. Refrence will be ignored!")
    args$ref <- NULL
  }
  if(length(med) > 1 & !stabilized){
    warning("Multiple mediators were provided, stabilized parameters will be ignored and stabilized weights will be used!")
    args$stabilized <- TRUE
  }
  if(length(unique(na.omit(data[, y]))) > 2 )
    stop("Function can only handles the binary outcome")

  # Bootstrap
  boot.res <- do.call(boot::boot, c(list(statistic = .estirow), args))

  out <- list(call      = cl,
              mediation = extract_boot(boot.res, conf.int = T))
  out
}

#' Inverse odds ratio weighting
#'
#' @description  Main function for Inverse odds ratio weighting, internal use.
#'
#' @param data Data set to be sued
#' @param trt Intervention/Exposure variable
#' @param med Name of the mediator(s).
#' @param y Name of the outcome variable.
#' @param cov A vector of the covariates' name.
#' @param time Survival time in survival analysis.
#' @param ref Only for one categorical mediator, set to NULL for numerical mediator or multiple mediator. Reference value of the mediator, where the mediator is evaluated at its refrence value.
#' @param family a description of the error distribution and link function to be used in the outcome model. "cox" and "aalen" can also be added for survival outcome.
#' @param stabilized Stabilized weights, TRUE(default) or FALSE.
#'
#'
.estirow <- function(data, index, trt, med, y, cov, time, ref, family, stabilized){

  d <- data[index, ]
  # Import Data
  dat <- d[, c(trt, med, y, time, cov)]
  colnames(dat) <- c("trt", "mediator", "y", "time", paste0("cov_", 1:length(cov)))
  cov <- paste0("cov_", 1:length(cov))
  form <- paste0("trt ~ mediator + ", paste(cov, collapse = " + "))

  if(length(unique(na.omit(dat$trt))) == 1 | length(unique(na.omit(dat$mediator))) == 1 | length(unique(na.omit(dat$y))) == 1){
    res <- rep(NA, 3)
    warning("Exposure (trt) or Mediator (med) or outcome (y) have only one level! NA will be returned!")
  }else{
    # Step 1: fit model for exposure
    Afit <- glm(as.formula(form), data = dat, family = "binomial")

    # Step 2: Compute IORW
    # For categorical
    if(!stabilized){
      newdata <- dat
      newdata$mediator <- ref
      p1 <- predict(Afit, newdata = dat, type="response")
      p2 <- predict(Afit, newdata = newdata, type="response")
      W <- ((1-p1)*p2)/(p1*(1-p2))
      W[dat$trt == 0] <- 1

    }else{
      newdata <- dat
      p1 <- predict(Afit, newdata = dat, type="response")
      W <- (1 - p1) / p1
    }

    dat$iorw <- W

    if(family == "aalen"){
      cov_c <- paste0("const(", cov, ")")
      fm <- paste0("Surv(time, y) ~ const(trt) + ",
                   paste(cov_c, collapse = " + "))
      # Step 3: estimate direct effect
      dat <- na.omit(dat) # Will be error if there are any missing values
      direct <- aalen(as.formula(fm), data=dat, weights = dat$iorw, robust = 0)

      # Estimate total effect
      total <- aalen(as.formula(fm), data=dat, robust = 0)

      #Natural indirect = total effect - natural direct:
      res <- c(coef.aalen(total)[1,1], coef.aalen(direct)[1,1],
               coef.aalen(total)[1,1] - coef.aalen(direct)[1,1])
    }else if (family == "cox"){

      fm <- paste0("Surv(time, y) ~ trt + ",
                   paste(cov, collapse = " + "))
      direct <- coxph(as.formula(fm), data=dat, weights = dat$iorw)
      # Estimate total effect
      total <- coxph(as.formula(fm), data=dat)

      #Natural indirect = total effect - natural direct:
      res <- c(coef(total)[1], coef(direct)[1],
               coef(total)[1] - coef(direct)[1])
    }else{
      fm <- paste0("y ~ trt + ", paste(cov, collapse = " + "))
      # Step 3: estimate direct effect
      direct <- glm(as.formula(fm), data = dat, weights = dat$iorw, family = family)

      # Estimate total effect
      total <- glm(as.formula(fm), data=dat, family = family)

      tot <- coef(total)[grep("trt", names(coef(total)))]
      dir <- coef(direct)[grep("trt", names(coef(direct)))]
      #Natural indirect = total effect - natural direct:
      res <- c(tot, dir, tot - dir)
    }
  }
  names(res) <- c("Total effect","Natural Direct effect","Natural Indirect effect")
  return(res)
}
