#' Mediation analysis for time-varying mediaiton, survival and non-survival
#'
#' @description Mediation analysis for time varying mediator, estimation based
#'  on g-formula and IPTW. Output contains total effect,
#' natrual direct effect and natural indirect effect. The confidence interval will
#' be calculated using 1000 bootstrap with normal approximation.
#'
#' @note Not that final outcome must be the same for in all rows per subject.
#'   If the dataset is  survival settings, records after the interested outcome
#'    must be deleted. The funciton it self do some data manupilation internally.
#'    Please prepare the data as longitudinal format.
#'
#' @param data Data set to be sued
#' @param id.var ID variable per subject.
#' @param exposure Intervention/Exposure variable
#' @param mediator Name of the time-varying mediator.
#' @param outcome Name of the outcome variable.
#' @param covariates A vector of the covariates' name.
#' @param time.var Time variable of the per row observed. Only for ordering the data.
#' @param m.family link function to be used in the mediator model.
#' @param y.family link function to be used in the outcome model.
#' @param estimator Use IPTW estimator or g-formula estimator or both (Default)
#' @param R The number of bootstrap replicates. Default is 1000.
#'
#' @references
#' Lin, S. H., Young, J. G., Logan, R., & VanderWeele, T. J. (2017). Mediation analysis for a survival outcome with time‐varying exposures, mediators, and confounders. \emph{Statistics in medicine}, 36(26), 4153-4166. DOI:10.1002/sim.7426
#' Zheng, W., & van der Laan, M. (2017). Longitudinal mediation analysis with time-varying mediators and exposures, with application to survival outcomes. \emph{Journal of causal inference}, 5(2). DOI:10.1515/jci-2016-0006
#'
#' @importFrom boot boot
#'
#' @examples
#' data(lipdat)
#' res <- medlong(data       = lipdat,
#' exposure   = "smoke",
#' mediator   = "hdl",
#' outcome    = "cvd",
#' id.var     = "id",
#' time.var   = "time",
#' covariates = c("bmi", "gender", "age0"),
#' m.family   = "gaussian",
#' y.family   = "binomial")
#'
#' summary(res)
#'
#'
#' @export
#'
#'

medlong <- function(data,
                    id.var,
                    exposure,
                    mediator,
                    outcome,
                    covariates,
                    time.var,
                    m.family  = c("binomial", "gaussian"),
                    y.family  = c("binomial", "gaussian"),
                    estimator = c("both", "gformula", "iptw"),
                    R         = 1000) {

  cl <- match.call()

  if (!("m.family" %in% names(cl)) | ("m.family" %in% names(cl) & !(cl$m.family %in% c("binomial", "gaussian")))) stop("No valid family specified for mediator (\"binomial\",  \"gaussian\")")
  if (!("y.family" %in% names(cl)) | ("y.family" %in% names(cl) & !(cl$y.family %in% c("binomial", "gaussian")))) stop("No valid family specified for outcome (\"binomial\",  \"gaussian\")")

  if(estimator[1] %in% c("gformula", "both") & is.null(y.family)){
    stop("Link function of for g-formula ")
  }
  if(estimator[1] == "iptw"){
    y.family <- NULL
  }

  # args <- mget(names(formals()),sys.frame(sys.nframe()))
  # Reformat Data
  data <- data[, c(id.var, exposure, outcome, mediator, time.var, covariates)]

  # Get the order number of observation per subject
  dfm <- transform(data, rank_ord = ave(1:nrow(data), eval(as.name(id.var)),
                                        FUN = function(x)
                                          order(eval(parse(text = paste0(time.var, "[x]"))))))

  data <- reshape(dfm,
                  drop = time.var,
                  idvar =  c(id.var, exposure, outcome, covariates),
                  v.names = mediator,
                  sep = "_",
                  timevar = "rank_ord", direction = "wide")


  # Bootstrap
  if(estimator[1] %in% c("gformula", "both"))
    res.g    <- boot::boot(data, calc_g,
                           R          = R,
                           exposure   = exposure,
                           mediator   = paste0(mediator, "_"),
                           outcome    = outcome,
                           covariates = covariates,
                           m.family   = m.family,
                           y.family   = y.family)


  if(estimator[1] %in% c("iptw", "both"))
    res.iptw <- boot::boot(data, calc_iptw,
                           R          = R,
                           exposure   = exposure,
                           mediator   = paste0(mediator, "_"),
                           outcome    = outcome,
                           covariates = covariates,
                           m.family   = m.family,
                           y.family   = y.family)

  if(estimator[1] == "gformula"){
    out <- list("call"      = cl,
                "g-formula" = res.g)
  }else if (estimator[1] == "iptw"){
    out <- list("call"      = cl,
                "iptw"      = res.iptw)
  }else{
    out <- list("call"      = cl,
                "g-formula" = res.g,
                "iptw"      = res.iptw)
  }

  class(out) <- "medlong"
  return(out)
}


# Calculate IPTW ======================
#' IPTW estimator
#'
#' @description IPTW estimator, internal use.
#'
#' @param data Data set to be sued
#' @param exposure Intervention/Exposure variable
#' @param mediator Prefix of the time-varying mediator.
#' @param outcome Name of the outcome variable.
#' @param covariates A vector of the covariates' name.
#' @param m.family link function to be used in the mediator model.
#' @param y.family link function to be used in the outcome model.


calc_iptw <- function(data, index, exposure, mediator, outcome,
                      covariates, m.family, y.family = NULL){

  dat <- data[index, ]

  # Make formula, regress on mediator
  make_form <- function(tm){
    if(tm > 1){
      fm <- c(exposure, covariates, paste0(mediator, 1:(tm - 1)))
    }else{
      fm <- c(exposure, covariates)
    }
    fm <- paste(fm, collapse = " + ")
    paste0(paste0(mediator, tm), " ~ ", fm)
  }

  # Extract maximum time
  tmax <- names(dat)[grep(mediator, names(dat))]
  tmax <- max(as.numeric(gsub(mediator, "", tmax)))

  # Fit and calculate weight
  w <- lapply(1:tmax, function(t){
    fm <- make_form(tm = t)
    fit <- glm(as.formula(fm), data = dat, family = m.family)

    dat1 <- dat
    dat1[, exposure] <- 0
    w1 <- predict(fit , type = "response" , newdata = dat1)

    dat2 <- dat
    dat2[, exposure] <- 1
    w2 <- predict(fit , type = "response" , newdata = dat2)

    w3 <- predict(fit , type = "response" , newdata = dat)

    #
    if(m.family == "binomial"){
      return(cbind(w1, w2, w3))
    }else{
      return(cbind(dnorm(dat[, paste0(mediator, t)], w1, sd(fit$residuals, na.rm = TRUE)),
                   dnorm(dat[, paste0(mediator, t)], w2, sd(fit$residuals, na.rm = TRUE)),
                   dnorm(dat[, paste0(mediator, t)], w3, sd(fit$residuals, na.rm = TRUE))))
    }
  })

  w1 <- apply(sapply(w, function(x)x[, 1]), 1, prod, na.rm = TRUE)
  w2 <- apply(sapply(w, function(x)x[, 2]), 1, prod, na.rm = TRUE)
  w3 <- apply(sapply(w, function(x)x[, 3]), 1, prod, na.rm = TRUE)

  fita <- glm(as.formula(paste0(exposure, " ~ ", covariates)), data = dat, family = "binomial")

  # Natural Indirect effect
  datnew <- dat
  datnew[, exposure] <- 1
  wa_nie <- predict(fita , type = "response" , newdata = datnew)
  w <- (1 - w1 / w2) / wa_nie

  # Natural Direct effect
  wa_nde <- predict(fita , type = "response" , newdata = dat)

  nie <- mean(w[dat[, exposure] == 1])
  nde <- mean((2*dat[, exposure] - 1) * (w1 / w3) /wa_nde)

  c("Total effect"            = nie + nde,
    "Natural Indirect effect" = nie,
    "Natural Direct effect"   = nde)
}



# Calculate g-formula  =================
#' G-formula estimator
#'
#' @description G-formula estimator, internal use.
#'
#' @param data Data set to be sued
#' @param exposure Intervention/Exposure variable
#' @param mediator Prefix of the time-varying mediator.
#' @param outcome Name of the outcome variable.
#' @param covariates A vector of the covariates' name.
#' @param m.family link function to be used in the mediator model.
#' @param y.family link function to be used in the outcome model.


calc_g <- function(data, index, exposure, mediator,
                   outcome, covariates, m.family, y.family){

  dat <- data[index, ]
  # Make formula
  make_form <- function(tm){
    if(tm > 1){
      fm.m <- c(exposure, covariates, paste0(mediator, 1:(tm - 1)))
      fm.y <- c(exposure, covariates, paste0(mediator, 1:tm))
    }else{
      fm.m <- c(exposure, covariates)
      fm.y <- c(exposure, covariates, paste0(mediator, 1))
    }

    fm.m <- paste0(paste0(mediator, tm), " ~ ", paste(fm.m, collapse = " + "))
    fm.y <- paste0(outcome, " ~ ", paste(fm.y, collapse = " + "))
    list("fm.m" = fm.m, "fm.y" = fm.y)
  }

  # Extract maximum time
  tmax <- names(dat)[grep(mediator, names(dat))]
  tmax <- max(as.numeric(gsub(mediator, "", tmax)))

  # Fit and calculate weight
  w <- lapply(1:tmax, function(t){
    fm <- make_form(tm = t)
    fit.m <- glm(as.formula(fm$fm.m), data = dat, family = m.family)
    fit.y <- glm(as.formula(fm$fm.y), data = dat, family = y.family)

    dat1 <- dat
    dat1[, exposure] <- 1
    dat0 <- dat
    dat0[, exposure] <- 0

    p.y1 <- predict(fit.y, type = "response" , newdata = dat1) # phi 1
    p.y0 <- predict(fit.y, type = "response" , newdata = dat0) # phi 0
    if(y.family == "gaussian"){
      p.y1 <- dnorm(dat[, outcome], p.y1, sd(fit.y$residuals, na.rm = TRUE))
      p.y0 <- dnorm(dat[, outcome], p.y0, sd(fit.y$residuals, na.rm = TRUE))
    }

    p.m1 <- predict(fit.m, newdata = dat1) # phi 1
    p.m0 <- predict(fit.m, newdata = dat0) # phi 0
    if(m.family == "gaussian"){
      p.m1 <- dnorm(dat[, paste0(mediator, t)], p.m1, sd(fit.m$residuals, na.rm = TRUE))
      p.m0 <- dnorm(dat[, paste0(mediator, t)], p.m0, sd(fit.m$residuals, na.rm = TRUE))
    }

    phi_11 <- apply(cbind(p.y1, p.m1), 1, prod, na.rm = TRUE)
    phi_10 <- apply(cbind(p.y1, p.m0), 1, prod, na.rm = TRUE)
    phi_00 <- apply(cbind(p.y0, p.m0), 1, prod, na.rm = TRUE)

    cbind(phi_11, phi_10, phi_00)
  })

  phi_11 <- apply(sapply(w, function(x)x[, 1]), 1, prod, na.rm = TRUE)
  phi_10 <- apply(sapply(w, function(x)x[, 2]), 1, prod, na.rm = TRUE)
  phi_00 <- apply(sapply(w, function(x)x[, 3]), 1, prod, na.rm = TRUE)

  nie <- sum(phi_11) - sum(phi_10)
  nde <- sum(phi_10) - sum(phi_00)

  c("Total effect"            = nie + nde,
    "Natural Indirect effect" = nie,
    "Natural Direct effect"   = nde)
}


#' @rdname medlong-methods
#' @method summary medlong
#' @export
summary.medlong <- function (object, ...){

  nams <- names(object)

  # Extract bootstrap results
  extrcF <- function(x, ..){
    coef.table <- extract_boot(x, ...)
    coeftab <- as.matrix(coef.table[, -1])
    dimnames(coeftab) <- list(coef.table$term,
                              c("Estimate", "Bias", "Std.error", "conf.low", "conf.high"))
    prop.med <- coef.table[coef.table$term == "Natural Indirect effect", "statistic"] /
      coef.table[coef.table$term == "Total effect", "statistic"]
    list(coeftab = coeftab, prop.med = prop.med)
  }

  if("iptw" %in% nams){
    iptw.res <- extrcF(object$iptw, ...)
  }else{
    iptw.res <- list(coeftab = NULL, prop.med = NULL)
  }

  if("g-formula" %in% nams){
    gform.res <- extrcF(object$`g-formula`, ...)
  }else{
    gform.res <- list(coeftab = NULL, prop.med = NULL)
  }

  summary <- list(call         = object$call,
                  coeff.iptw   = iptw.res$coeftab,
                  coeff.gform  = gform.res$coeftab,
                  prop.iptw    = 100 * iptw.res$prop.med,
                  prop.gform   = 100 * gform.res$prop.med)

  class(summary) <- "summary.medlong"
  attr(summary, "class_object") <- class(object)
  return(summary)
}

#' @method print summary.medlong
#' @export
print.summary.medlong <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)

  cat("------\n")
  cat("Exposure:", x$call$exposure, "\nMediator:", paste(x$call$mediator,
                                                            collapse = ", "))

  cat("Estimation of standard errors based on the non-parametric bootstrap\n---\n")
  if (!is.null(x$coeff.gform)){
    cat("\n------\n")
    cat("Natural effect parameter estimates with G-formula estimation:\n")
    print(round(x$coeff.gform, digits = digits))
    cat(paste0("------\nProportion Mediated: ", round(x$prop.gform, digits), "%"))
  }

  if (!is.null(x$coeff.iptw)){
    cat("\n------\n")
    cat("\nNatural effect parameter estimates with IPTW estimation:\n")
    print(round(x$coeff.iptw, digits = digits))
    cat(paste0("------\nProportion Mediated: ", round(x$prop.iptw, digits), "%"))
  }
}





