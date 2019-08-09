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
#' @param id ID variable per subject.
#' @param trt Intervention/Exposure variable
#' @param med Name of the time-varying mediator.
#' @param y Name of the outcome variable.
#' @param cov A vector of the covariates' name.
#' @param time Time variable of the per row observed. Only for ordering the data.
#' @param m.family link function to be used in the mediator model.
#' @param y.family link function to be used in the outcome model.
#' @param estimator Use IPTW estimator or g-formula estimator or both (Default)
#' @param R The number of bootstrap replicates. Default is 1000.
#'
#' @importFrom boot boot
#'
#' @example medlong(data = dat, trt = "a", med = "mt", y = "status", id = "id", time = "tij", cov = c("w1", "w2"), m.family  = "gaussian", y.family  = "binomial")
#'
#' @export
#'
#'

medlong <- function(data,
                    id,
                    trt,
                    med,
                    y,
                    cov,
                    time,
                    m.family  = c("binomial", "gaussian"),
                    y.family  = c("binomial", "gaussian"),
                    estimator = c("both", "gformula", "iptw"),
                    R         = 1000) {

  data <- format_data(data, id, trt, time, med, y, cov)

  cl <- match.call()
  # args <- mget(names(formals()),sys.frame(sys.nframe()))

  if(estimator[1] %in% c("gformula", "both") & is.null(y.family)){
    stop("Link function of for g-formula ")
  }
  if(estimator[1] == "iptw"){
    args$y.family <- NULL
  }

  # Bootstrap
  if(estimator[1] %in% c("gformula", "both"))
    res.g    <- boot::boot(data, .calc_g,
                           R        = R,
                           trt      = trt,
                           med      = "med_",
                           y        = y,
                           cov      = cov,
                           m.family = m.family,
                           y.family = y.family)


  if(estimator[1] %in% c("iptw", "both"))
    res.iptw <- boot::boot(data, .calc_iptw,
                           R        = R,
                           trt      = trt,
                           med      = "med_",
                           y        = y,
                           cov      = cov,
                           m.family = m.family,
                           y.family = y.family)

  if(estimator[1] == "gformula"){
    out <- list("call"      = cl,
                "g-formula" = extract_boot(res.g, conf.int = T))
  }else if (estimator[1] == "iptw"){
    out <- list("call"      = cl,
                "iptw"      = extract_boot(res.iptw, conf.int = T))
  }else{
    out <- list("call"      = cl,
                "g-formula" = extract_boot(res.g, conf.int = T),
                "iptw"      = extract_boot(res.iptw, conf.int = T))
  }

  out
}


# Calculate IPTW ======================
#' IPTW estimator
#'
#' @description IPTW estimator, internal use.
#'
#' @param data Data set to be sued
#' @param trt Intervention/Exposure variable
#' @param med Prefix of the time-varying mediator.
#' @param y Name of the outcome variable.
#' @param cov A vector of the covariates' name.
#' @param m.family link function to be used in the mediator model.
#' @param y.family link function to be used in the outcome model.


.calc_iptw <- function(data, index, trt, med, y, cov, m.family, y.family = NULL){

  dat <- data[index, ]

  # Make formula, regress on mediator
  make_form <- function(tm){
    if(tm > 1){
      fm <- c(trt, cov, paste0(med, 1:(tm - 1)))
    }else{
      fm <- c(trt, cov)
    }
    fm <- paste(fm, collapse = " + ")
    paste0(paste0(med, tm), " ~ ", fm)
  }

  # Extract maximum time
  tmax <- names(dat)[grep(med, names(dat))]
  tmax <- max(as.numeric(gsub("\\D", "", tmax)))

  # Fit and calculate weight
  w <- lapply(1:tmax, function(t){
    fm <- make_form(tm = t)
    fit <- glm(as.formula(fm), data = dat, family = m.family)

    dat1 <- dat
    dat1$a <- 0
    w1 <- predict(fit , type = "response" , newdata = dat1)

    dat2 <- dat
    dat2$a <- 1
    w2 <- predict(fit , type = "response" , newdata = dat2)

    w3 <- predict(fit , type = "response" , newdata = dat)
    #
    if(all(dplyr::between(w1, 0, 1))){
      cbind(w1, w2, w3)
    }else{
      cbind(dnorm(dat[, paste0(med, t)], w1, sd(fit$residuals, na.rm = TRUE)),
            dnorm(dat[, paste0(med, t)], w2, sd(fit$residuals, na.rm = TRUE)),
            dnorm(dat[, paste0(med, t)], w3, sd(fit$residuals, na.rm = TRUE)))
    }

  })
  w1 <- apply(sapply(w, function(x)x[, 1]), 1, prod, na.rm = TRUE)
  w2 <- apply(sapply(w, function(x)x[, 2]), 1, prod, na.rm = TRUE)
  w3 <- apply(sapply(w, function(x)x[, 3]), 1, prod, na.rm = TRUE)

  fita <- glm(as.formula(paste0(trt, " ~ ", cov)), data = dat, family = "binomial")

  # Natural Indirect effect
  datnew <- dat
  datnew$a <- 1
  wa_nie <- predict(fita , type = "response" , newdata = datnew)
  w <- (1 - w1 / w2) / wa_nie

  # Natural Direct effect
  wa_nde <- predict(fita , type = "response" , newdata = dat)

  nie <- mean(w[dat$a == 1])
  nde <- mean((2*dat[, y] - 1) * (w1 / w3) /wa_nde)

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
#' @param trt Intervention/Exposure variable
#' @param med Prefix of the time-varying mediator.
#' @param y Name of the outcome variable.
#' @param cov A vector of the covariates' name.
#' @param m.family link function to be used in the mediator model.
#' @param y.family link function to be used in the outcome model.


.calc_g <- function(data, index, trt, med, y, cov, m.family, y.family){

  dat <- data[index, ]
  # Make formula
  make_form <- function(tm){
    if(tm > 1){
      fm.m <- c(trt, cov, paste0(med, 1:(tm - 1)))
      fm.y <- c(trt, cov, paste0(med, 1:tm))
    }else{
      fm.m <- c(trt, cov)
      fm.y <- c(trt, cov, paste0(med, 1))
    }

    fm.m <- paste0(paste0(med, tm), " ~ ", paste(fm.m, collapse = " + "))
    fm.y <- paste0(y, " ~ ", paste(fm.y, collapse = " + "))
    list("fm.m" = fm.m, "fm.y" = fm.y)
  }

  # Extract maximum time
  tmax <- names(dat)[grep(med, names(dat))]
  tmax <- max(as.numeric(gsub("\\D", "", tmax)))

  # Fit and calculate weight
  w <- lapply(1:tmax, function(t){
    fm <- make_form(tm = t)
    fit.m <- glm(as.formula(fm$fm.m), data = dat, family = m.family)
    fit.y <- glm(as.formula(fm$fm.y), data = dat, family = y.family)

    dat1 <- dat
    dat1$a <- 1
    dat0 <- dat
    dat0$a <- 0

    p.y1 <- predict(fit.y, type = "response" , newdata = dat1) # phi 1
    p.y0 <- predict(fit.y, type = "response" , newdata = dat0) # phi 0
    if(y.family == "gaussian"){
      p.y1 <- dnorm(dat[, y], p.y1, sd(fit.y$residuals, na.rm = TRUE))
      p.y0 <- dnorm(dat[, y], p.y0, sd(fit.y$residuals, na.rm = TRUE))
    }

    p.m1 <- predict(fit.m, newdata = dat1) # phi 1
    p.m0 <- predict(fit.m, newdata = dat0) # phi 0
    if(m.family == "gaussian"){
      p.m1 <- dnorm(dat[, paste0(med, t)], p.m1, sd(fit.m$residuals, na.rm = TRUE))
      p.m0 <- dnorm(dat[, paste0(med, t)], p.m0, sd(fit.m$residuals, na.rm = TRUE))
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






