#' Inverse odds ratio weighting.
#'
#' @description Mediation analysis for time fixed mediator using inverse odds ratio.
#' Output contains total effect, natrual direct effect and natural indirect effect.
#'  The confidence interval will be calculated using 1000 bootstrap with normal approximation.
#'
#' @param data Data set to be sued
#' @param trt Intervention/Exposure variable
#' @param med Name of the mediator(s).
#' @param time Survival time in survival analysis.
#' @param ref Only for one categorical mediator, set to NULL for numerical mediator or multiple mediator. Reference value of the mediator, where the mediator is evaluated at its refrence value.
#' @param family a description of the error distribution and link function to be used in the outcome model. "cox" and "aalen" can also be added for survival outcome.
#' @param stabilized Stabilized weights, TRUE(default) or FALSE.
#' @param R The number of bootstrap replicates. Default is 1000.
#'
#' @references
#' Tchetgen Tchetgen, E. J. (2013). Inverse odds ratio‐weighted estimation for causal mediation analysis. \emph{Statistics in medicine}, 32(26), 4567-4580. DOI:10.1002/sim.5864
#'
#' @example iorw(data = dt, trt = "a", med = "mt", y = "status", time = "eventtime", family = "cox", cov = c("w1", "w2"))
#'
#' @importFrom boot boot
#'
#' @export
#'

iorw <- function(fitY,
                 data       = NULL,
                 exposure,
                 mediator,
                 family,
                 ref        = NULL,
                 stabilized = TRUE,
                 R          = 1000) {

  #save input
  tempcall <- match.call()
  #some basic input checks
  if (!("exposure" %in% names(tempcall))) stop("No exposure variable specified")
  if (!("mediator" %in% names(tempcall))) stop("No mediator variable(s) specified")
  if (!("family" %in% names(tempcall)) | ("family" %in% names(tempcall) & !(tempcall$family %in% c("binomial", "gaussian")))) stop("No valid family specified (\"binomial\",  \"gaussian\")")
  if (!is.null(tempcall$data)) message("Data from the model will be used!")
  if (tempcall$mediator %in% all.vars(formula(fitY))) stop("No mediators should be included in the outcome model!")
  if (!tempcall$exposure %in% all.vars(formula(fitY))) stop("Exposure variable should be included in the outcome model!")

  args <- mget(names(formals()),sys.frame(sys.nframe()))

  if(length(args$mediator) > 1 & !is.null(args$ref)){
    warning("Reference values were provided along with multiple mediators. Refrence will be ignored!")
    args$ref <- NULL
  }
  if(length(args$mediator) > 1 & !args$stabilized){
    warning("Multiple mediators were provided, stabilized parameters will be ignored and stabilized weights will be used!")
    args$stabilized <- TRUE
  }

  args$data <- if (missing(data)) {
    extrCall(fitY)$data
  }

  # Model for exposure
  medform <- update(formula(fitY), paste0(exposure, " ~ . + ", paste(mediator, collapse = " + ")))

  args$fitA <- medform
  # Bootstrap
  bootFunc <- function(data, index, ...){
    d <- data[index, ]
    cl <- list(...)
    cl$data <- substitute(d)
    cl$R <- NULL
    coef(do.call(estirow, cl))
  }
  boot.res <- do.call(boot::boot, c(list(statistic = bootFunc), args))

  # Original
  if(!is.null(args$R)){
    args$R <- NULL
  }

  res <- do.call(estirow, args)

  out <- c(call  = match.call(iorw),
           res,
           boots = list(boot.res))
  class(out) <- "iorw"
  return(out)
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
#' @importFrom stats as.formula coef dnorm glm na.omit predict qnorm sd var
#'
estirow <- function(fitA,
                     fitY,
                     data,
                     exposure,
                     mediator,
                     family,
                     ref        = NULL,
                     stabilized = TRUE){

  # Import Data
  tempcall <- match.call()
  # Step 1: fit model for exposure
  Afit <- glm(fitA, data = data, family = family)

  # Step 2: Compute IORW
  # For categorical
  if(!stabilized){
    newdata <- data
    newdata$mediator <- ref
    p1 <- predict(Afit, newdata = data, type="response")
    p2 <- predict(Afit, newdata = newdata, type="response")
    W <- ((1-p1)*p2)/(p1*(1-p2))
    W[data$trt == 0] <- 1

  }else{
    p1 <- predict(Afit, newdata = data, type="response")
    W <- (1 - p1) / p1
  }
  if(family == "gaussian"){
    p1 <- predict(Afit, newdata = data, type="response")
    pt <- predict(Afit, newdata = data, type = "term")
    pre <- rowSums(p1 * pt[mediator], na.rm = TRUE)/sd(Afit$residuals, na.rm = TRUE)
    W <- 1 / exp(pre)
  }

  fitY <- extrCall(fitY)
  fitY$data <- substitute(data)
  # Estimate total effect
  total <- eval(fitY)

  # Estimate direct effect
  fitY$weights <- substitute(W)
  direct <- eval(fitY)

  tot <- coef(total)
  dir <- coef(direct)
  if(is.matrix(tot) | is.data.frame(tot)){
    tot <- tot[, 1]
  }
  if(is.matrix(dir) | is.data.frame(dir)){
    dir <- dir[, 1]
  }

  tot <- tot[grep(exposure, names(tot))]
  dir <- dir[grep(exposure, names(dir))]

  #Natural indirect = total effect - natural direct:
  res <- c(tot, dir, tot - dir)
  names(res) <- c("Total effect","Natural Direct effect","Natural Indirect effect")

  res <- list(Afit    = Afit$call,
              weights = W,
              effect  = res)
  class(res) <- "irow"
  return(res)
}


#' @rdname iorw-methods
#' @export
coef.irow <- function (object, ...)
{
  object$effect
}

#' @rdname iorw-methods
#' @method summary iorw
#' @export
summary.iorw <- function (object, ...){
  coef.table <- extract_boot(object$boots, ...)

  coeftab <- as.matrix(coef.table[, -1])
  dimnames(coeftab) <- list(coef.table$term,
                               c("Estimate", "Bias", "Std.error", "conf.low", "conf.high"))

  prop.med <- coef.table[coef.table$term == "Natural Indirect effect", "statistic"] /
    coef.table[coef.table$term == "Total effect", "statistic"]

  summary <- list(call = object$call, coefficients = coeftab,
                  prop = 100 * prop.med)
  class(summary) <- "summary.iorw"
  attr(summary, "class_object") <- class(object)
  return(summary)
}

#' @method print summary.iorw
#' @export
print.summary.iorw <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)
  cat("------\n")
  cat("Natural effect model\n")
  if (attr(x, "class_object")[1] == "neModelBoot") cat("with standard errors based on the non-parametric bootstrap\n---\n")
  cat("Exposure:", x$call$exposure, "\nMediator(s):", paste(x$call$mediator,
                                                      collapse = ", "), "\n------\n")
  cat("Parameter estimates:\n")
  print(round(x$coefficients, digits = digits))

  cat(paste0("------\nProportion Mediated: ", round(x$prop, digits), "%"))

}













