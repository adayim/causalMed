#' Inverse odds ratio weighting.
#'
#' @description Mediation analysis for time fixed mediator using inverse odds ratio.
#' Output contains total effect, natrual direct effect and natural indirect effect.
#'  The confidence interval will be calculated using 1000 bootstrap with normal approximation.
#'
#' @param fitY model object of the final outcome, all variables of interest should be included except mediators.
#' @param data Data set to be sued, if NULL, the data from outcome model will be used.
#' @param exposure Intervention/Exposure variable, control or unexposure must be set to 0.
#' @param mediator Name of the mediator(s).
#' @param ref Only for one categorical mediator, set to NULL for numerical mediator or
#'  multiple mediator. Reference value of the mediator, where the mediator is
#'  evaluated at its refrence value.
#' @param family a description of the error distribution and link function to be used in the exposure model. Only binomial, multinomial and gaussian supported now.
#' @param stabilized Stabilized weights, TRUE(default) or FALSE.
#' @param R The number of bootstrap replicates. Default is 1000.
#'
#' @references
#' Tchetgen Tchetgen, E. J. (2013). Inverse odds ratio‐weighted estimation for causal mediation analysis. \emph{Statistics in medicine}, 32(26), 4567-4580. DOI:10.1002/sim.5864
#'
#' @examples
#'
#' data(lipdat)
#' dtbase <- lipdat[lipdat$time == 0, ]
#' out <- iorw(coxph(Surv(os, cvd) ~ bmi + age0 + smoke, data = dtbase),
#' exposure   = "smoke",
#' mediator   = "hdl",
#' family     = "binomial")
#'
#' summary(out)
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

  # Setting seeds
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } ) # Set to it's roginal seed after exit
  set.seed(2)

  #save input
  tempcall <- match.call()
  #some basic input checks
  if (!("exposure" %in% names(tempcall))) stop("No exposure variable specified")
  if (!("mediator" %in% names(tempcall))) stop("No mediator variable(s) specified")

  dt_name <- extrCall(fitY)$data
  if (deparse(substitute(dt_name)) != deparse(substitute(data)) & !is.null(tempcall$data))
    stop("Data set name of the outcome model and the specified data are not the same! Make sure correct dataset is used.")

  if (!("family" %in% names(tempcall)) | ("family" %in% names(tempcall) & !(tempcall$family %in% c("binomial", "gaussian", "multinomial")))) stop("No valid family specified (\"binomial\",  \"gaussian\",  \"multinomial\")")
  if (!is.null(tempcall$data)) message("Data from the model will be used!")
  if (any(mediator %in% all.vars(formula(fitY)))) stop("No mediators should be included in the outcome model!")
  if (!exposure %in% all.vars(formula(fitY))) stop("Exposure variable should be included in the outcome model!")

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

  if (tempcall$family %in% c("binomial", "multinomial") & !0 %in% unique(eval(args$data)[, exposure])) stop("Exposure must include level 0 and set as control/unexposure")

  # Model for exposure
  medform <- update(formula(fitY), paste0(exposure, " ~ . + ", paste(mediator, collapse = " + ")))

  args$fitA <- substitute(medform)
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
#' @param fitA Exposure model formula, glm will be used in this step.
#' @param fitY model object of the final outcome, all variables of interest should be included except mediators.
#' @param data Data set to be sued, if NULL, the data from outcome model will be used.
#' @param exposure Intervention/Exposure variable
#' @param mediator Name of the mediator(s).
#' @param ref Only for one categorical mediator, set to NULL for numerical mediator or multiple mediator. Reference value of the mediator, where the mediator is evaluated at its refrence value.
#' @param family a description of the error distribution and link function to be used in the exposure model. Only binomial and gaussian supported now.
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

  #weights binomial
  if(family == "binomial"){
    Afit <- do.call("glm", eval(list(formula = fitA, data = substitute(data), family = family)))
    # Step 2: Compute IORW
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
  }

  #weights multinomial
  if(family == "multinomial"){
    tempdat <- data.frame(exposure = data[,as.character(exposure)])

    Afit <- do.call("nnet::multinom", eval(list(formula = fitA, data = substitute(data))))
    # Step 2: Compute IORW
    if(!stabilized){
      newdata <- data
      newdata$mediator <- ref
      p1 <- as.data.frame(predict(Afit, newdata = data, type="probs"))
      for (i in 1:length(unique(tmpdat$exposure))){
        tmpdat$p1[with(tmpdat, exposure == sort(unique(tmpdat$exposure))[i])] <-
          p1[tmpdat$exposure == sort(unique(tmpdat$exposure))[i],i]
      }
      p2 <- as.data.frame(predict(Afit, newdata = newdata, type="probs"))
      for (i in 1:length(unique(tmpdat$exposure))){
        tmpdat$p2[with(tmpdat, exposure == sort(unique(tmpdat$exposure))[i])] <-
          p2[tmpdat$exposure == sort(unique(tmpdat$exposure))[i],i]
      }

      W <- ((1-tmpdat$p1)*tmpdat$p2)/(tmpdat$p1*(1-tmpdat$p2))
      W[data$trt == 0] <- 1

    }else{
      p1 <- as.data.frame(predict(Afit, newdata = data, type="probs"))
      for (i in 1:length(unique(tmpdat$exposure))){
        tmpdat$p1[with(tmpdat, exposure == sort(unique(tmpdat$exposure))[i])] <-
          p1[tmpdat$exposure == sort(unique(tmpdat$exposure))[i],i]
      }
      W <- (1 - tmpdat$p1) / tmpdat$p1
    }
  }

  # weights gaussian
  if(family == "gaussian"){
    Afit <- do.call("glm", eval(list(formula = fitA, data = substitute(data), family = family)))

    #p1 <- predict(Afit, newdata = data, type="response")
    pt <- predict(Afit, newdata = data, type = "term")
    pre <- rowSums(data[,as.character(exposure)] * pt[, mediator], na.rm = TRUE)/sd(Afit$residuals, na.rm = TRUE)
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
              Yfit    = extrCall(tempcall$fitY),
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

  summary <- list(call         = object$call,
                  Acall        = object$Afit,
                  Ycall        = object$Yfit,
                  coefficients = coeftab,
                  prop         = 100 * prop.med)

  class(summary) <- "summary.iorw"
  attr(summary, "class_object") <- class(object)
  return(summary)
}

#' @method print summary.iorw
#' @export
print.summary.iorw <- function (x, digits = max(3, getOption("digits") - 3), ...){
  cat("Call:\n")
  print(x$call)

  cat("\nOutcome Model Call:\n")
  print(x$Ycall)

  cat("\nExposure Model Call:\n")
  print(x$Acall)

  cat("------\n")
  cat("Natural effect model\n")
  cat("with standard errors based on the non-parametric bootstrap\n---\n")
  cat("Exposure:", x$call$exposure, "\nMediator(s):", paste(x$call$mediator,
                                                      collapse = ", "), "\n------\n")
  cat("Parameter estimates:\n")
  print(round(x$coefficients, digits = digits))

  cat(paste0("------\nProportion Mediated: ", round(x$prop, digits), "%"))

}














