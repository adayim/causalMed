
#' @method print gmed_boot
#'
#' @export
#'
print.gmed_boot <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Call:\n")
  print(x$call)

  cat("------\n")
  cat("Natural effect model\n")
  cat("with standard errors based on the non-parametric bootstrap\n---\n")
  cat("Parameter estimates:\n")
  out <- cbind(x$estimated, x$confint)
  colnames(out)[1] <- "Estimated"
  print(round(out, digits = digits))

  # cat(paste0("------\nProportion Mediated: ", round(x$prop_med, digits), "%"))
}

#' Print
#'
#' Print method for \code{`gformula`}
#'
#' @param x gformula formula output of \code{\link{gformula}} or \code{\link{mediation}}.
#' @param models logical value indicating if model model details printed.
#' @param ... Not used.
#'
#' @export
print.gformula <- function(x, models = FALSE, ...){

  cat("Call:\n")
  print(x$call)

  # cat(sprintf("\nID variable: %s\n", x$all.args$id_var))
  # cat(sprintf("Baseline variables: %s\n",
  #             paste(x$all.args$base_vars, collapse = ", ")))
  # cat(sprintf("Exposure variable: %s\n", x$all.args$exposure))
  # cat(sprintf("Time variable: %s\n", x$all.args$time_var))

  text1 <- "Estimated effect size of risk ratio/difference"
  cat("\n", text1, "\n")
  row <- paste(rep("=", nchar(text1)), collapse = "")
  cat(row, "\n")

  if(is.null(x$estimate)){
    cat("\nNote: The risk ratio/difference can not be estimated with only 1 intervention.\n")
  }else{
    print(x$estimate)
  }

  text2 <- "Estimated effect"
  cat("\n", text2, "\n")
  row <- paste(rep("=", nchar(text1)), collapse = "")
  cat(row, "\n")
  print(x$effect_size)

  if(models){
    text3 <- "Models specified"
    row <- paste(rep("=", nchar(text3)*2), collapse = "")
    cat(row, "\n")
    cat(text3, "\n")
    cat(row, "\n")
    for(fitmod in x$fitted_models){

      attr_lst <- attributes(fitmod)

      text4 <- sprintf("Model type: %s\nResponse variable type: %s\n",
                       attr_lst$mod_type, attr_lst$var_type)
      row <- paste(rep("-", nchar(text4)), collapse = "")
      cat(row, "\n")
      cat(text4, "\n")
      # cat("Model call:\n")
      if(!is.null(attr_lst$recodes)){
        cat("\nData recodes:\n")
        print(attr_lst$recodes)
      }
      # if(!is.null(attr_lst$subset)){
      #   cat("\nData subset in the fitting process:\n")
      #   print(attr_lst$subset)
      # }

      # cat(row, "\n")
      if(x$all.args$return_fitted){
        cat("Model summary :\n")
        summary(fitmod)
      }else{
        print(fitmod$call)
        cat("Model coefficients :\n")
        print(fitmod$coeff)
      }

    }
  }

  invisible(x)

}
