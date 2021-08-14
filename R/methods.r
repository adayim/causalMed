
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
