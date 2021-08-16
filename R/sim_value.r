
#' Random data simulation from predicted value.
#'
#' @description
#'  Internal use only, predict response and simulate random data. The simulated
#'  value will be restricted within the observed value range for numeric values.
#'
#' @param model fitted objects defined in the `spec_model`.
#' @param newdt a data frame in which to look for variables with which to predict.
#'
#' @return A simulated random vector using the predicted value from model and newdt.
#' 
#'
#' @export
#'
sim_value <- function(model, newdt) {
  var_type <- model$var_type
  fit <- model$fitted
  rmse <- sqrt(mean((stats::fitted(fit) - fit$y)^2, na.rm = TRUE))

  # If the cunstom simulation function is defined.
  if (!is.null(model$custom_sim)) {
    out <- model$custom_sim(fit, newdt)

    # Set the simulated values within observed range.
    if (is.numeric(model$val_ran)) {
      out <- pmax(out, min(model$val_ran))
      out <- pmin(out, max(model$val_ran))
    }
    return(out)
  }

  # Random number generation for Monte Carlo simulation
  if (var_type == "categorical") {
    pred <- predict(fit, newdata = newdt, type = "probs")
    rMultinom(pred, 1)
  } else if (var_type == "binary") {
    pred <- predict(fit, newdata = newdt, type = "response")
    rbinom(nrow(newdt), 1, pred)
  } else {
    pred <- predict(fit, newdata = newdt, type = "response")
    out <- rnorm(nrow(newdt), pred, rmse)

    # Set the simulated values within observed range.
    out <- pmax(out, min(model$val_ran))
    out <- pmin(out, max(model$val_ran))

    return(out)
  }
}
