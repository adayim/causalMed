
#' Random data simulation from predicted value.
#'
#' @description
#'  Internal use only, predict response and simulate random data. The simulated
#'  value will be restricted within the observed value range for numeric values.
#' 
#' @param models fitted objects defined in the `spec_model`.
#' @param newdt a data frame in which to look for variables with which to predict.
#' 
#' @return A simulated random vector using the predicted value from models and newdt.
#'
#' @export
#'
sim_value <- function(models, newdt){

  var_type <- models$var_type
  fit <- models$fitted
  rmse <- sqrt(mean((stats::fitted(fit) - fit$y)^2, na.rm = TRUE))

  # If the cunstom simulation function is defined.
  if(!is.null(models$custom_sim)){
    out <- models$custom_sim(fit, newdt)

    # Set the simulated values within observed range.
    if(is.numeric(models$val_ran)){
      out <- pmax(out, min(models$val_ran))
      out <- pmin(out, max(models$val_ran))
    }
    return(out)
  }

  # Random number generation for Monte Carlo simulation
  if(var_type == "categorical"){
    pred <- predict(fit, newdata = newdt, type =  "probs")
    rMultinom(pred, 1)

  }else if(var_type == "binomial"){
    pred <- predict(fit, newdata = newdt, type = "response")
    rbinom(nrow(newdt), 1, pred)

  }else{
    pred <- predict(fit, newdata = newdt, type = "response")
    out <- rnorm(nrow(newdt), pred, rmse)

    # Set the simulated values within observed range.
    out <- pmax(out, min(models$val_ran))
    out <- pmin(out, max(models$val_ran))

    return(out)
  }
}







