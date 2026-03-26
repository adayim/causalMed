
causalmed_env <- new.env(parent = emptyenv())

# Suppress R CMD check NOTEs for data.table column names used with := and .SD
utils::globalVariables(c(
  ".",
  "Effect",
  "Est",
  "Estimate",
  "Pred_Y",
  "S",
  "Sc",
  "Sd",
  "new_ID",
  "new_id"
))

#' @importFrom stats coef delete.response family formula model.matrix
#'   na.omit predict quantile rbinom rnorm sd sigma terms
NULL
