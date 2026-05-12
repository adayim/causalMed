
causalmed_env <- new.env(parent = emptyenv())

# Suppress R CMD check NOTEs for data.table column names used with := and .SD
utils::globalVariables(c(
  ".",
  "Effect",
  "Est",
  "Estimate",
  "Intervention",
  "Pred_Y",
  "RD",
  "RR",
  "Risk_type",
  "S",
  "Sc",
  "Sd",
  "Sd_RR",
  "new_ID",
  "new_id",
  "norm_lcl",
  "norm_lcl_RR",
  "norm_ucl",
  "norm_ucl_RR",
  "perct_lcl",
  "perct_lcl_RR",
  "perct_ucl",
  "perct_ucl_RR"
))

#' @importFrom stats coef delete.response family formula model.matrix
#'   na.omit predict quantile rbinom rnorm sd sigma terms
NULL
