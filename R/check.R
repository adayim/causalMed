
#' Check variables in data
#'
#' Check if the variables in the data, throw an error with variable names if not.
#'
#' @param vars Variables to check
#' @param data Data set to be checked
#'
#' @keywords internal
#'
#' @examples
check_var_in <- function(vars, data){
  diff_vars <- setdiff(vars, names(data))
  if(!identical(diff_vars, character(0)))
    stop("The following variables cannot be found in the data: ",
         paste(diff_vars, collapse = ", "), domain = "causalMed")

}

#' Check intervention settings
#'
#' @param intervention Intervention settings pass from the function.
#' @param time Time index from the data.
#'
#' @keywords internal
#'
#' @examples
check_intervention <- function(intervention, time){
  time_len <- length(unique(time))

  if(!is.list(intervention) & !is.null(names(intervention)))
    stop("Intervention must be a named list object", domain = "causalMed")

  # Intervention should have the same length as time
  if(!is.null(intervention) & all(sapply(tmp, length) == time_len))
    stop("Length of each elements in the intervention should be the same as the `time.var`",
         domain = "causalMed")

}





