#' Define parameters for recoding
#' 
#' @description Capture expressions for variable assignment.
#' @param ... Named expressions (e.g., daysq = day^2).
#' @return A list of expressions with class 'causalMed_recodes'.
#' @export
recodes <- function(...) {
  # 1. Capture arguments as unevaluated expressions
  # alist() ensures we get 'day^2' without trying to run it
  exprs <- eval(substitute(alist(...)))

  class(exprs) <- c("causalMed_recodes", "list")
  
  return(exprs)
}

#' @import data.table
apply_recodes <- function(data, recode_params) {
  
  # Handle Custom Class (from recodes())
  for (i in seq_along(recode_params)) {
    var_name <- names(recode_params)[i]
    expr <- recode_params[[i]]
    
    # Evaluate expression within the data.table environment
    # 'var_name' is wrapped in parentheses to use the string as the column name
    data[, (var_name) := eval(expr)]
  }
  
  invisible(data)
}
