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
  
  # 2. Assign a custom class so we can verify this object later
  class(exprs) <- c("causalMed_recodes", "list")
  
  return(exprs)
}

#' @import data.table
apply_recodes <- function(data, recode_params) {
  
  # Handle Custom Class (from recodes())
  if (inherits(recode_params, "causalMed_recodes")) {
    for (i in seq_along(recode_params)) {
      var_name <- names(recode_params)[i]
      expr <- recode_params[[i]]
      
      # Determine if it is a new variable or an update
      # data.table handles both with :=
      data[, (var_name) := eval(expr)]
    }
  } 
  
  invisible(data)
}
