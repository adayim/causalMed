
# Internal function to intialise the warning
init_warn <- function(){
  causalmed_env$warning <- c()
}

# Collect warnings
run_withwarning_collect <- function(expr, msg) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      warn_messg <- sprintf("%s;\nWarninig message:\n%s", 
                            msg, 
                            paste(w$message, collapse = "\n"))
      causalmed_env$warning <- c(causalmed_env$warning, warn_messg)
    }                  
  )
}

