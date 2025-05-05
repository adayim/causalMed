# Derive parameters needed for the function from current environment
# Ref: https://stackoverflow.com/a/51002887
get_args_for <- function(fun, env = parent.frame(), ..., dots = NULL) {
  potential <- names(formals(fun))

  if ("..." %in% potential) {
    if (missing(dots)) {
      # return everything from parent frame
      return(as.list(env))
    } else if (!is.list(dots)) {
      stop("If provided, 'dots' should be a list.")
    }

    potential <- setdiff(potential, "...")
  }

  # get all formal arguments that can be found in parent frame
  args <- mget(potential, env, ..., ifnotfound = list(NULL), inherits = FALSE)
  # remove not found
  args <- args[sapply(args, Negate(is.null))]
  # return found args and dots
  c(args, dots)
}

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

