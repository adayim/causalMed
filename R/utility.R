

## Code is derived from Hmisc package
rMultinom <- function(probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) {
    lev <- 1:k
  }

  ran <- matrix(lev[1], ncol = m, nrow = n)
  z <- apply(probs, 1, sum)
  if (any(abs(z - 1) > .00001)) {
    stop("error in multinom: probabilities do not sum to 1")
  }

  U <- apply(probs, 1, cumsum)
  for (i in 1:m)
  {
    un <- rep(runif(n), rep(k, n))
    ran[, i] <- lev[1 + apply(un > U, 2, sum)]
  }

  ran
}

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


# Extract data
extrCall <- function(x) {
  if (isS4(x)) x@call else x$call
}


#' Extract bootstrap Confidence interval
extract_boot <- function(x,
                         conf.int = TRUE,
                         conf.level = 0.95,
                         conf.method = c("perc", "bca", "basic", "norm"),
                         ...) {
  conf.method <- match.arg(conf.method)

  # calculate the bias and standard error
  # this is an adapted version of the code in print.boot, where the bias
  # and standard error are calculated
  boot.out <- x
  index <- 1:ncol(boot.out$t)
  sim <- boot.out$sim
  cl <- boot.out$call
  t <- matrix(boot.out$t[, index], nrow = nrow(boot.out$t))
  allNA <- apply(t, 2L, function(t) all(is.na(t)))
  index <- index[!allNA]
  t <- matrix(t[, !allNA], nrow = nrow(t))
  rn <- paste("t", index, "*", sep = "")

  if (is.null(t0 <- boot.out$t0)) {
    if (is.null(boot.out$call$weights)) {
      op <- cbind(
        apply(t, 2L, mean, na.rm = TRUE),
        sqrt(apply(t, 2L, function(t.st) var(t.st[!is.na(t.st)])))
      )
    } else {
      op <- NULL
      for (i in index) op <- rbind(op, boot::imp.moments(boot.out, index = i)$rat)
      op[, 2L] <- sqrt(op[, 2])
    }
    colnames(op) <- c("estimate", "std.error")
  } else {
    t0 <- boot.out$t0[index]
    if (is.null(boot.out$call$weights)) {
      op <- cbind(t0, apply(t, 2L, mean, na.rm = TRUE) -
        t0, sqrt(apply(t, 2L, function(t.st) var(t.st[!is.na(t.st)]))))
      colnames(op) <- c("statistic", "bias", "std.error")
    } else {
      op <- NULL
      for (i in index) {
        op <- rbind(op, boot::imp.moments(boot.out,
          index = i
        )$rat)
      }
      op <- cbind(t0, op[, 1L] - t0, sqrt(op[, 2L]), apply(t,
        2L, mean,
        na.rm = TRUE
      ))
      colnames(op) <- c("statistic", "bias", "std.error", "estimate")
    }
  }

  # bring in rownames as "term" column, and turn into a data.frame
  ret <- cbind.data.frame(op, stringsAsFactors = FALSE)

  if (conf.int) {
    ci.list <- lapply(seq_along(x$t0),
      boot::boot.ci,
      boot.out = x,
      conf = conf.level, type = conf.method
    )

    ## boot.ci uses c("norm", "basic", "perc", "stud") for types
    ## stores them with longer names
    ci.pos <- pmatch(conf.method, names(ci.list[[1]]))

    if (conf.method == "norm") {
      ci.tab <- cbind(ci.list[[1]][ci.pos][[1]][2:3], ci.list[[2]][ci.pos][[1]][2:3])
    } else {
      ci.tab <- t(sapply(ci.list, function(x) x[[ci.pos]][4:5]))
    }

    colnames(ci.tab) <- c("conf.low", "conf.high")
    ret <- cbind(ret, ci.tab)
  }
  ret
}

bound <- function(obj, trunc = 0.01) {
  if (!is.matrix(obj)) {
    quant1 <- quantile(obj, 0 + trunc, na.rm = TRUE)
    quant2 <- quantile(obj, 1 - trunc, na.rm = TRUE)
    obj[obj <= quant1] <- quant1
    obj[obj > quant2] <- quant2
    obj
  } else {
    apply(obj, 2, function(x) {
      quant1 <- quantile(x, 0 + trunc, na.rm = TRUE)
      quant2 <- quantile(x, 1 - trunc, na.rm = TRUE)
      x[x <= quant1] <- quant1
      x[x > quant2] <- quant2
      x
    })
  }
}


#' Higher order function to display progress.
#'
#' @param total Number of total runs in bootstrap
#'
#' @param nBars Number of bars to generate, default is 100.
#'
#' @param f Function used for progress bar.
#'
#' @param ... functions to be passed to function f.
#'
#' @references
#' https://stackoverflow.com/questions/27019190/txtprogressbar-for-parallel-bootstrap-not-displaying-properly
#'
#' @examples
#' \dontrun{
#' boot_res <- boot::boot(
#'   data = data, statistic =
#'     progressReporter(500, nBars = 100, f = boot_func), R = 500
#' )
#' }
#'
progressReporter <- function(total, nBars = 100, f, ...) {
  count <- 1
  step <- ceiling(total / nBars)
  cat(paste(rep("|", nBars), collapse = ""), "\r")
  flush.console()
  function(...) {
    if (count %% step == 0) {
      cat(".")
    }
    count <<- count + 1
    f(...)
  }
}
