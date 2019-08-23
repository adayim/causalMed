

## Code is from Hmisc package
rMultinom <- function(probs, m)
{
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if(!length(lev))
    lev <- 1:k

  ran <- matrix(lev[1], ncol=m, nrow=n)
  z <- apply(probs, 1, sum)
  if(any(abs(z-1) > .00001))
    stop('error in multinom: probabilities do not sum to 1')

  U <- apply(probs, 1, cumsum)
  for(i in 1:m)
  {
    un <- rep(runif(n), rep(k,n))
    ran[,i] <- lev[1 + apply(un > U, 2, sum)]
  }

  ran
}


# Extract data
extrCall <- function (x)
{
  if (isS4(x)) x@call else x$call
}


#' Extract bootstrap Confidence interval
extract_boot <- function(x,
                         conf.int = T,
                         conf.level = 0.95,
                         conf.method = "norm") {
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
    }
    else {
      op <- NULL
      for (i in index) op <- rbind(op, boot::imp.moments(boot.out,
                                                         index = i
      )$rat)
      op <- cbind(t0, op[, 1L] - t0, sqrt(op[, 2L]), apply(t,
                                                           2L, mean,
                                                           na.rm = TRUE
      ))
      colnames(op) <- c("statistic", "bias", "std.error", "estimate")
    }
  }
  # bring in rownames as "term" column, and turn into a data.frame
  ret <- data.frame(term = row.names(op),
                    as.data.frame(op),
                    row.names = NULL)

  if (conf.int) {
    ci.list <- lapply(seq_along(x$t0),
                      bootci,
                      boot.out = x,
                      conf = conf.level)
    ## stores them with longer names
    ci.tab <- do.call("rbind", ci.list)
    ci.tab <- data.frame(term = row.names(ci.tab),
                         as.data.frame(ci.tab),
                         row.names = NULL)

    ret <- merge(ret, ci.tab, by = "term", sort = FALSE)
  }
  ret
}

# Extact bootstrap confidence interval
bootci <- function (boot.out = NULL, conf = 0.95, index = 1, var.t0 = NULL,
                    t0 = NULL, t = NULL, L = NULL, h = function(t) t, hdot = function(t) 1,
                    hinv = function(t) t) {
  if (is.null(t0)) {
    if (!is.null(boot.out))
      t0 <- boot.out$t0[index]
    else stop("bootstrap output object or 't0' required")
  }
  if (!is.null(boot.out) && is.null(t))
    t <- boot.out$t[, index]
  if (!is.null(t)) {
    fins <- seq_along(t)[is.finite(t)]
    t <- h(t[fins])
  }
  if (is.null(var.t0)) {
    if (is.null(t)) {
      if (is.null(L))
        stop("unable to calculate 'var.t0'")
      else var.t0 <- sum((hdot(t0) * L/length(L))^2)
    }
    else var.t0 <- var(t)
  }
  else var.t0 <- var.t0 * hdot(t0)^2
  t0 <- h(t0)
  if (!is.null(t))
    bias <- mean(t) - t0
  else bias <- 0
  ci_len <- 2*qnorm((1+conf)/2)*sqrt(var.t0)
  merr <- sqrt(var.t0) * qnorm((1 + conf)/2)
  out <- cbind(conf.low = hinv(t0 - bias - merr),
               conf.high = hinv(t0 - bias + merr))
  out
}
