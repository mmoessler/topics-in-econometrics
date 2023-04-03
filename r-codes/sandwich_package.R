
library(sandwich)

# sandwich::vcovHAC(lm.res.02)
# # Error in umat - res : non-conformable arrays

sandwich___bwAndrews <- function (x, order.by = NULL,
                                  kernel = c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
                                  approx = c("AR(1)", "ARMA(1,1)"), weights = NULL, prewhite = 1, ar.method = "ols", data = list(), ...) {

  # # inputs
  # x <- lm.res.02
  # order.by <- NULL
  # kernel <- c("Quadratic Spectral")
  # approx <- c("AR(1)")
  # weights <- NULL
  # prewhite <- 1
  # ar.method <- "ols"
  # data = list()



  if (is.list(x) && !is.null(x$na.action)) {
    class(x$na.action) <- "omit"
  }
  kernel <- match.arg(kernel)
  approx <- match.arg(approx)
  prewhite <- as.integer(prewhite)
  umat <- if (inherits(x, "matrix")) {
    x
  } else {
    estfun(x)[, , drop = FALSE]
  }
  if (zoo::is.zoo(umat)) {
    umat <- as.matrix(coredata(umat))
  }
  n <- nrow(umat)
  k <- ncol(umat)
  if (!is.null(order.by)) {
    if (inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[, ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }
  umat <- umat[index, , drop = FALSE]

  # problem umat (n x m*k) - res ()
  if (is.null(weights)) {
    weights <- rep(1, k)
    unames <- colnames(umat)
    if (!is.null(unames) && "(Intercept)" %in% unames) {
      weights[which(unames == "(Intercept)")] <- 0
    } else {
      # res <- try(as.vector(rowMeans(estfun(x)/model.matrix(x), na.rm = TRUE)), silent = TRUE)
      # MM: added for m=2!
      tmp <- model.matrix(x)
      for (i in 1:ncol(residuals(x))) { tmp <- cbind(tmp, model.matrix(x)) }
      res <- try(as.vector(rowMeans(estfun(x)/cbind(model.matrix(x), model.matrix(x)), na.rm = TRUE)), silent = TRUE)
      if (inherits(res, "try-error")) {
        res <- try(residuals(x), silent = TRUE)
      }
      if (!inherits(res, "try-error")) {
        weights[which(colSums((umat - res)^2) < 1e-16)] <- 0
      }
    }
    if (isTRUE(all.equal(weights, rep(0, k)))) {
      weights <- rep(1, k)
    }
  } else {
    weights <- rep(weights, length.out = k)
  }
  if (length(weights) < 2) {
    weights <- 1
  }
  if (prewhite > 0) {
    var.fit <- ar(umat, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method)
    if (inherits(var.fit, "try-error")) {
      stop(sprintf("VAR(%i) prewhitening of estimating functions failed", prewhite))
    }
    umat <- as.matrix(na.omit(var.fit$resid))
    n <- n - prewhite
  }
  if (approx == "AR(1)") {
    fitAR1 <- function(x) {
      rval <- ar(x, order.max = 1, aic = FALSE, method = "ols")
      rval <- c(rval$ar, sqrt(rval$var.pred))
      names(rval) <- c("rho", "sigma")
      return(rval)
    }
    ar.coef <- apply(umat, 2, fitAR1)
    denum <- sum(weights * (ar.coef["sigma", ]/(1 - ar.coef["rho", ]))^4)
    alpha2 <- sum(weights * 4 * ar.coef["rho", ]^2 * ar.coef["sigma", ]^4/(1 - ar.coef["rho", ])^8)/denum
    alpha1 <- sum(weights * 4 * ar.coef["rho", ]^2 * ar.coef["sigma", ]^4/((1 - ar.coef["rho", ])^6 * (1 + ar.coef["rho", ])^2))/denum
  } else {
    fitARMA11 <- function(x) {
      rval <- arima(x, order = c(1, 0, 1), include.mean = FALSE)
      rval <- c(rval$coef, sqrt(rval$sigma2))
      names(rval) <- c("rho", "psi", "sigma")
      return(rval)
    }
    arma.coef <- apply(umat, 2, fitARMA11)
    denum <- sum(weights * ((1 + arma.coef["psi", ]) * arma.coef["sigma", ]/(1 - arma.coef["rho", ]))^4)
    alpha2 <- sum(weights * 4 * ((1 + arma.coef["rho", ] * arma.coef["psi", ]) * (arma.coef["rho", ] + arma.coef["psi", ]))^2 * arma.coef["sigma", ]^4/(1 - arma.coef["rho", ])^8)/denum
    alpha1 <- sum(weights * 4 * ((1 + arma.coef["rho", ] * arma.coef["psi", ]) * (arma.coef["rho", ] + arma.coef["psi", ]))^2 * arma.coef["sigma", ]^4/((1 - arma.coef["rho", ])^6 * (1 + arma.coef["rho", ])^2))/denum
  }

  rval <- switch(kernel,
                 Truncated = { 0.6611 * (n * alpha2)^(1/5) },
                 Bartlett = { 1.1447 * (n * alpha1)^(1/3) },
                 Parzen = { 2.6614 * (n * alpha2)^(1/5) },
                 `Tukey-Hanning` = { 1.7462 * (n * alpha2)^(1/5) },
                 `Quadratic Spectral` = { 1.3221 * (n * alpha2)^(1/5) } )

  return(rval)

}

sandwich__kweights <- function (x,
                                kernel = c("Truncated", "Bartlett", "Parzen", "Tukey-Hanning", "Quadratic Spectral"),
                                normalize = FALSE) {

  kernel <- match.arg(kernel, c("Truncated", "Bartlett", "Parzen", "Tukey-Hanning", "Quadratic Spectral"))
  if (normalize) {
    ca <- switch(kernel,
                 Truncated = 2,
                 Bartlett = 2/3,
                 Parzen = 0.539285,
                 `Tukey-Hanning` = 3/4,
                 `Quadratic Spectral` = 1)
  } else {
    ca <- 1
  }
  switch(kernel,
         Truncated = { ifelse(ca * abs(x) > 1, 0, 1) },
         Bartlett = { ifelse(ca * abs(x) > 1, 0, 1 - abs(ca * x)) },
         Parzen = { ifelse(ca * abs(x) > 1, 0, ifelse(ca * abs(x) < 0.5, 1 - 6 * (ca * x)^2 + 6 * abs(ca * x)^3, 2 * (1 - abs(ca * x))^3)) },
         `Tukey-Hanning` = { ifelse(ca * abs(x) > 1, 0, (1 + cos(pi * ca * x))/2) },
         `Quadratic Spectral` = {
           qs <- function(x) {
             y <- 6 * pi * x/5
             3 * (1/y)^2 * (sin(y)/y - cos(y))
           }
           w <- qs(x)
           if (length(ix <- which(abs(x) < 0.001)) > 0L) {
             cf <- 1e+06 * log(qs(0.001))
             w[ix] <- exp(cf * x[ix]^2)
           }
           w}
  )

}

sandwich__weightsAndrews <- function (x, order.by = NULL, bw = bwAndrews,
                                      kernel = c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
                                      prewhite = 1, ar.method = "ols", tol = 1e-07, data = list(), verbose = FALSE, ...) {

  # # inputs
  # x <- lm.res.02
  # order.by <- NULL
  # bw <- bwAndrews
  # kernel <- NULL
  # prewhite <- 1
  # ar.method <- "ols"
  # tol <- 1e-07
  # verbose <- FALSE
  # data = list()



  if (is.list(x) && !is.null(x$na.action)) {
    class(x$na.action) <- "omit"
  }
  kernel <- match.arg(kernel)

  # call bandwidth function bw(), bwAndrews()
  if (is.function(bw)) {
    # bw <- bw(x, order.by = order.by, kernel = kernel, prewhite = prewhite, data = data, ar.method = ar.method, ...)
    # bw <- bw(x, order.by = order.by, kernel = kernel, prewhite = prewhite, data = data, ar.method = ar.method)
    bw <- sandwich___bwAndrews(x, order.by = order.by, kernel = kernel, prewhite = prewhite, data = data, ar.method = ar.method)
  }

  if (verbose) {
    cat(paste("\nBandwidth chosen:", format(bw), "\n"))
  }

  n <- NROW(estfun(x)) - as.integer(prewhite)
  weights <- kweights(0:(n - 1)/bw, kernel = kernel)
  weights <- weights[1:max(which(abs(weights) > tol))]

  return(weights)

}

sandwich__meatHAC <- function(x, order.by = NULL, prewhite = FALSE, weights = weightsAndrews,
                              adjust = TRUE, diagnostics = FALSE, ar.method = "ols",
                              data = list(), ...) {

  # # inputs
  # x <- lm.res.02
  # order.by <- NULL
  # prewhite <- FALSE
  # weights <- weightsAndrews
  # adjust <- TRUE
  # diagnostics <- FALSE
  # ar.method <- "ols"
  # data <- list()



  if (is.list(x) && !is.null(x$na.action)) {
    class(x$na.action) <- "omit"
  }
  prewhite <- as.integer(prewhite)
  # umat <- estfun(x, ...)[, , drop = FALSE]
  umat <- estfun(x)[, , drop = FALSE]
  if (zoo::is.zoo(umat)) {
    umat <- as.matrix(coredata(umat))
  }
  # # MM: check
  # head(xe.01)
  # head(umat) # (t x m*k)

  n.orig <- n <- nrow(umat)
  k <- ncol(umat)
  if (!is.null(order.by)) {
    if (inherits(order.by, "formula")) {
      z <- model.matrix(order.by, data = data)
      z <- as.vector(z[, ncol(z)])
    } else {
      z <- order.by
    }
    index <- order(z)
  } else {
    index <- 1:n
  }
  umat <- umat[index, , drop = FALSE]
  if (prewhite > 0) {
    var.fit <- try(ar(umat, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method))
    if (inherits(var.fit, "try-error")) {
      stop(sprintf("VAR(%i) prewhitening of estimating functions failed", prewhite))
    }
    if (k > 1) {
      D <- solve(diag(ncol(umat)) - apply(var.fit$ar, 2:3, sum))
    } else {
      D <- as.matrix(1/(1 - sum(var.fit$ar)))
    }
    umat <- as.matrix(na.omit(var.fit$resid))
    n <- n - prewhite
  }

  # call weight function weights(), weightsAndrews()
  if (is.function(weights)) {
    # weights <- weights(x, order.by = order.by, prewhite = prewhite, ar.method = ar.method, data = data)
    weights <- sandwich__weightsAndrews(x, order.by = order.by, prewhite = prewhite, ar.method = ar.method, data = data)
  }
  if (length(weights) > n) {
    warning("more weights than observations, only first n used")
    weights <- weights[1:n]
  }
  utu <- 0.5 * crossprod(umat) * weights[1]
  wsum <- n * weights[1]/2
  w2sum <- n * weights[1]^2/2
  if (length(weights) > 1) {
    for (ii in 2:length(weights)) {
      utu <- utu + weights[ii] * crossprod(umat[1:(n - ii + 1), , drop = FALSE], umat[ii:n, , drop = FALSE])
      wsum <- wsum + (n - ii + 1) * weights[ii]
      w2sum <- w2sum + (n - ii + 1) * weights[ii]^2
    }
  }
  utu <- utu + t(utu)
  if (adjust) {
    utu <- n.orig/(n.orig - k) * utu
  }
  if (prewhite > 0) {
    utu <- crossprod(t(D), utu) %*% t(D)
  }
  wsum <- 2 * wsum
  w2sum <- 2 * w2sum
  bc <- n^2/(n^2 - wsum)
  df <- n^2/w2sum
  rval <- utu/n.orig
  if (diagnostics) {
    attr(rval, "diagnostics") <- list(bias.correction = bc, df = df)
  }

  return(rval)

}
