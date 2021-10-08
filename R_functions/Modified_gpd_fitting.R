
library(maxLik)

## Modified Peaks Over Threshold GPD fitting, for box and linear constraints

gpd_constraints_fit_optim <- function (xdat, threshold, npy = 365, ydat = NULL, sigl = NULL, 
                                       shl = NULL, siglink = identity, shlink = identity, siginit = NULL, 
                                       shinit = NULL, show = TRUE, method = "L-BFGS-B", maxit = 100000, lower = -Inf, upper = Inf, 
                                       ...){
  ## Modified from original implementation in the ismev R package: https://cran.r-project.org/web/packages/ismev/
  ## to allow box constraints during fitting
  z <- list()
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  n <- length(xdat)
  z$trans <- FALSE
  if (is.function(threshold)) 
    stop("`threshold' cannot be a function")
  u <- rep(threshold, length.out = n)
  if (length(unique(u)) > 1) 
    z$trans <- TRUE
  xdatu <- xdat[xdat > u]
  xind <- (1:n)[xdat > u]
  u <- u[xind]
  in2 <- sqrt(6 * var(xdatu))/pi
  in1 <- mean(xdatu, na.rm = TRUE) - 0.57722 * in2
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(siginit)) 
      siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
    if (is.null(siginit)) 
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(shinit)) 
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
    if (is.null(shinit)) 
      shinit <- c(0.1, rep(0, length(shl)))
  }
  init <- c(siginit, shinit)
  z$model <- list(sigl, shl)
  z$link <- deparse(substitute(c(siglink, shlink)))
  z$threshold <- threshold
  z$nexc <- length(xdatu)
  z$data <- xdatu
  gpd.lik <- function(a) {
    sc <- siglink(sigmat %*% (a[seq(1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npsc + 1, length = npsh)]))
    y <- (xdatu - u)/sc
    y <- 1 + xi * y
    if (min(sc) <= 0) 
      l <- 10^6
    else {
      if (min(y) <= 0) 
        l <- 10^6
      else {
        l <- sum(log(sc)) + sum(log(y) * (1/xi + 1))
      }
    }
    l
  }
  x <- optim(init, gpd.lik, hessian = TRUE, method = method, lower = lower, upper = upper,
             control = list(maxit = maxit, ...))
  sc <- siglink(sigmat %*% (x$par[seq(1, length = npsc)]))
  xi <- shlink(shmat %*% (x$par[seq(npsc + 1, length = npsh)]))
  z$conv <- x$convergence
  z$nllh <- x$value
  z$vals <- cbind(sc, xi, u)
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdatu - u))/sc)^(-1/xi)))
  }
  z$mle <- x$par
  z$rate <- length(xdatu)/n
  z$cov <- solve(x$hessian)
  z$se <- sqrt(diag(z$cov))
  z$n <- n
  z$npy <- npy
  z$xdata <- xdat
  if (show) {
    if (z$trans) 
      print(z[c(2, 3)])
    if (length(z[[4]]) == 1) 
      print(z[4])
    print(z[c(5, 7)])
    if (!z$conv) 
      print(z[c(8, 10, 11, 13)])
  }
  z$OPTIMDEBUG <- x
  class(z) <- "gpd.fit"
  invisible(z)
}


gpd_constraints_fit_maxLik <- function (xdat, threshold, npy = 365, ydat = NULL, sigl = NULL, 
                                        shl = NULL, siglink = identity, shlink = identity, siginit = NULL, 
                                        shinit = NULL, show = TRUE, method=NULL, maxit = 100000,
                                        grad=NULL, constraints=NULL, ...){
  ## Modified from original implementation in the ismev R package: https://cran.r-project.org/web/packages/ismev/
  ## to allow linear constraints during fitting
  z <- list()
  npsc <- length(sigl) + 1
  npsh <- length(shl) + 1
  n <- length(xdat)
  z$trans <- FALSE
  if (is.function(threshold)) 
    stop("`threshold' cannot be a function")
  u <- rep(threshold, length.out = n)
  if (length(unique(u)) > 1) 
    z$trans <- TRUE
  xdatu <- xdat[xdat > u]
  xind <- (1:n)[xdat > u]
  u <- u[xind]
  in2 <- sqrt(6 * var(xdatu))/pi
  in1 <- mean(xdatu, na.rm = TRUE) - 0.57722 * in2
  if (is.null(sigl)) {
    sigmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(siginit)) 
      siginit <- in2
  }
  else {
    z$trans <- TRUE
    sigmat <- cbind(rep(1, length(xdatu)), ydat[xind, sigl])
    if (is.null(siginit)) 
      siginit <- c(in2, rep(0, length(sigl)))
  }
  if (is.null(shl)) {
    shmat <- as.matrix(rep(1, length(xdatu)))
    if (is.null(shinit)) 
      shinit <- 0.1
  }
  else {
    z$trans <- TRUE
    shmat <- cbind(rep(1, length(xdatu)), ydat[xind, shl])
    if (is.null(shinit)) 
      shinit <- c(0.1, rep(0, length(shl)))
  }
  init <- c(siginit, shinit)
  z$model <- list(sigl, shl)
  z$link <- deparse(substitute(c(siglink, shlink)))
  z$threshold <- threshold
  z$nexc <- length(xdatu)
  z$data <- xdatu
  gpd.lik <- function(a) {
    sc <- siglink(sigmat %*% (a[seq(1, length = npsc)]))
    xi <- shlink(shmat %*% (a[seq(npsc + 1, length = npsh)]))
    y <- (xdatu - u)/sc
    y <- 1 + xi * y
    if (min(sc) <= 0) 
      l <- -10^6
    else {
      if (min(y) <= 0) 
        l <- -10^6
      else {
        l <- -sum(log(sc)) - sum(log(y) * (1/xi + 1))
      }
    }
    l
  }
  if (is.null(method)) {
    if (is.null(constraints)) {
      method <- "nr"
    }
    else if (identical(names(constraints), c("ineqA", "ineqB"))) {
      if (is.null(grad)) 
        method <- "Nelder-Mead"
      else method <- "BFGS"
    }
    else method <- "nr"
  }
  x <- maxLik(gpd.lik, grad=grad, start=init, method=method, constraints=constraints, control=list(iterlim = maxit, ...))
  z$OPTIMDEBUG <- x
  #return(z)
  sc <- siglink(sigmat %*% (coef(x)[seq(1, length = npsc)]))
  xi <- shlink(shmat %*% (coef(x)[seq(npsc + 1, length = npsh)]))
  z$conv <- returnCode(x)
  z$nllh <- maxValue(x)
  z$vals <- cbind(sc, xi, u)
  if (z$trans) {
    z$data <- -log(as.vector((1 + (xi * (xdatu - u))/sc)^(-1/xi)))
  }
  z$mle <- coef(x)
  z$rate <- length(xdatu)/n
  z$cov <- tryCatch(solve(hessian(x)), error=function(cond) {
    #warning("Singular covariance matrix in GPD fit: ", cond)#not corrected in the constrained case anyways (avoid spam)
    return(NULL)
  }, warning=function(cond) {
    warning("Warning with covariance matrix in GPD fit: ", cond)
    return(suppressWarnings(solve(hessian(x))))
  }, finally = {})
  if(!is.null(z$cov)){
    z$se <- sqrt(diag(z$cov))
  }else{
    z$se <- NULL
  }
  z$n <- n
  z$npy <- npy
  z$xdata <- xdat
  if (show) {
    if (z$trans) 
      print(z[c(2, 3)])
    if (length(z[[4]]) == 1) 
      print(z[4])
    print(z[c(5, 7)])
    if (!z$conv) 
      print(z[c(8, 10, 11, 13)])
  }
  z$returnMessage <- returnMessage(x)
  class(z) <- "gpd.fit"
  invisible(z)
}
#optim: x $ par convergence value hessian
#makLik: coef returnCode,returnMessage maxValue hessian (x)
