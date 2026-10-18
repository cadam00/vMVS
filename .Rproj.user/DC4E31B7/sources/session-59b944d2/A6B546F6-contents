vm.gomp <- function(y, x, rads = TRUE, xstand = TRUE, alpha = 0.05, R = 499) {

  runtime <- proc.time()
  if ( !is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    y <- cbind( cos(y), sin(y) )
  }
  if ( xstand )  x <- Rfast::standardise(x)
  dm <- dim(x)   ;   n <- dm[1]   ;   p <- dm[2]

  ind <- 1:p
  lik <- Rfast::vmf.mle(y, tol = 1e-6)$loglik
  ela <- Rfast::eachcol.apply(x, y[, 1])^2 + Rfast::eachcol.apply(x, y[, 2])^2
  sel <- which.max(ela)
  sela <- NULL
  rho <- c( 0, .perm.reg0(y, x[, sel], R = R) )

  i <- 2
  while ( rho[i] < alpha  &  length(sela) < p - 1  &  i < n ) {
    i <- i + 1
    sela <- c(sela, sel)
    ind[sela] <- 0
    mod <- vm.reg(y, x[, sela], xnew = x[, sela])
    lik[i - 1] <- mod$loglik
    res <- y - mod$est * mod$A1
    r <- numeric(p)
    ii <- ind[ind > 0]
    r[ind] <- Rfast::eachcol.apply(x, res[, 1], indices = ii)^2 +
              Rfast::eachcol.apply(x, res[, 2], indices = ii)^2
    sel <- which.max(r)
    rho[i] <- .perm.reg1(y, x = x[, sel], xsel = x[, sela], R = R)
  }

  mod <- vm.reg(y, x[, sela])
  runtime <- proc.time() - runtime
  len <- length(sela)
  res <- cbind( c(0, sela), lik[1:(len + 1)], c(NA, rho[seq_len(len) + 1]) )
  colnames(res) <- c("Selected", "log-likelihood", "p-value")
  list(runtime = runtime, res = res, mod = mod)
}




.perm.reg0 <- function(y, x, tol = 1e-6, maxiters = 500, R = 999) {

  n <- dim(y)[1]
  lik0 <- Rfast::vmf.mle(y, tol = tol)$loglik
  stat <- vm.reg(y = y, x = x, tol = tol, maxiters = maxiters)$loglik - lik0

  bstat <- numeric(R)
  for ( i in 1:R ) {
    yb <- y[sample.int(n, n), ]
    bstat[i] <- vm.reg(y = yb, x = x, tol = tol, maxiters = maxiters)$loglik - lik0
  }
  ( sum(bstat >= stat) + 1 ) / (R + 1)
}


.perm.reg1 <- function(y, x, xsel, tol = 1e-6, maxiters = 500, R = 999) {
  n <- dim(y)[1]
  xf <- cbind(xsel, x)

  mod0 <- vm.reg(y = y, x = xsel, xnew = xsel, tol = tol, maxiters = maxiters)
  u <- mod0$est
  stat <- vm.reg(y = y, x = xf, tol = tol, maxiters = maxiters)$loglik - mod0$loglik

  e <- cbind( y[, 1] * u[, 1] + y[, 2] * u[, 2], y[, 2] * u[, 1] - y[, 1] * u[, 2] )

  bstat <- numeric(R)
  for ( i in 1:R ) {
    eb <- e[sample.int(n, n), ]
    yb <- cbind( eb[, 1] * u[, 1] - eb[, 2] * u[, 2], eb[, 1] * u[, 2] + eb[, 2] * u[, 1] )
    bstat[i] <- vm.reg(y = yb, x = xf, tol = tol, maxiters = maxiters)$loglik -
      vm.reg(y = yb, x = xsel, tol = tol, maxiters = maxiters)$loglik
  }
  ( sum(bstat >= stat) + 1 ) / (R + 1)
}





