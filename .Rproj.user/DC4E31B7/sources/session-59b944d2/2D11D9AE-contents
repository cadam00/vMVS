vm.reg <- function(y, x, rads = TRUE, xnew = NULL, tol = 1e-6, maxiters = 500) {

  if ( !is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    y <- cbind( cos(y), sin(y) )
  }

  x <- model.matrix( ~., data.frame(x) )
  dm <- dim(x)
  n <- dm[1]  ;  p <- dm[2] - 1

  be <- solve( crossprod(x), crossprod(x, y) )  ## initial values for the beta
  mu <- x %*% be
  ki <- sqrt( Rfast::rowsums(mu^2) )
  lik1 <- sum(mu * y) - sum( log(besselI(ki, 0, expon.scaled = TRUE) ) + ki )

  A1  <- besselI(ki, 1, expon.scaled = TRUE) / besselI(ki, 0, expon.scaled = TRUE)
  A1d <- 1 - A1 / ki - A1^2
  ## --- gradient: p x 2 ---
  resid <- y - (A1 / ki) * mu
  grad <- crossprod(x, resid)
  ## --- Hessian weights ---
  mu_ki  <- mu / ki
  d11 <- A1d * mu_ki[, 1]^2 + (A1 / ki) * (1 - mu_ki[, 1]^2)
  d22 <- A1d * mu_ki[, 2]^2 + (A1 / ki) * (1 - mu_ki[, 2]^2)
  d12 <- (A1d - A1/ki) * mu_ki[, 1] * mu_ki[, 2]
  ## --- Hessian blocks: each is p x p ---
  H11 <-  -crossprod(x * d11, x)
  H22 <-  -crossprod(x * d22, x)
  H12 <-  -crossprod(x * d12, x)
  ## --- full 2p x 2p Hessian ---
  H <- rbind( cbind(H11, H12), cbind(H12, H22) )

  be <- be - solve(H, as.vector(grad) )
  be <- matrix(be, ncol = 2)
  mu <- x %*% be
  ki <- sqrt( Rfast::rowsums(mu^2) )
  lik2 <- sum(mu * y) - sum( log(besselI(ki, 0, expon.scaled = TRUE) ) + ki )

  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters ) {
    i <- i + 1
    lik1 <- lik2

    A1  <- besselI(ki, 1, expon.scaled = TRUE) / besselI(ki, 0, expon.scaled = TRUE)
    A1d <- 1 - A1 / ki - A1^2
    ## --- gradient: p x 2 ---
    resid <- y - (A1 / ki) * mu
    grad <- crossprod(x, resid)
    ## --- Hessian weights ---
    mu_ki  <- mu / ki
    d11 <- A1d * mu_ki[, 1]^2 + (A1 / ki) * (1 - mu_ki[, 1]^2)
    d22 <- A1d * mu_ki[, 2]^2 + (A1 / ki) * (1 - mu_ki[, 2]^2)
    d12 <- (A1d - A1/ki) * mu_ki[, 1] * mu_ki[, 2]
    ## --- Hessian blocks: each is p x p ---
    H11 <-  -crossprod(x * d11, x)
    H22 <-  -crossprod(x * d22, x)
    H12 <-  -crossprod(x * d12, x)
    ## --- full 2p x 2p Hessian ---
    H <- rbind( cbind(H11, H12), cbind(H12, H22) )

    be <- be - solve(H, as.vector(grad) )
    be <- matrix(be, ncol = 2)
    mu <- x %*% be
    ki <- sqrt( Rfast::rowsums(mu^2) )
    lik2 <- sum(mu * y) - sum( log(besselI(ki, 0, expon.scaled = TRUE) ) + ki )
  }

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew))
    est <- xnew %*% be
    est <- est / sqrt( Rfast::rowsums(est^2) )
  }
  loglik <- lik2 - n * log(2 * pi)
  list(be = be, loglik = loglik, ki = ki, A1 = A1, est = est)
}






