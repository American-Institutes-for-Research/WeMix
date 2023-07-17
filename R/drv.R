# calculate a numerical central derivative 
d <- function(f, par, delta=1E-5) {
  res <- rep(NA, length(par))
  for(i in 1:length(par)) {
    d <- rep(0,length(par))
    d[i] <- delta/2
    fp <- f(par + d)
    fm <- f(par - d)
    res[i] <- (fp - fm) / delta
  }
  res
}


# efficently calculates the gradient and hessian
# per https://en.wikipedia.org/wiki/Finite_difference#Multivariate_finite_differences
getGradHess <- function(func, x, inputs=1:length(x), f0=NULL, diag=TRUE) {
  # use a wide net
  k <- length(inputs)
  if(is.null(f0)){
    f0 <- func(x)
  }
  # do not let h get smaller than the fourth root of the machine precision
  # but let it scale up when x is large
  h <- (abs(x) + 1) * (.Machine$double.eps)^0.25 
  # derivatives at a larger x (fp, for f plus h) and smaller x (fm for f minus h)
  hess <- matrix(NA, nrow=k, ncol=k)
  fxmi <- fxpi <- grad <- rep(NA, k)
  for(i in 1:k) {
    xm <- xp <- x
    xp[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
    # xm = (x[1], ..., x[i] - h[i], ...)
    xm[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
    fxpi[i] <- func(xp)
    fxmi[i] <- func(xm)
    grad[i] <- (fxpi[i] - fxmi[i])/(2*h[inputs[i]])
    hess[i,i] <- (fxpi[i] - 2 *f0 + fxmi[i])/(h[inputs[i]]*h[inputs[i]])
  }
  for(i in 1:(k-1)) {
    for(j in (i+1):k) {
      xpp <- xmm <- x
      # e.g. xpm is x with (x[1], x[2], ..., x[i] + h[i], ..., x[j] - h[j], ...)
      # so the first index (p, in the example) indicates x[i] is incrimented
      # by h[i] while the second index (m in the example) indicates x[j]
      # is decrimented by h[j].
      # xpp = (x[1], ..., x[i] + h[i], ..., x[j] + h[j], ...)
      xpp[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
      xpp[inputs[j]] <- x[inputs[j]] + h[inputs[j]]
      # xmm = (x[1], ..., x[i] - h[i], ..., x[j] - h[j], ...)
      xmm[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
      xmm[inputs[j]] <- x[inputs[j]] - h[inputs[j]]
      # find \partial^2 f(x) / \partial x[i] \partial x[j]
      res <- (func(xpp) - fxpi[i] - fxpi[j] + 2*f0 - fxmi[i] - fxmi[j] + func(xmm)) / (2* h[inputs[i]]*h[inputs[j]])
      hess[j,i] <- hess[i,j] <- res
    }
  }
  return(list(grad=grad,hess=hess))
}


# This function returns the second derivative of func with
# respect to the elements of x as a matrix
getHessian <- function(func, x, inputs=1:length(x), f0=NULL) {
  # use a wide net
  k <- length(inputs)
  if(is.null(f0)){
    f0 <- func(x)
  }
  if(length(f0) > 1) {
    return(lapply(1:length(f0),
                  function(i) {
                    f2 <- function(x) { 
                      func(x)[i]
                    }
                    getHessian(f2, x, inputs, f0=f0[i])
                  }))
  }
  # do not let h get smaller than the fourth root of the machine precision
  # but let it scale up when x is large
  h <- (abs(x) + 1) * (.Machine$double.eps)^0.25 
  # derivatives at a larger x (fp, for f plus h) and smaller x (fm for f minus h)
  hess <- matrix(NA, nrow=k, ncol=k)
  for(i in 1:k) {
    for(j in i:k) {
      xpp <- xmm <- xmp <- xpm <- xm <- xp <- x
      if(i!=j) {
        # e.g. xpm is x with (x[1], x[2], ..., x[i] + h[i], ..., x[j] - h[j], ...)
        # so the first index (p, in the example) indicates x[i] is incrimented
        # by h[i] while the second index (m in the example) indicates x[j]
        # is decrimented by h[j].
        # xpp = (x[1], ..., x[i] + h[i], ..., x[j] + h[j], ...)
        xpp[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
        xpp[inputs[j]] <- x[inputs[j]] + h[inputs[j]]
        # xmm = (x[1], ..., x[i] - h[i], ..., x[j] - h[j], ...)
        xmm[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
        xmm[inputs[j]] <- x[inputs[j]] - h[inputs[j]]
        # xmp = (x[1], ..., x[i] + h[i], ..., x[j] - h[j], ...)
        xpm[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
        xpm[inputs[j]] <- x[inputs[j]] - h[inputs[j]]
        # xmp = (x[1], ..., x[i] - h[i], ..., x[j] + h[j], ...)
        xmp[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
        xmp[inputs[j]] <- x[inputs[j]] + h[inputs[j]]
        # find \partial^2 f(x) / \partial x[i] \partial x[j]
        # Start with u(x) = \partial f(x) / \partial x[i] = ( f(xp_i) - f(xm_i)) / (2h_i)
        # and apply that same formula to hess[i,j] = \partial u(x) / \partial x[j]
        hess[j,i] <- hess[i,j] <- (func(xpp) - func(xpm) - func(xmp) + func(xmm)) / (4* h[inputs[i]]*h[inputs[j]])
      } else {
        # xp = (x[1], ..., x[i] + h[i], ...)
        xp[inputs[i]] <- x[inputs[i]] + h[inputs[i]]
        # xm = (x[1], ..., x[i] - h[i], ...)
        xm[inputs[i]] <- x[inputs[i]] - h[inputs[i]]
        # xpp = (x[1], ..., x[i] + 2*h[i], ...)
        xpp[inputs[i]] <- x[inputs[i]] + 2*h[inputs[i]]
        # xmm = (x[1], ..., x[i] - 2*h[i], ...)
        xmm[inputs[i]] <- x[inputs[i]] - 2*h[inputs[i]]
        # the above formula would suggest
        # (func(xpp) - 2*f0 + func(xmm)) / (4*h[inputs[i]]^2)
        # but this is more o(h^4) as opposed to o(h^2) above,
        # and the diagonal is critical, so use it
        hess[i,i] <- (-1*func(xmm) + 16*func(xp) - 30*f0 + 16*func(xm)-1*func(xpp)) / (12*h[inputs[i]]^2)
      }
    }
  }
  return(hess)
}
