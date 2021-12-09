# This function is meant to be called internally for optimization.
# It is called recursively (it calls itself) with each recursion decreasing
# the level of the requested likelihood by 1 (so, l decreases by one).
# @description
# calculates the log-likelihood of an l level mixed model using GH quadrature.
# the genral model is y = Xb + ZU + e
# but that values of beta and U are not included in the call. Instead this information
# is contained in yhat which incorporates Xb at the top level and all relevant 
# Zu information at lower levels.
# @param y numeric vector of y values 
# @param X numeric matrix of X values 
# @param Z a list of the Z matricies where the index of Z indicates the level.
#           of the Z matricies. Z[[1]] is NULL because there is no individual lavel Zs.
# @param ZFull Z expended such that each Z[[i]] matrix contains  one row for each observation (ie each element has same number of rows as original data)
# @param Qi  a list of the scaling factosr for the adaptive quadrature points (per group) where index indicated level 
# @param QiFull the scaling factors for adaptive quarature points duplicated so each element has one row per observation (like Zfull)
# @param omega numeric, list of the b estimates for each group
# @param omegaFull numeric, the b estimates for each group duplicated so each element has one row per observation (like Zfull)
# @param W list where each element is a numeric matrix of weights at that level, must have `w` and `index` columns
# @param k integer, the number of fixed effects 
# @param qp numeric vector (length= nQuad) the original quadrature points
# @param cConstructor function, the function used to construct the covariance matrices at each level 
# @param bobyqa boolean, indicating whether bobyqa should be used for optimization, default false
# @param verbose boolean, should additional output ususally used for debugging be printed, default true
# @param acc0 numeric, integer the accuracy for mfpr 
# @param mappedDefault boolean, when mapped is TRUE then variances are used as is, otherwise small (less than one) variances are exponentiated. This argument sets the default for a parameter of the returned function.
param.lnl.quad <- function(y, X, levels, Z, ZFull, Qi, QiFull, omega, omegaFull,
                           W, k, qp, cConstructor, bobyqa=FALSE, verbose=TRUE,
                           acc0, mappedDefault=FALSE, family=NULL) {
  function(par, acc=acc0, top=TRUE, integralMultiplierExponent=NULL,
           integralZColumn=1, mapped=mappedDefault, varFloor=NULL) {
    # par- the random and fixed effects and variance 
    # acc -  numeric, accuracy of the mpfr
    # top - whether the top level (ie overall likelihood rather than group likelihood) is to be returned. 
    # integralMultiplierExponent - a single integer, the function evaluates the
    #                                   integral times the random effect to this power
    #                                   when set to 0, it is just the log-likelhood
    #                                   when set to 1, this can be used to estimate the
    #                                   expected value.
    # integralZColumn is the column index of Z to use integralMultiplierExponent on.
    #                        only one random effect at a time can be integrated over, and the integration happens at the top level
    # mapped- parameter to be passed to covariance constructor which indicates wether very small variances are supposed to
    # mapped to exp(var-1) space 
    #
    
    parBeta <- par[1:k]
    if(is.null(varFloor)) {
      par <- par[ -(1:k)]
    } else {
      par <- pmax(varFloor, par[ -(1:k)])
    }
    parC <- cConstructor(par, mapped=mapped)
    yyh0 <- X %*% parBeta

    res <- calc.lin.lnl.quad(y=y, yhat=yyh0, level=levels, Z=Z, ZFull=ZFull,
                             Qi=Qi, QiFull=QiFull,
                             omega=omega, omegaFull=omegaFull, W=W,
                             C=parC, qp = qp, top=top,
                             verbose=verbose, acc=acc,
                             integralMultiplierExponent=integralMultiplierExponent,
                             integralZColumn=integralZColumn,
                             family=family)
    return(res)
  } # ends internal function call - this is what the top level function returns
}


# This function calcuates the liklihood of the model using integration by adaptive quadrature 
# @param y a numeric vector, the response.
# @param yhat the current predicted values of the response
# @param level an integer that respresents the number of levels in the likelihood
#          that is desired. In a two level model this function will be called
#          with l=2 and it will recurse and call itself with l=1
# @param Z a list of the Z matricies where the index of Z indicates the level.
#           of the Z matri. Z[[1]] is NULL because there is no individual lavel Zs.
# @param ZFull Z expended such that each Z[[i]] matrix contains  one row for each observation (ie each element has same number of rows as original data)
# @param Qi the scaling factor for the adaptive quadratures (per group)
# @param QiFull the scaling factors for adaptive quarature points duplicated so each element has one row per observation (like Zfull)
# @param omega a list of the b estimates for ech group
# @param omegaFull numeric, the b estimates for each group duplicated so each element has one row per observation (like Zfull)
# @param C a list of Cholesky decompositions of the Sigma matricies.
#          C[[1]] is simply the residual variance (a scalar) while C[[l]] for 
#          l > 1 is a matrix with the name number of rows and columns as the 
#          Z matrix for that level.
# @param qp Gaussian quadrature result from statmod::gauss.quad.
# @param W list of weight matricies. must have `w` and `index` columns
# @param top boolean set to TRUE to return a single scalar, otherwise returns a vector
# @param verbose boolean set to TRUE to get verbose output
# @param acc numeric, accuracy of the mpfr
# @param atPoint boolean, indicates likelihood should be calculated at single point
#                at the top level and then integrated below that. This is useful
#                for finding the maximum posterior (or likelihood) extimate for 
#                the random effects
# @param integralMultiplierExponent a single integer, the function evaluates the
#                                   integral times the randome effect to this power
#                                   when set to 0, it is just the log-likelhood
#                                   when set to 1, this can be used to estimate the
#                                   expected value.
# @param integralZColumn is the column index of Z to use integralMultiplierExponent on
#                        only one random effect at a time can be integrated over, and the integration happens at the top level
#
# @description
# calculates the log-likelihood of an l level mixed model using adaptive quadrature.
# the genral model is y = Xb + ZU + e
# but that values of beta and U are not included in the call. Instead this information
# is contained in yhat which incorporates Xb at the top level and all relevant 
# Zu information at lower levels.
# Applies the methodlogy of Rabe-Hesketh et al. 2005 (see pg 304-305)
#' @importFrom Rmpfr mpfr
#' @importFrom stats dnorm aggregate
#' @importFrom statmod gauss.quad
#' @importFrom NPflow mvnpdfC
calc.lin.lnl.quad <- function(y, yhat, level, Z, Qi, omega, W, C, qp,
                              top=TRUE, atPoint=FALSE, integralMultiplierExponent=NULL,
                              integralZColumn=1L,
                              verbose=TRUE, acc,
                              omegaFull=omega, QiFull=Qi, ZFull=Z, family) {
  ##################################
  ######### Outline ##############
  ### 1) Set up initial values 
  ### 2) Create grid of integration points
  ### 3) Evaluate likelihood at each point and apply weights 
  ### 4) sum likelihood and return results 

  ### 1) Set up initial values 
  #Bind y to the index present in the lowest level weights 
  data <- cbind(y, W[[1]])
  
  # grab matrixes for this level 
  # C stands for Covariance.
  #W stands for Weights and 
  # Z represents the random effects. 
  # Q represents the scaling factors 
  # the l indicates that these are lists with each element  having the covariance, weights or random effects for level l 
  # Notation follows Rabe-Hesketh, 2006 and Hartzel, 2001 (see main vignette for citations)
  Cl <- C[[level]] # covariance 
  Wl <- W[[level]] # weights
  Wlm1 <- W[[level-1]] # weights for next level down
  Zl <- Z[[level]] 
  ZFulll <- ZFull[[level]]
  Qil <- Qi[[level]]
  QiFulll <- QiFull[[level]]
  
  # replace aggregate() calls with matrix modification for speed. declaring arrays and matrices 
  ni_Wl <- table((Wlm1$indexp1)) 
  ntot <- length(Wlm1$indexp1)
  Wl_indices <- unique(Wlm1$indexp1)
  k_indices <- length(Wl_indices)
  # create block diagonal matrix for later multiplication
  diagM1 <- matrix(0, k_indices, sum(ni_Wl))
  runtot <- 0
  # fill block diagonal matrix 
  for (ii in 1:k_indices) { 
    diagM1[ii,] <- c(rep(0, runtot), rep(1, ni_Wl[ii]), rep (0, ntot - ni_Wl[ii] - runtot)) 
    runtot <- runtot + ni_Wl[ii]
  }
  
  ### 2) Create grid of integration points
  
  #Calcuate traditional quadrature points
  if(atPoint) {
    # if option atPoint is true whe use only one quadrature point 
    qp <- gauss.quad(1,"hermite")
    qp$weights <- 1
    nz <- nrow(Cl)
  }
  
  grd <- as.data.frame(qp) # create data frame of quadrature point(s) 
  colnames(grd)[1:2] <- paste0(c("v","w"),1)
  grd$w <- grd$w1
  
  # this creates a grid of integration points over which the likelihood will be evaluated
  if(nrow(Cl)>1) {
    for(i in 2:nrow(Cl)) {
      dfi <- grd
      colnames(dfi) <- paste0(c("v","w"),i)
      # by=NULL makes the Cartesian product 
      grd <- merge(grd, dfi, by=NULL)
      grd$w <- grd$w * grd[,paste0("w",i)]
    }
  }
  
 
  ### 3) Evaluate likelihood at each point and apply weights 
  Wl$l <- 0 #set inital value of likelihood to 0, this will accumulate the likelihood as we iterate through grid points

  # iterate over grid of integration points over which the likelihood will be evaluated
  for(i in 1:nrow(grd)) {
    # v is the IID N(0,I) vector at this integration point
    v <- t(grd[i, paste0("v", 1:ncol(Zl))])
   
    # this result is already weighted at the individual level
    if(level == 2) {
      # this calcuation follows the methodology described in Hartzel 2001, pg 87
      Wlm1$pts <- (omega[[level]] + sqrt(2) * matrix(t(v) %*% Qil, ncol=length(v), byrow=TRUE))
      # new predicted value of y 
      yyh <- yhat + rowSums(Zl * Wlm1$pts)
      # calculate the weigthts at this level based on the original wight and the  probabiilty 
      if(is.null(family)) {
        Wlm1$ll <- Wlm1$w * dnorm(y, mean=yyh, sd=C[[1]], log=TRUE)
      } else {
        Wlm1$ll <- family$lnl(y, family$linkinv(yyh), Wlm1$w, sd=C[[1]])
      }
    } else {
      # similar calcualation to level 2, but at he individual observaiton level to handle level 1 probabilty 
      # used below to get "prior" lnl
      Wlm1$pts <- (omega[[level]] + sqrt(2) * matrix(t(v) %*% Qil, ncol=length(v), byrow=TRUE))
      # used here to get lnl at these integration points
      pts <- (omegaFull[[level]] + sqrt(2) * matrix(t(v) %*% QiFulll, ncol=length(v), byrow=TRUE))
      yyh <- yhat + rowSums(ZFulll * pts)
      # set top=FALSE because it will not be the top level
      Wlm1$ll <- calc.lin.lnl.quad(y=y, yhat=yyh, level=level-1, Z=Z, Qi=Qi,
                                   omega=omega, W=W, C=C, qp=qp, top=FALSE,
                                   atPoint=FALSE, integralMultiplierExponent=NULL,
                                   verbose=verbose, acc=acc,
                                   omegaFull=omegaFull, QiFull=QiFull, ZFull=ZFull,
                                   family=family)
    } # ends the else part of the if (level==2 ) conditional 
    if(atPoint) {
      # replace aggregate() calls with matrix modification for speed. actual matrix multiplication
      agg2 <- diagM1 %*% Wlm1$ll
      agg <- cbind.data.frame(Wl_indices, agg2)
      #the above two lines replicate what was originally done with this line below
      #agg <- aggregate(ll ~ indexp1, Wlm1, sum) 
      colnames(agg) <- c("index", "ll")
      if(!top) {
        return(agg$ll)
      } else {
        return(sum(agg$ll))
      }
    } else {
      # evaluation of quadrature if not evaluating at a single point
      Wlm1NonD <- Wlm1[!duplicated(Wlm1$indexp1),] 
      # get normal probabilty of points under normal distribution  using multivariate normal pdf 
      Wlm1NonD$g_weight <- apply(Wlm1NonD$pts, MARGIN=1, function(p) {
        mvnpdfC(as.matrix(p), rep(0, length = length(p)), varcovM=Cl%*%t(Cl), Log=TRUE)
      })
      #For BLUE the integralMultiplierExponent is 1 whcih calcuates the expected value rather than the log likelihood 
      if(!is.null(integralMultiplierExponent)) {
        Wlm1NonD$x <- Wlm1NonD$pts[,integralZColumn]
        aggPrior <- Wlm1NonD[,c("indexp1", "g_weight", "x")]
      } else {
        aggPrior <- Wlm1NonD[,c("indexp1", "g_weight")]
      }
      
      # add likelihood for all observation in each group
      
      # replace aggregate() calls with matrix modification for speed. actual matrix multiplication
      agg2 <- diagM1 %*% Wlm1$ll
      agg <- cbind.data.frame(Wl_indices, agg2)
      #the above two lines replicate what was originally done with this line below
      #agg <- aggregate(ll ~ indexp1, Wlm1, sum) 
      colnames(agg) <- c("indexp1", "ll")
      names(agg)[!names(agg) %in% "indexp1"] <- "lli"
      
      agg <- merge(agg, aggPrior, by="indexp1")
      # exponetiate to reverse the log manipulations done earlier 
      # using log manipulation is important for maintaining precision when likelihood is very small 
      agg$li <- exp(mpfr(log(grd$w[i]) + agg$lli + agg$g_weight + sum(v*v), acc))
      if(!is.null(integralMultiplierExponent)) {
        agg$li <- agg$li * agg$x^integralMultiplierExponent
      }
      agg <- agg[,c("indexp1", "li")]
      Wl <- merge(Wl, agg, by.x="index", by.y="indexp1")
      Wl$l <- Wl$l + Wl$li * 2^(ncol(Zl)/2) * Wl$detQ
      Wl$li <- NULL
    }# ends the else part of the if(atPoint) conditional starting in line 192
  }
  # this turns the mpfr back to a double
  

  ### 4) sum likelihood and return results 
  if(!is.null(integralMultiplierExponent)) {
    res <- Wl$l
    return(res)
  }
  
  # log likelihood is weight * log(likelihood). As.numeric converts back from mfpr number to numeric
  Wl$ll <- Wl$w * as.numeric(log(Wl$l))
  # reset values that became - infintity 
  Wl$ll[Wl$ll== -Inf] <- -.Machine$double.xmax
  
  
  if(top) {
    # if user wants top level result (ie total likelihood sum all likelihood and return)
    return(sum(Wl$ll, na.rm=TRUE))
  }
  
  # other wise return likelihood of each group 
  res <- Wl$ll
  names(res) <- Wl$index
  return(res)
}
