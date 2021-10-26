# constructs a list where each element is the variance covariance matrix at that level
# @param covMat data frame of variance and covariences  in the shape of the dataframe of the summary of the varcov term of results of lmer 
# @description
# Returns a function that takes the current parameter values and uses them to construct a list, 
# each element index i of the list corresponds to the variance covariance matrix at level i. 
# Variances less than .01 are transformed using the function exp(var-1) in order to prevent evaluation problems caused by variaces close to 0. 
# If internal parameter "mapped" is set to TRUE then this variance transformation is not performed.  
# @keywords internal
# @author Claire Kelley 
covMat2Cov <- function(covMat) {
  covMat$sdcor <- NULL
  # returnFull- if true return both the variance-covariance matrix and the data frame (in the shape of the variance df from lme) of variance and covariances
  #             used only in debugging
  # mapped - if mapped is TRUE then variances <.01 are NOT transformed using function exp(var-1)
  # this transformation is continuous and allows more accurate estimation of variances close to 0
  # so mapped defaults  FALSE 
  function(par, returnFull=FALSE, mapped=FALSE) {
    ############ Outline ############
    ## 1) handle mapping very small variances to exp(var-1) in order to get more accurate ests
    ## 2) reshape data from the pairwise data frame (ie rows are  Var1,Var2,Cov) provided by lme into 
    ##    list of square variance covariance matrixes for each level of the model 

    covMat$vcov <- par
    
    ## 1) handle mapping very small variances to exp(var-1) in order to get more accurate ests
    # any variances (ie has a blank in var 2) are exponentiatied if <- 0 to avoid problems around 0
    if(!mapped) {
      covMat$vcov <- ifelse(is.na(covMat$var2) & (covMat$vcov<1),exp(pmax(covMat$vcov-1,-4.6)),covMat$vcov)
    } else {
      if(any(is.na(covMat$var2) & (covMat$vcov < 0))) {
        stop("All variances must be positive.")
      }
      if(any(is.na(covMat$var2) & (covMat$vcov < exp(-4.6)))) {
        warning("Result may not be exact, variance too small. Try setting variance to a value over 0.01.")
      }
    }
    
    ## 2) reshape data
    # add the residual variance to the result list
    res <- list(s=sqrt(abs(covMat$vcov[covMat$grp=="Residual"]))) #level one variance is the residual variance 
    ulvl <- unique(covMat$level)
    ulvl <- ulvl[ulvl!=1] # remove level 1, we have already done that
    for(lvli in sort(ulvl)) {
      covMati <- covMat[covMat$level==lvli,]
      ki <- length(unique(covMati$var1))
      vci <- matrix(0, nrow=ki, ncol=ki)
      colnames(vci) <- rownames(vci) <- unique(covMati$var1)
      for(ii in 1:nrow(covMati)) {
        if(is.na(covMati$var2[ii])) {
          # when there is no var2 it is a diagonal element (ie variance rather than covariance)
          vci[rownames(vci)==covMati$var1[ii],colnames(vci)==covMati$var1[ii]] <- sqrt(abs(covMati$vcov[ii]))
        } else {
          # this is an off diagonal element (ie covariance), find its matrix spot by matching row and column names
          vci[rownames(vci)==covMati$var1[ii],colnames(vci)==covMati$var2[ii]] <- sign(covMati$vcov[ii]) * sqrt(abs(covMati$vcov[ii]))
          vci[rownames(vci)==covMati$var2[ii],colnames(vci)==covMati$var1[ii]] <- sign(covMati$vcov[ii]) * sqrt(abs(covMati$vcov[ii]))
        }
      }
      res <- c(res, list(vci))
    }
    if(returnFull) {
      return(list(cov=res, covMat=covMat))
    }
    return(res)
  }
}
