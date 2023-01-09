#' Survey Weighted Mixed-Effects Models 
#' @param formula  a formula object in the style of \code{lme4} that creates the model.
#' @param data  a data frame containing the raw data for the model.
#' @param weights a character vector of names of weight variables found in the data frame starts with units (level 1) and increasing (larger groups).
#' @param cWeights logical, set to \code{TRUE} to use conditional weights. Otherwise, \code{mix} expects unconditional weights.
#' @param center_group a list where the name of each element is the name of the aggregation level, and the element is a formula of
#'  variable names to be group mean centered; for example to group mean center gender and age within the group student:
#'   \code{list("student"= ~gender+age)}, default value of NULL does not perform any group mean centering. 
#' @param center_grand  a formula of variable names  to be grand mean centered, for example to center the variable education by overall mean of 
#' education: \code{~education}. Default is NULL which does no centering. 
#' @param nQuad an optional integer  number of quadrature points to evaluate models solved by adaptive quadrature.
#'  Only non-linear models are evaluated with adaptive quadrature. See notes for additional guidelines. 
#' @param max_iteration a optional integer, for non-linear models fit by adaptive quadrature which limits number of iterations allowed
#'   before quitting. Defaults  to 10. This is used because if the liklihood surface is flat, 
#'   models may run for a very  long time without converging.
#' @param run  logical; \code{TRUE} runs the model while \code{FALSE} provides partial output for debugging or testing. Only applies to non-linear
#'  models evaluated by adaptive quadrature. 
#' @param keepAdapting logical, set to \code{TRUE} when the adaptive quadrature should adapt after every Newton step. Defaults to \code{FALSE}. 
#' \code{FALSE} should be used for faster (but less accurate) results. Only applies to non-linear models. 
#' @param verbose logical, default \code{FALSE}; set to \code{TRUE} to print results of intermediate steps of adaptive quadrature. Only applies to non-linear models.
#' @param start optional numeric vector representing the point at which the model should start optimization; takes the shape of c(coef, vars) 
#' from results (see help). 
#' @param family the family; optionally used to specify generalized linear mixed models. Currently only \code{binomial(link="logit")} is supported.
#' @param acc0 deprecated; ignored. 
#' @param fast logical; deprecated
#' @description
#' Implements a survey weighted mixed-effects model using the provided formula. 
#' @details
#' Linear models are solved using a modification of the analytic solution developed by Bates and Pinheiro (1998).
#' Non-linear models are solved using adaptive quadrature following the method in STATA's GLAMMM (Rabe-Hesketh & Skrondal, 2006).
#' For additional details, see the vignettes \emph{Weighted Mixed Models: Adaptive Quadrature} and  \emph{Weighted Mixed Models: Analytical Solution} 
#' which provide extensive examples as well as a description of the mathematical basis of the estimation procedure and comparisons to model 
#' specifications in other common software. 
#' Notes: 
#' \itemize{
#' \item Standard errors of random effect variances are robust; see vignette for details. 
#' \item To see the function that is maximized in the estimation of this model, see the section on "Model Fitting" in the
#'  \emph{Introduction to Mixed Effect Models With WeMix} vignette.
#' \item When all weights above the individual level are 1, this is similar to a \code{lmer} and you should use \code{lme4} 
#' because it is much faster.
#' \item If  starting coefficients are not provided they are estimated using \code{lme4}. 
#' \item For non-linear models, when the variance of a random effect is very low (<.1), WeMix doesn't estimate it, because very 
#' low variances create problems with  numerical evaluation. In these cases, consider estimating without that random effect. 
#'  \item The model is estimated by maximum likelihood estimation.
#' \item To choose the number of quadrature points for non-linear model evaluation, a balance is needed between accuracy and
#' speed; estimation time increases quadratically with the number of points chosen. In addition, an odd number of points is 
#' traditionally used. We recommend starting at 13 and increasing or decreasing as needed. 
#' }
#' @importFrom lme4 getME lmer glmer lFormula GHrule
#' @importFrom stats dnorm aggregate terms dpois dgamma dbinom ave model.matrix terms.formula as.formula sigma complete.cases update
#' @importFrom numDeriv genD hessian grad
#' @importFrom minqa bobyqa 
#' @importFrom Matrix nearPD sparse.model.matrix Cholesky Matrix .updateCHMfactor
#' @importFrom matrixStats rowLogSumExps
#' @return object of class \code{WeMixResults}. 
#' This is a list with elements: 
#' \item{lnlf}{function, the likelihood function.} 
#' \item{lnl}{numeric, the log-likelihood of the model. }
#' \item{coef}{numeric vector, the estimated coefficients of the model. }
#' \item{ranefs}{the group-level random effects.}
#' \item{SE}{the standard errors of the fixed effects, robust if robustSE was set to true.}
#' \item{vars}{numeric vector, the random effect variances.}
#' \item{theta}{the theta vector.}
#' \item{call}{the original call used.}
#' \item{levels}{integer, the number of levels in the model.}
#' \item{ICC}{numeric, the intraclass correlation coefficient.}
#' \item{CMODE}{the conditional mean of the random effects.}
#' \item{invHessian}{inverse of the second derivative of the likelihood function.}
#' \item{ICC}{the interclass correlation.}
#' \item{is_adaptive}{logical, indicates if adaptive quadrature was used for estimation.}
#' \item{sigma}{the sigma value.}
#' \item{ngroups}{the number of observations in each group.}
#' \item{varDF}{the variance data frame in the format of the variance data frame returned by lme4.}
#' \item{varVC}{the variance-covariance matrix of the random effects.}
#' \item{cov_mat}{the variance-covariance matrix of the fixed effects.}
#' \item{var_theta}{the variance covariance matrix of the theta terms.}
#' \item{wgtStats}{statistics regarding weights, by level.}
#' \item{ranefMat}{list of matrixes; each list element is a matrix of random effects by level with IDs in the rows and random effects in the columns.}
#' @example \man\examples\mix.R
#' @author Paul Bailey, Blue Webb, Claire Kelley, and Trang Nguyen 
#' @export
mix <- function(formula, data, weights, cWeights=FALSE, center_group=NULL,
                center_grand=NULL, max_iteration=10, nQuad=13L, run=TRUE,
                verbose=FALSE, acc0=120, keepAdapting=FALSE, start=NULL,
                fast=FALSE, family=NULL) {
  #############################################
  #                   Outline:                #   
  #     1) data preparation and reshaping     #
  #    2) Identify integration parameters     #
  #     3) Maximum Likelihood Estimation      #
  #            4) Post Estimation             #
  #############################################
  
  #############################################
  #              1) Data Preparation          #
  #############################################
  call <- match.call()
  # check class and some requirements of arguments
  if(!inherits(formula, "formula")) stop(paste0("The argument ", sQuote("formula"), " must be a formula."))
  if(!inherits(data, "data.frame")) stop(paste0("The argument ", sQuote("data"), " must be a data.frame."))  
  if(length(class(data)) > 1) { 
    data <- as.data.frame(data)
  }
  if(!missing(acc0)) message("The argument acc0 is now deprecated and will be ignored.")
  if(nQuad <= 0) stop(paste0("The argument ", sQuote("nQuad"), " must be a positive integer."))
  if(!inherits(run, "logical")) stop(paste0("The argument ", sQuote("run"), " must be a logical."))
  if(!inherits(verbose, "logical")) stop(paste0("The argument ", sQuote("verbose"), " must be a logical."))
  if(!inherits(weights, "character")) stop(paste0("The argument ", sQuote("weights"), " must be a character vector of weight column names in ", sQuote("data"), "."))
  if(any(!weights %in% colnames(data))) stop(paste0("The argument ", sQuote("weights"), " must specify valid columns in ", sQuote("data"), ". Could not find column(s): ", paste(dQuote(weights[!weights %in% colnames(data)]), collapse=", "), "."))
  if(!missing(fast)) warning(paste0("The ", sQuote("fast"), " argument is deprecated."))
  if(any(grepl("[|].*[.]",attributes(terms(formula))$term.labels))) stop("The formula is not valid for mix. The name of conditioning variables must not contain a dot. Try renaming variables after | in the fomrula so they do not contain a dot.")
  #currently wemix can only use complete cases
  #this removes any  incomplete cases if they exist 
  if(any(is.na(data[ , c(all.vars(formula), weights)]))) {
    cl <- call("model.frame",
               formula=formula(paste0("~", paste0(unique(c(all.vars(formula), weights)),collapse=" + "))),
               data=data)
    dt <- eval(cl, parent.frame(1L))
    warning(paste0("There were ", sum(nrow(data) - nrow(dt)), " rows with missing data. These have been removed."))
    data <- dt 
    rm(dt)
  }
  weights0 <- weights # keep track of incomming weights column names
  data <- dropNonPositiveWeights(data, weights)

  family <- setup_family(family)

  # set up initial values 
  adapter <- "MAP" # initial function evaluation through MAP, BLUE estimator can be used post hoc
  nQuad <- round(nQuad) #nquad must be an integer
  
  # get the group names (ie level 2+ variables) from the formula
  # start by getting the lme4 parsed formula
  lformula <- lFormula(formula=formula, data=data)
  # get the unparsed group names, this could include, e.g. a:b
  unparsedGroupNames <- names(lformula$reTrms$cnms)
  # a function to parse group names
  groupParser <- function(groupi) {
    # have base::all.vars parse the group name
    all.vars(formula(paste0("~", groupi)))
  }
  
  # apply the parser to each random effect, and unique them
  groupNames <- rev(unique(unlist(lapply(unparsedGroupNames, groupParser))))
  # reorder data by groups (in order to make Z matrix obtained from lme easier to work with)
  data <- data[do.call(order, lapply(rev(groupNames), function(colN) data[ , colN])), ]
  
  data <- do_center_group(center_group, data, groupNames, weights0)
  data <- do_center_grand(center_grand, data, weights0)

  # remove row names so that resorted order is used in lme model 
  row.names(data) <- NULL
  
  # run lmer to get a ballpark fit and model structure information
  lme <- fit_unweighted_model(formula, data, family, verbose)

  # get y, test for validity
  mf <- model.frame(lme)
  responseCol <- attributes(attributes(mf)$terms)$response
  y <- as.numeric(mf[ , responseCol])
  if(!is.null(family) && family$family == "binomial") {
    if(length(unique(y)) == 2) {
      y <- ifelse(y == min(y), 0, 1)
    }
    if(any(!y %in% c(0,1))) {
      bady <- unique(y)
      bady <- bady[1:(min(length(bady), 5))]
      stop("For a binomial model the outcomes must be 0 or 1. Examples of values found:", paste(bady, collapse=", "))
    }
  }
  # Get the Z (random effects) matrix from LME 
  model_matrix <- getME(lme,"mmList")
  
  z_groups <- names(model_matrix) #get full names random effects whcih include both variable and group level
  
  #now capture interactions
  all_groups <- names(summary(lme)$ngrps)
  groupNames <- all_groups

  # create columns for any interactions coming from the formula
  # this will be mainly used in 3 level models and is included for forward compatability 
  missingGroupVars <- all_groups[!all_groups %in% names(data)] #find which group names are not in data set (ie composite groups)

  # drop levels of factors in present vars
  presentVars <- all_groups[all_groups %in% names(data)]
  for(i in seq_along(presentVars)) {
    if(inherits(data[ , presentVars[i]], "factor")) {
      data[ , presentVars[i]] <- droplevels(data[, presentVars[i]])
    }
  }

  # the lowest level, assumed to be the grouping var,
  # but updated when not in "if (length(missingGroupVars)>0){" loop
  all_groups_lowest_level <- all_groups

  # for each missing group var
  for(mgi in seq_along(missingGroupVars)) {
    #transform interactions into thier componenet parts
    vars <- rownames(attr(terms.formula(as.formula(paste(". ~", paste(missingGroupVars[mgi], collapse="+"))) ), "factors"))[-1]
    # drop levels on vars in missing groups
    for(i in seq_along(vars)) {
      if(inherits(data[, vars[i]], "factor")) {
        data[, vars[i]] <- droplevels(data[, vars[i]])
      }
    }
    data[ , missingGroupVars[mgi]] <- apply(data[ , vars], 1, function(x) {
                                        paste(x, collapse=":")
                                      }) 
    # this is a composite term, use the final 
    vtab <- lapply(vars, function(x) {
              tab <- table(data[ , x])
              sum(tab>0)
            })
    all_groups_lowest_level[all_groups_lowest_level == all_groups[mgi]] <- vars[which.max(unlist(vtab))]
  }

  # prepare Z matrices for each level 
  Z <- list(NULL)
  # Z is de-duplicated. ZFull is always the size of the outcome vector

  ZFull <- list(NULL)
  n_rows <- nrow(data)
  for (i in 1:length(all_groups)){
    z_to_merge <- grepl(paste0("[|]\\W", all_groups[i], "$"), z_groups)
    Z_i <- matrix( unlist(model_matrix[z_to_merge], use.names=FALSE), nrow=n_rows)
    ZFull <- c(ZFull, list(Z_i))
    if(i > 1) {
      Z_i <- Z_i[!duplicated(data[ , all_groups[i-1]]), , drop=FALSE]
    }
    Z <- c(Z, list(Z_i))
  }

  # find number of random effects classes
  nz <- list(0) # there are not REs at level 1
  for(i in 1:length(Z)) {
    if(!is.null(Z[[i]])) {
      nz[[i]] <- ncol(Z[[i]])
    }
  }
  #calculate n levels from Z and warn if not the same dimension as weights 
  levels <- length(Z)
  if(length(weights) != levels) {
    stop(paste0("The argument ", dQuote("weights"),
                " must be a vector of column names with length equal to levels. length ",
                dQuote("weights"), " is ", length(weights), ", expecting ", levels))
  }
  
  # store the full sample weights (or cWeights) in wgts0
  # we also want conditional weights
  wgts0 <- wgtsC <- data[ , weights]
  if(cWeights) {
    for(i in (ncol(wgts0)-1):1) {
      wgts0[ , i] <- wgts0[ , i] * wgts0[ , i+1]
    }
  }else{
    for(i in (ncol(wgtsC)-1):1) {
      wgtsC[ , i] <- wgtsC[ , i]/wgtsC[ , i+1]
    }
  }
  
  # check if weights are potentially conditional
  # must happen after sorting so it is sorted correctly
  wgts0 <- getWgts0(data, weights, cWeights)
  wgtsC <- getWgtsC(data, weights, cWeights)
  
  nRE <- list()
  for(i in 1:(levels-1)){
    nRE[[i]] <- ncol(Z[[i+1]])
  }
  # transform weights into a list of  dataframes where each data frame has
  # columns representing unit weight at that level and index of group names 
  weights <- weights_cond <- list()
  for(i in 1:length(nz)) {
    df <- data.frame(w=unname(wgts0[ , i]), stringsAsFactors=FALSE)
    dfC <- data.frame(w=unname(wgtsC[ , i]), stringsAsFactors=FALSE)
    # add the index for this level, when there is one
    if(i < length(nz)) {
      df$indexp1 <- dfC$indexp1 <- data[ , all_groups[i]]
      
    }
    if(i > 1) {
      # for levels >1 data frame indexed by group name and values represent weights 
      df$index <- dfC$index <- data[ , all_groups[i-1]]
      agg <- aggregate(w ~ index, data=df, FUN=rvar)
      aggC <- aggregate(w ~ index, data=dfC, FUN=rvar)
      if(any(agg$w > sqrt(.Machine$double.eps))) {
        # filter to just non-zero SD cases
        agg <- agg[agg$w > sqrt(.Machine$double.eps), ]
        colnames(agg) <- c("index", "Std. Dev.")
        message("Cases with non-zero standard deviation of group weight, within group:")
        print(agg, row.names=FALSE)
        stop(paste0("Some level-", i, " weights vary within group."))
      }
      df <- df[!duplicated(df$index), ]
      dfC <- dfC[!duplicated(dfC$index), ]
    }
    weights[[i]] <- df
    weights_cond[[i]] <- dfC
  }



  #get the y variable name from the formula 
  y_label <- as.character(formula[[2]])
  
  # get lmer fixed effects and covariance terms
  k <- length(lmeb <- getME(lme, "fixef"))
  parlme <- c(lmeb)
  lmesummary <- summary(lme)
  # find number of unique groups 
  ngrp <- lmesummary$ngrps
  if(length(unique(ngrp)) != length(ngrp)) {
    # if non-nested groups are present (ie not all obs at lower levels are part of upper level groups) 
    # then there will be non unique entries in ngrp 
    message("Groups n-sizes:")
    print(ngrp)
    stop("This does not appear to be a nested model. Some levels of this model have the same number of subject/groups as the level above them.")
  }
  ngrpW <- lapply(weights, function(wdf) {
    return(list(mean=mean(wdf$w), sum=sum(wdf$w), min=min(wdf$w), max=max(wdf$w))) 
  })

  # set up variance and coefficients 
  lmeVarDF <- data.frame(lmesummary$varcor)
  parlme <- c(parlme, lmeVarDF$vcov)

  # use initial values for coefficients if they are provied, otherwise default to lme4
  if(is.null(start)) {
    est0 <- parlme
  } else {
    if(length(start) != length(parlme)) {
      stop(paste0("Expecting argument ", sQuote("start"), " to have ", length(est0), " elements, found ", length (start), " elements."))
    }
    est0 <- start
    names(est0) <- names(parlme)
  } 
  
  # prepare variance data frame
  # rename groups in var-cov matrix to be all distinct 
  ind <- 1
  while(sum(grepl(paste0("\\.", ind, "$"), lmeVarDF$grp)) > 0) {
    lmeVarDF$grp <- sub(paste0("\\.", ind, "$"), "", lmeVarDF$grp)
    ind <- ind + 1
  }
  lmeVarDF$sdcor <- NULL # remove extraneous column
  lmeVarDF$ngrp <- NA
  lmeVarDF$grp <- gsub(".", ":", lmeVarDF$grp, fixed=TRUE)
  
  # add column showing how many elements are in each group
  for(vari in 1:nrow(lmeVarDF)) {
    if(lmeVarDF$grp[vari] == "Residual") {
      lmeVarDF$ngrp[vari] <- nrow(data)
    } else {
      lmeVarDF$ngrp[vari] <- ngrp[names(ngrp) == lmeVarDF$grp[vari]]
    }
  }
  
  # add column indicating which model level each group belongs to
  ngrp2 <- rev(sort(unique(lmeVarDF$ngrp)))
  for(grpi in 1:length(ngrp2)) {
    lmeVarDF$level[ngrp2[grpi] == lmeVarDF$ngrp] <- grpi + ifelse("Residual" %in% lmeVarDF$grp, 0, 1) # add 1 if no residual
  }
  varCorrect <- is.na(lmeVarDF$var2) & lmeVarDF$vcov < 1
  if(any(varCorrect)) {
    lmeVarDF$vcov[varCorrect] <- pmax(log(lmeVarDF$vcov[varCorrect]) + 1, -3.59)
  }
  # add names to the variance terms of the paramter vector
  names(est0)[-(1:k)] <- lmeVarDF$grp
  
  #add full variable name to lmeVarDF for later use 
  lmeVarDF$fullGroup <- paste0(lmeVarDF$grp, ifelse(!is.na(lmeVarDF$var1), paste0(".", lmeVarDF$var1), ""))

  # use helper funtion to create a covariance matrix from the data frame with variance and covariance information
  covarianceConstructor <- covMat2Cov(lmeVarDF)


  # C is the realization of the covariance constructor
  C <- covarianceConstructor(est0[-(1:k)])
  # these are the realized y/X vector/matrix
  X <- getME(lme, "X")
  
  ##############################################################
  #   2a) Pass linear models to symbolic integration method    #
  ##############################################################
  if(is.null(family)){
    #get the raw Z matrix
    Z <- getME(lme, "Z")
    #extract Zt from lme and transpose to Z  
    temp_Z <- getME(lme, "Ztlist")
    #find out which level each applies to 
    z_levels <- unique(lmeVarDF[lmeVarDF$fullGroup%in%names(temp_Z),c("fullGroup","level")])
    Zlist <- list()
    for (i in 2:levels){
      z_names <- z_levels[z_levels$level==i,"fullGroup"]
      Zlist[[i-1]] <- Matrix::t(do.call(rbind, temp_Z[z_names]))
      #Zlist[[i-1]] <- t(as.matrix((Reduce(rbind, temp_Z[z_names]))))
    }
    #seperating into single random effects using the pointers
    pointers <- getME(lme, "Gp")

    #find levels at which z applies
    grp_level <- lmeVarDF$level
    names(grp_level) <- lmeVarDF$grp
    
    ref_comps <- names(getME(lme, "cnms"))
    Zlevels <- unique(grp_level[ref_comps])
    
   
    # group id needs to be incremetally increasing starting at 1
    # as factor then as numeric should take care of this; also store a crosswalk
    group_id_list <- lapply(all_groups, FUN=function(x){
      res <- data.frame(data[,x], as.numeric(as.factor(data[,x])))
      colnames(res) <- c(x, "index")
      res
    })
    # this as one big matrix
    group_id <- do.call(cbind, group_id_list)
    cn <- c()
    # give group_id good names; necessary for : and / to work
    names(all_groups) <- make.names(all_groups)
    for(i in 1:length(all_groups)) {
      cn <- c(cn, all_groups[i], paste0(all_groups[i], "_index"))
    }
    colnames(group_id) <- cn
    group_id <- group_id[ , c(all_groups, paste0(all_groups, "_index"))]
    weights_list <- lapply(1:length(weights), FUN=function(wi) {
      if(wi == 1) {
        # sort by incoming index, so no resorting
        return(weights[[wi]]$w)
      }
      # not lowest level, sort by numeric(factor())
      x <- weights[[wi]]
      x <- x[order( as.numeric(as.factor(x$index)) ), ]
      x$w
    })
    
    # if level >2 then set up conditional weigths
    weights_list_cond <- weights_list
    if(levels > 2 ){
      cWeights <- cbind(group_id, wgts0)
      for (level in 1:(levels-1)){
        cWeights[ , weights0[level]] <- cWeights[ , weights0[level]] / cWeights[ , weights0[level + 1]]
      }
      weights_list_cond[[1]] <- cWeights[ , weights0[1]] #first level weights dont get grouped 
      for (level in 2:levels){
        # grab the first element, sorted as the data came
        weights_list_cond[[level]] <- cWeights[!duplicated(cWeights[,all_groups[level-1]]), weights0[level]]
      }
    }
    theta <- getME(lme, "theta")
    theta1 <- theta
    for(i in 1:length(theta1)) {
      theta1[i] <- 1
    }
    group_id <- group_id[ , c(paste0(all_groups, "_index"), all_groups)]
    bsqG <- devG(y, X, Zlist=Zlist, Zlevels=Zlevels, weights=weights_list, weightsC = weights_list_cond,
                 groupID = group_id,
                 lmeVarDF = lmeVarDF,
                 v0=theta1)
    if(verbose) {
      message("Fitting weighted model.")
    }
    opt <- bobyqa(fn=bsqG, par=theta)
    names(opt$par) <- names(theta)
    
    #opt$par are the theta estimates
    bsq <- bsqG(opt$par, getBS=TRUE) 

    if(verbose) {
      message("Estimating covariance.")
    }
    bhatq <- bsq(opt$par, robustSE=TRUE) #calculate robust SE
    uMatList <- makeUMatList(bhatq, Zlist, theta)
    b2 <- function(f, optpar, b, sigma0, inds) {
      function(x) {
        sigma <- x[length(x)]
        x <- x[-length(x)]
        xp <- optpar
        xp[inds] <- x
        names(xp) <- names(optpar)
        f(v=xp, sigma=sigma, beta=b)$lnl
      }
    }
    varDF <- lmeVarDF[,c("grp", "var1", "var2", "vcov", "ngrp", "level")]
    varVC <- list(Residual=bhatq$sigma^2)
    varDF$vcov <- 0
    varDF$SEvcov <- NA
    #set up list jaccobian of gradient for post hoc wald test
    j_mat_list <- list() #set up to hold jaccobian 
    for(li in 2:levels) {
      iDelta <- bhatq$iDelta[[li]]
      iDeltai <- bhatq$sigma^2 * (iDelta %*% t(iDelta))
      # varDF for just this level, this is a temp df that is used just in this loop and then not stored
      varDFi <- varDF[varDF$level %in% c(li,1),]
      # names of the thetas, minus the sigma term which we add back manually when needed
      thetaNamesi <- ifelse(is.na(varDFi$var2), paste0(varDFi$grp,".", varDFi$var1), paste0(varDFi$grp, ".", varDFi$var2, ".", varDFi$var1))[-nrow(varDFi)]
      inds <- names(opt$par) %in% thetaNamesi
      ihes <- -1*getHessian(b2(f=bsq, optpar=opt$par, b=bhatq$b, sigma0=bhatq$sigma, inds=inds),
                                           x=c(opt$par[inds], sigma=bhatq$sigma))
      eihes <- eigen(ihes)
      if(min(eihes$values) <= 0 || max(eihes$values)/min(eihes$values) >= 1/((.Machine$double.eps)^0.25)) {
        warning("Numerical instability in estimating the standard error of variance terms. Consider the variance term standard errors approximate.")
        ihes <- nearPD(ihes,  posd.tol=400*sqrt(.Machine$double.eps))$mat
      }
      theta_cov_mat <- solve(ihes)
      colnames(theta_cov_mat) <- rownames(theta_cov_mat) <- c(names(opt$par[inds]),"sigma")
      J <- bhatq$Jacobian[rownames(theta_cov_mat), colnames(theta_cov_mat)]
      preVCi <- theta_cov_mat %*% J %*% theta_cov_mat
      # reorder to be in the same order as varDFi
      preVCi <- preVCi[c(thetaNamesi, "sigma"), c(thetaNamesi, "sigma")]
      cn <- colnames(iDeltai)
      sigma2 <- bhatq$sigma^2
      j_list <-  list()
      for(ii in 1:nrow(iDeltai)) {
        for(jj in ii:ncol(iDeltai)) {
          varDFi$grad <- 0
          if(ii==jj) {
            # diagonal element, so it is a variance
            # at this level, var1 is the column, and var2 is NA means it is a variance
            varDF[varDF$level==li & varDF$var1==cn[ii] & is.na(varDF$var2),"vcov"] <- iDeltai[ii,ii]
            # claculate the variance of the variance. let g = grad(iDeltai[ii,ii]) wrt the elements or iDelta and sigma2
            # diagonal element
            varDFi$grad[varDFi$var1 %in% rownames(iDelta)[ii] & is.na(varDFi$var2)] <- sigma2 * 2 * iDelta[ii,ii]
            # off-diagonal elements are needed when ii > 1
            if(ii > 1){
              for(iii in 1:(ii-1)) {
                varDFi$grad[(varDFi$var1 %in% rownames(iDelta)[ii] | varDFi$var2 %in% rownames(iDelta)[ii]) & (varDFi$var1 %in% rownames(iDelta)[iii] | varDFi$var2 %in% rownames(iDelta)[iii])] <- sigma2 * 2 * iDelta[ii,iii]
              }
            }
            # part of grad that is sigma
            varDFi$grad[nrow(varDFi)] <- 2 * iDeltai[ii,ii]/sqrt(sigma2)
            # this is the Taylor series gT %*% [VC of iDelta] %*% g
            varDF[varDF$level==li & varDF$var1==cn[ii] & is.na(varDF$var2),"SEvcov"] <- sqrt(t(varDFi$grad) %*% preVCi %*% varDFi$grad)
            j_list <- c(j_list,list(varDFi$grad))
          } else {
            # off-diagonal element, so it is a covariance
            varDF[varDF$level %in% li & varDF$var1 %in% cn[ii] & varDF$var2 %in% cn[jj],"vcov"] <- iDeltai[ii,jj]
            varDF[varDF$level %in% li & varDF$var1 %in% cn[jj] & varDF$var2 %in% cn[ii],"vcov"] <- iDeltai[ii,jj]
            # only if we calculate this term, this could also be iDeltai[ii,jj] == 0, but that could happen when fit
            if(any(varDF$level==li & (( varDF$var1==cn[ii] & varDF$var2 %in% cn[jj]) | (varDF$var1==cn[jj] & varDF$var2 %in% cn[ii])))) {
              # claculate the variance of the variance. let g = grad(iDeltai[ii,ii]) wrt the elements or iDelta and sigma2
              # here the term is sigma2 * sum_{k=1}^{min(ii,jj)} (iDelta[ii,k] * iDelta[jj,k])
              for(iii in 1:min(ii, jj)) {
                # the derivative wrt iDelta[ii,k] is sigma2 * iDelta[jj,k]
                # since diagonal elements are on varDFi with var2=NA, not var2=var1, they need special handeling
                if(ii == iii) {
                  varDFi$grad[(varDFi$var1 %in% rownames(iDelta)[ii] | varDFi$var2 %in% rownames(iDelta)[ii]) & is.na(varDFi$var2)] <- sigma2 * iDelta[jj,iii]
                } else {
                  varDFi$grad[(varDFi$var1 %in% rownames(iDelta)[ii] | varDFi$var2 %in% rownames(iDelta)[ii]) & (varDFi$var1 %in% rownames(iDelta)[iii] | varDFi$var2 %in% rownames(iDelta)[iii])] <- sigma2 * iDelta[jj,iii]
                }
                # the derivative wrt iDelta[jj,k] is sigma2 * iDelta[ii,k]
                if(jj == iii) {
                  varDFi$grad[(varDFi$var1 %in% rownames(iDelta)[jj] | varDFi$var2 %in% rownames(iDelta)[jj]) & is.na(varDFi$var2)] <- sigma2 * iDelta[ii,iii]
                } else {
                  varDFi$grad[(varDFi$var1 %in% rownames(iDelta)[jj] | varDFi$var2 %in% rownames(iDelta)[jj]) & (varDFi$var1 %in% rownames(iDelta)[iii] | varDFi$var2 %in% rownames(iDelta)[iii])] <- sigma2 * iDelta[ii,iii]
                } 
              }
              # part of grad that is sigma
              varDFi$grad[nrow(varDFi)] <- 2 * iDeltai[ii,jj]/sqrt(sigma2)
              # this is the Taylor series gT %*% [VC of iDelta] %*% g
              varDF[varDF$level==li & (( varDF$var1==cn[ii] & varDF$var2 %in% cn[jj]) | (varDF$var1==cn[jj] & varDF$var2 %in% cn[ii])),"SEvcov"] <- sqrt(t(varDFi$grad) %*% preVCi %*% varDFi$grad)
              #build Jaccobian of gradient for use in post hoc wald test 
              j_list <- c(j_list,list(varDFi$grad))
            }
          }
        }
      }
      
      jacobian <- matrix(unlist(j_list),ncol=length(j_list[[1]]),byrow=T)

      var_mat_var <-  jacobian %*% preVCi %*% t(jacobian)
      # add names - theta has names in   correct order so subset based on theta names also in this level
      rownames(var_mat_var) <-colnames(var_mat_var)   <-  names(theta)[names(theta) %in% rownames(preVCi)]
      j_list <-  list()  #re set for next level
      j_mat_list[[li-1]] <- var_mat_var

      # one could calculate var(residual) at every level. Prefer level-2
      if(li==2) {
        varDF[varDF$level==1,"SEvcov"] <- sqrt((2*sqrt(sigma2))^2*preVCi[nrow(preVCi),ncol(preVCi)])
      }
      varVC <- c(varVC, list(iDeltai))
      names(varVC)[li] <- (varDF$grp[varDF$level %in% li])[1]
    }
    var_of_var <- bdiag(j_mat_list)

    rownames(var_of_var) <- colnames(var_of_var)  <- unlist(sapply(j_mat_list,FUN=rownames))
    #get names from names of j_mat_list 
    varDF$vcov[varDF$grp=="Residual"] <- bhatq$sigma^2
    varDF$fullGroup <- paste0(varDF$grp,ifelse(!is.na(varDF$var1),paste0(".",varDF$var1),""))
    
    vars <- varDF$vcov[is.na(varDF$var2)]
    names(vars) <- varDF$fullGroup[is.na(varDF$var2)]
    
    # other output stats
    nobs <- nrow(X)
    names(nobs) <- "Number of obs"
    ngroups <- c(nobs, ngrp)

    #Calculate ICC 
    var_between <- sum(varDF[which(!is.na(varDF$var1) & is.na(varDF$var2)),"vcov"])
    var_within <- varDF$vcov[varDF$grp=="Residual"]
    ICC <- var_between/(var_between+var_within)
    
    # for backwards compatibility with EdSurvey 2.2.2
    env <- environment(bsq)
    #covCon <- get("cConstructor", env)
    #lmeVarDf <- get("covMat", environment(covCon))
    covMat <- env$lmeVarDF
    cc <- function() {
    }
    assign("cConstructor", value=cc, envir=env)


    # now build the results
    res <-list(lnlf=bsq, lnl= bhatq$lnl, coef = bhatq$b, ranefs=bhatq$ranef,
               SE = bhatq$seBetaRobust, 
               vars= vars,
               theta=bhatq$theta, call=call,
               levels=levels, CMODE=bhatq$ranef,
               invHessian=bhatq$cov_mat, ICC=ICC,
               is_adaptive=FALSE, sigma=bhatq$sigma, cov_mat=bhatq$varBetaRobust,
               ngroups=ngroups, varDF=varDF, varVC=varVC,var_theta=var_of_var,
               wgtStats=ngrpW, ranefMat = uMatList)
    class(res) <- "WeMixResults"
    return(res)
  }
  
  ##############################################################
  # 2b) Identify integration parameter for non linear models  #
  ##############################################################
  
  if(verbose) {
    cat("Identifying initial integration locations estimates for random effects.\n")
  }
  
  n_l2 <- unique(lmeVarDF$ngrp[lmeVarDF$level == 2])
  n_top <- unique(lmeVarDF$ngrp[lmeVarDF$level == max(lmeVarDF$level)])

  l2_group_sizes <- as.vector(table(data[[unique(lmeVarDF$grp[lmeVarDF$level==2])]]))
  group_sizes <- as.vector(table(aggregate(.~data[[groupNames[1]]],data,unique)[[unique(lmeVarDF$grp[lmeVarDF$level==max(lmeVarDF$level)])]]))

  q2 <- nrow(lmeVarDF[lmeVarDF$level == 2 & is.na(lmeVarDF$var2),])
  q3 <- nrow(lmeVarDF[lmeVarDF$level == 3 & is.na(lmeVarDF$var2),])
  q_tot <- q3*n_top + q2*n_l2
  
  Z_mat <- getME(lme,"Z")
  lambda <- as(Cholesky(Matrix(constructSigma(exp(est0[-(1:ncol(X))]),nRE,n_l2,n_top))),"sparseMatrix")
  Whalf <- Diagonal(x = sqrt(1 / y))
  mu_eta <- Whalf
  PsiVec <- rep(weights[[2]]$w,each=q2)
  if(levels == 3){
    PsiVec <- c(PsiVec, rep(weights[[3]]$w,each=q3))
  }
  Psi <- Diagonal(x=PsiVec)
  
  #############################################
  #     3) Maximum Likelihood estimation      #
  #############################################
  
  # this gets used in optimization
  main_lnl <- main_lnl_container(q2,q3,X,y,Z_mat,lambda,Whalf,mu_eta,n_l2,n_top,
                                 l2_group_sizes,group_sizes,levels,weights,
                                 weights_cond,nQuad,family)
  
  # this gets returned to the user/used in summary methods
  main_lnlR <- lnl_by_group(q2,q3,X,y,Z_mat,lambda,Whalf,mu_eta,n_l2,n_top,
                            l2_group_sizes,group_sizes,levels,weights,
                            weights_cond,nQuad,family,parameterization = "var")
  tol <- 1e-8
  not_converged <- TRUE
  iteration <- 0
  est <- est0
  vcov.tmp <- est[-(1:ncol(X))]
  est <- c(est[1:ncol(X)],sqrt(abs(vcov.tmp)))

  while(iteration < max_iteration & not_converged){
    step_size <- 1
    
    if(iteration == 0){
      gradhess <- getGradHess(func=main_lnl, x=est)
    } else {
      gradhess <- gradhess_new
    }
    grad <- gradhess$grad
    hess <- -1 * gradhess$hess
    
    # if hessian isn't PD, transform it to be
    inc <- pt_inverse(hess)%*%grad
    
    going <- TRUE
    while(going & main_lnl(est) >
          main_lnl(est + step_size*inc)){
      step_size <- step_size * 0.5
      if(abs(step_size) < 0.001){
        if(step_size < 0) {
          going <- FALSE
        } else {
          step_size <- -1
        }
      }
    }
    
    est <- est + step_size*inc
    
    gradhess_new <- getGradHess(func=main_lnl, x=est)
    g <- gradhess_new$grad
    
    if(max(pmin(abs(est*g), abs(g))) < tol){
      not_converged <- FALSE
    }
    
    iteration <- iteration + 1
  }

  lambda <- updateLambda(lambda,est[-(1:ncol(X))],q2,q3,n_l2,n_top)
  Sigma <- lambda%*%Matrix::t(lambda)
  
  if(verbose) {
    message("Iterations complete.")
  }
  if (iteration >= max_iteration){
    #warning("Model exceeded maximum number of iterations and may not have converged.")
  }
  #############################################
  #            4) Post Estimation             #
  #############################################
  
  # need to add this
  tmp <- pirls_u(y,X,Z_mat,lambda,u=rep(0,length(PsiVec)),beta=est[1:ncol(X)],Whalf,mu_eta,
                 family,w1=weights[[1]]$w,PsiVec,Psi,l2_group_sizes,group_sizes,n_l2,n_top,
                 levels,q2,q3,by_group=FALSE)
  u_vector <- tmp$modes
  
  # make est numeric
  est.chol <- as.numeric(est)
  order <- vector(mode="numeric")
  if(q2==1){
    order <- c(1)
  }else{
    order <- c(1,3,2)
  }
  if(levels==3){
    if(q3 == 1){
      order <- c(order,max(order)+1)
    }else{
      order <- c(order,max(order)+1,max(order)+3,max(order)+2)
    } 
  }
  est[-(1:k)] <- unique(Sigma@x)[order]
  # make sure the names agree with lmer
  names(est) <- names(parlme)
  
  lnl_final <- main_lnl(est.chol)
  
  main_lnl <- main_lnl_container(q2,q3,X,y,Z_mat,lambda,Whalf,mu_eta,n_l2,n_top,
                                 l2_group_sizes,group_sizes,levels,weights,
                                 weights_cond,nQuad,family,parameterization = "var")

  hessian <- getHessian(main_lnl,est)
  # fix vars that are less than one to be mapped
  covs_and_vars <- est[-(1:k)]
  vars <- covs_and_vars[which(is.na(lmeVarDF$var2))] #select only vars not covars 
  need_fix_vars <- which(vars < 1)
  #covs_and_vars[need_fix_vars] <- exp(covs_and_vars[need_fix_vars] - 1)
  
  # get vars and name them
  vars <- covs_and_vars
  names(vars) <- gsub(":NA", "", paste(lmeVarDF$grp, lmeVarDF$var1, lmeVarDF$var2, sep=":"))

  # If any of variances got re mapped, re calculate hessian just inside the line, and 
  # share the limitation on WeMix estimation and the model wit the user.
  if (length(need_fix_vars) > 0){
    warning(paste0("Group variances too small to estimate accurately. The estimated variance in the group level terms(s) ", paste(dQuote(names(vars)[need_fix_vars]), collapse=", "), " is near zero.",
                   " Very low variance suggests that the data is not hierarchical and that a model without these levels should be considered.",
                   " If this removes all groups then a non-hierarchical model, such as logistic regression, should be considered."))
    # calculate hessian, moving covariances to just larger values so the Hessian is internal
    hessian <- getGradHess(main_lnl, c(est[1:k], covs_and_vars+0.0002*need_fix_vars))$hess
  }

  #calculate interclass corelation (ICC) 
  var_between <- sum(vars[which(!is.na(lmeVarDF$var1) & is.na(lmeVarDF$var2))])
  var_within <- vars[which(lmeVarDF$grp=="Residual")]
  ICC <- var_between/(var_between+var_within)
  
  nobs <- nrow(X)
  names(nobs) <- "Number of obs"
  ngroups <- c(nobs, ngrp)
  
  #set up the variance covariance matrix 
  varDF <- lmeVarDF[,c("grp", "var1", "var2", "vcov", "ngrp", "level")]
  
  varDF$vcov <- 0
  varDF$fullGroup <- paste0(varDF$grp,ifelse(!is.na(varDF$var1),paste0(".",varDF$var1),""))
  
  varDF$vcov <- vars #re assign in variance from mix.
  
  res <- list(lnlf=main_lnlR, lnl=lnl_final, coef=est[1:k], vars=vars,
              call=call, levels=levels, ICC=ICC, CMODE=u_vector,
              invHessian=hessian, is_adaptive=TRUE, ngroups=ngroups, varDF=varDF,
              wgtStats=ngrpW)
  class(res) <- "WeMixResults"
  return(res)
}


pirls_u <- function(y,X,Z_mat,lambda,u,beta,Whalf,mu_eta,family,w1,PsiVec,Psi,
                    l2_group_sizes,group_sizes,n_l2,n_top,levels,q2,q3,by_group=FALSE){
  
  useCholUpdate <- FALSE
  pirls <- TRUE
  iter <- 0
  Xb <- as.vector(X%*%beta)
  LtZt <- Matrix::t(lambda)%*%Matrix::t(Z_mat)
  u <- u*1
  u0 <- u
  eta <- Xb + as.vector(Matrix::crossprod(LtZt, u))
  mu <- family$linkinv(eta)
  
  L <- Cholesky(Matrix(Matrix::tcrossprod(LtZt)), perm=FALSE, LDL=FALSE, Imult=1)
  updatemu <- function(uu) {
    eta[] <<- as.numeric(Xb + as.vector(Matrix::crossprod(LtZt, uu)))
    mu[] <<- family$linkinv(eta)
    sum(family$dev.resids(y, mu, wt=w1)) + sum(PsiVec*uu^2)
  }

  olducden <- updatemu(u*0)
  
  while(pirls & iter <= 500){
    mu_eta@x <- family$mu.eta(eta)
    Whalf@x <- sqrt(w1 / (family$variance(mu)))
    # update weighted design matrix - this is U
    LtZtMWhalf <- as(LtZt%*%(mu_eta*diag(Whalf)),"sparseMatrix")
    
    # update Cholesky decomposition - if fitting unweighted model,
    # we can use the more efficient update function

    if(useCholUpdate){
      L <- update(L, LtZtMWhalf, 1) 
    }else{
      L <- Cholesky(Matrix::tcrossprod(LtZtMWhalf) + Psi, perm=FALSE, LDL=FALSE, Imult=0)
    }
    # update weighted residuals
    wtres <- diag(Whalf)*(y - mu)
    
    # solve for the increment
    delu <- as.vector(Matrix::solve(L, LtZtMWhalf%*%wtres - Psi%*%u))
    step_size <- 1
    if(iter == 0){
      # the first step can be overly large
      # if it remains a good idea, the second step can always take it
      step_size <- 0.01
    }
    
    ucden <- updatemu(u + delu)
    
    while(step_size > 0.001 && ucden > olducden){
      step_size <- step_size/2
      ucden <- updatemu(u + step_size*delu)
    }
    
    if(abs((olducden - ucden) / ucden) < 1e-8){
      pirls <- FALSE
    }
    
    olducden <- ucden
    u <- u + step_size*delu
    iter <- iter + 1
  }

  ucden <- updatemu(u)
  Whalf@x <- sqrt(w1 / family$variance(mu))
  mu_eta@x <- family$mu.eta(eta)
  LtZtMWhalf <- as(LtZt%*%(mu_eta*diag(Whalf)),"sparseMatrix")
  # update Cholesky decomposition
  if(useCholUpdate){
    L <- update(L,LtZtMWhalf,1)
  }else{
    L <- Cholesky(Matrix::tcrossprod(LtZtMWhalf) + Psi, perm=FALSE, LDL=FALSE, Imult=0)
  }
  L <- as(L,"sparseMatrix")

  llh <- -0.5*family$aic(y,rep.int(1,length(y)),mu,wt=w1,dev=NULL) - 0.5*sum(PsiVec*u^2) - 
    Matrix::determinant((L/sqrt(PsiVec))^PsiVec)$modulus
  
  det_by_group <- NULL
  # new args i need to add: l2_group_sizes,group_sizes,n_l2,n_top,levels
  if(by_group){

    llh <- llh_l2 <- l2_group_L <- l3_group_L <- vector(mode="numeric")
    log_y_u <- rowsum(family$lnl(y,mu,w1,sd=NULL),rep(1:n_l2,l2_group_sizes))

    for(i in 1:n_l2){
      group_Psi <- PsiVec[((i-1)*q2 + 1):(i*q2)]
      l2_group_L[i] <- Matrix::determinant(Matrix((L[((i-1)*q2 + 1):(i*q2),((i-1)*q2 + 1):(i*q2)]/sqrt(group_Psi))^group_Psi))$modulus
      group_u <- u[((i-1)*q2 + 1):(i*q2)]
      llh_l2[i] = log_y_u[i] - 0.5*sum(group_Psi*group_u^2) 
    }
    
    if(levels == 2){
      llh <- llh_l2 - l2_group_L
      det_by_group <- l2_group_L
    }else{
      l3_Psi <- PsiVec[((q2*n_l2)+1):length(PsiVec)]
      l2_Psi <- PsiVec[1:(q2*n_l2)]
      llh_tmp <- rowsum(llh_l2,rep(1:n_top,group_sizes))
      l3_u <- u[((q2*n_l2)+1):length(u)]
      l2_L <- L[1:(n_l2*q2),1:(n_l2*q2)]
      l3_L <- L[((n_l2*q2)+1):nrow(L),((n_l2*q2)+1):nrow(L)]
      l2_l3_L <- L[((n_l2*q2)+1):nrow(L),(1:n_l2*q2)]

      for(n in 1:n_top){
        l3_comp <- l3_L[((n-1)*q3 + 1):(n*q3),((n-1)*q3 + 1):(n*q3)]
        l2_comp <- l2_L[l2_l3_L[n,]!=0,l2_l3_L[n,]!=0]
        l2_l3_comp <- l2_l3_L[n,][l2_l3_L[n,]!=0]
        col.tmp <- matrix(0,ncol=q3,nrow=group_sizes[n])
        group_L <- Matrix(cbind(rbind(l2_comp,l2_l3_comp),rbind(col.tmp,l3_comp)))
        l3_group_Psi <- c(l2_Psi[l2_l3_L[n,]!=0],l3_Psi[((n-1)*q3 + 1):(n*q3)])
        l3_group_L[n] <- Matrix::determinant((group_L/sqrt(l3_group_Psi))^l3_group_Psi)$modulus
        l3_group_u <- l3_u[((n-1)*q3 + 1):(n*q3)]
        llh[n] <- llh_tmp[n] - 0.5*sum(l3_Psi[((n-1)*q3 + 1):(n*q3)]*l3_group_u^2) - l3_group_L[n]
        
      }
      det_by_group <- l3_group_L
    }
  }
  
  return(list(modes=u,hess=L,logLik=llh,det_by_group=det_by_group))
}

makeAdaptedPts <- function(muHat,R,nQuad,levels,q2,q3,n_l2,n_top){
  
  quadpts <- GHrule(nQuad)
  
  # extending univariate quadrature rule to dimension of random effects 
  # at each level
  # based on Pinheiro & Chao (2006) p. 71
  l2_points <- t(as.matrix(expand.grid(lapply(seq_len(q2),function(k,u) u[,1], u = quadpts))))
  l2_weights_tmp <- as.matrix(expand.grid(lapply(seq_len(q2),function(k,u) u[,2], u = quadpts)))
  l2_weights <- exp(rowSums(t(l2_points)*t(l2_points))/2)*apply(l2_weights_tmp,1,prod)
  
  muHat2 <- Matrix(muHat[1:(n_l2*q2)])
  
  R_l2 <- R[(1:(n_l2*q2)),(1:(n_l2*q2))]

  if(levels == 3){
    l3_points <- t(as.matrix(expand.grid(lapply(seq_len(q3),function(k,u) u[,1], u = quadpts))))
    l3_weights_tmp <- as.matrix(expand.grid(lapply(seq_len(q3),function(k,u) u[,2], u = quadpts)))
    l3_weights <- exp(rowSums(t(l3_points)*t(l3_points))/2)*apply(l3_weights_tmp,1,prod)
    muHat3 <- Matrix(muHat[(n_l2*q2+1):length(muHat)])
    
    R_l3 <- R[((n_l2*q2 +1):ncol(R)),((n_l2*q2 +1):ncol(R))]
    R_l2_l3 <- R[(1:(n_l2*q2)),((n_l2*q2 +1):ncol(R))]
    new_l3 <- Matrix::solve(R_l3,l3_points[rep(1:nrow(l3_points),times=n_top),]) + muHat3
    new_l3_expanded <- new_l3[,rep(1:ncol(new_l3),each=nQuad^q2)]
    
    new_l2 <- Matrix::solve(R_l2)%*%(l2_points[rep(1:nrow(l2_points),times=n_l2),
                                               rep(1:ncol(l2_points),times=nQuad^q3)] - 
                                       (R_l2_l3%*%Matrix::solve(R_l3)%*%l3_points[rep(1:nrow(l3_points),each=n_top),
                                                                                  rep(1:ncol(l3_points),each=nQuad^q2)])) + muHat2
    combined <- rbind(new_l2,new_l3_expanded)
    
    output <- list(adapt_l2 = new_l2,adapt_l3=new_l3,adapt_comb=combined,
                   l2_weights=l2_weights,l3_weights=l3_weights)
    
  }else{
    new_l2 <- Matrix::solve(R_l2,l2_points[rep(1:nrow(l2_points),times=n_l2),]) + muHat2
    combined <- new_l2
    output <- list(adapt_l2 = new_l2,l2_weights=l2_weights,adapt_comb=combined)
  }
  
  return(output)
}

constructSigma <- function(vcov,nRE,n_l2,n_top){
  
  mat_list <- list()
  
  # there will always be at least 2 levels
  q2 <- nRE[[1]]
  l2_var <- vcov[1:q2]
  l2_mat <- Diagonal(x=l2_var)
  l2_cov <- NULL
  if(q2 > 1){
    l2_cov <- vcov[(q2+1):(q2+choose(q2,2))]
    l2_mat[row(l2_mat) < col(l2_mat)] <- l2_mat[row(l2_mat) > col(l2_mat)] <- l2_cov
  }

  if(length(nRE) == 1){
    mat_list <- replicate(n_l2,l2_mat,simplify = FALSE)
    group_sigma <- bdiag(mat_list)
  }
  
  # for three levels, we need to create a block-diagonal matrix
  if(length(nRE) > 1){
    
    vcov <- vcov[-(1:(length(l2_var)+length(l2_cov)))]
    
    q3 <- nRE[[2]]
    l3_var <- vcov[1:q3]
    l3_mat <- Diagonal(x=l3_var)
    
    if(q3 > 1){
      l3_cov <- vcov[(q3+1):(q3+choose(q3,2))]
      l3_mat[row(l3_mat) < col(l3_mat)] <- l3_mat[row(l3_mat) > col(l3_mat)] <- l3_cov
    }

    mat_list <- c(replicate(n_l2,l2_mat,simplify = FALSE),replicate(n_top,l3_mat,simplify = FALSE))
    
    group_sigma <- bdiag(mat_list)
  }
  
  return(group_sigma)
  
}

updateLambda <- function(lambda,vcov,q2,q3,n_l2,n_top){

  if(q2 == 1){
    l2_theta <- rep(vcov[1],n_l2)
    if(length(vcov)>1){
      vcov_l3 <- vcov[2:length(vcov)]
    }
  }else{
    l2_theta <- rep(c(vcov[1],vcov[3],vcov[2]),n_l2)
    if(length(vcov)>3){
      vcov_l3 <- vcov[4:length(vcov)]
    }
  }
  
  if(q3 != 0){
    if(q3 == 1){
      l3_theta <- rep(vcov_l3[1],n_top)
    }else{
      l3_theta <- rep(c(vcov_l3[1],vcov_l3[3],vcov_l3[2]),n_top)
    }
    lambda@x <- c(l2_theta,l3_theta)
  }else{
    lambda@x <- l2_theta
  }
  
  lambda
}

pt_inverse <- function(mat,m=1e-5){
  e_dec <- eigen(mat)
  lambda <- diag(e_dec$values)
  Q <- e_dec$vectors
  for(i in 1:ncol(lambda)){
    if(abs(lambda[i,i]) >= m){
      lambda[i,i] <- abs(lambda[i,i])
    }else{
      lambda[i,i] <- m
    }
  }
  Q%*%solve(lambda)%*%t(Q)
}

main_lnl_container <- function(q2,q3,X,y,Z_mat,lambda,Whalf,mu_eta,n_l2,n_top,
                               l2_group_sizes,group_sizes,levels,weights,weights_cond,
                               nQuad,family,
                               parameterization="cholesky") {
  w1 <- weights[[1]]$w
  w1c <- weights_cond[[1]]$w
  w2 <- weights[[2]]$w
  w2c <- weights_cond[[2]]$w
  PsiVec <- rep(weights[[2]]$w,each=q2)
  if(levels == 3){
    PsiVec <- c(PsiVec, rep(weights[[3]]$w,each=q3))
  }
  Psi <- Diagonal(x=PsiVec)
  q_tot <- q3*n_top + q2*n_l2
  if(levels == 2){
    nRE <- list(q2)
  }else{
    nRE <- list(q2,q3)
  }
  old_u <- rep(0,q_tot)
  function(thetabeta){
    beta <- thetabeta[1:ncol(X)]
    vcov <- thetabeta[-(1:ncol(X))]
    if(parameterization != "cholesky"){
      sigma_mat <- as(nearPD(Matrix(constructSigma(vcov,nRE,n_l2,n_top)))$mat,"sparseMatrix")
      lambda <- as(Cholesky(sigma_mat),"sparseMatrix")
    }else{
      lambda <- updateLambda(lambda,vcov,q2,q3,n_l2,n_top)
    }

    postModeVar <- pirls_u(y,X,Z_mat,lambda,u=old_u,beta,Whalf,mu_eta,family,w1,PsiVec,Psi,
                           l2_group_sizes,group_sizes,n_l2,n_top,levels,q2,q3,by_group=FALSE)
    #old_u <<- 
    muHat <- postModeVar$modes
    
    L <- postModeVar$hess
    R <- Matrix::t(as(L,"sparseMatrix"))
    
    if(nQuad == 1){
      llh <- postModeVar$logLik
    } else {
      # get the adapted quadrature locations for levels 2 and 3
      adapted_points <- makeAdaptedPts(muHat,R/sqrt(PsiVec),nQuad,levels,q2,q3,n_l2,n_top)
      new_l2 <- adapted_points$adapt_l2
      l2_weights <- adapted_points$l2_weights
      Z_lambda <- Z_mat%*%lambda
      
      combined <- adapted_points$adapt_comb
      eta <- as.matrix(as.vector(X%*%beta) + Z_lambda%*%combined)
      mu <- family$linkinv(eta)
      
      # Don't do all the weighting here - do it once the point values have been
      # summed. Only weighting that should be done here is applying conditional
      # level 1 weights. 
      glm_llh <- family$lnl(y,mu,w1c,sd=NULL)
      norm_llh_l2 <- -0.5*rowsum(as.matrix(new_l2^2),group=rep(1:n_l2,each=q2))
      llh_l2_point <- rowsum(glm_llh,group=rep(1:n_l2,times=l2_group_sizes)) + norm_llh_l2
      
      # TO-DO: figure out more elegant way to do this - idea is summing up the point
      # values within a single level 3 group (i.e. sum up the 'l' level 2 points at
      # level 3 point 'k')
      
      if(levels == 3){
        w3 <- weights[[3]]$w
        new_l3 <- adapted_points$adapt_l3
        l3_weights <- adapted_points$l3_weights
        
        llh_l2_inner_sum <- vector()

        for(i in 1:nQuad^q3){
          llh_l2_inner_sum <- cbind(llh_l2_inner_sum, 
                                    rowLogSumExps(llh_l2_point[,((i-1)*(nQuad^q2) + 1):((nQuad^q2)*i)] + 
                                                    rep(log(l2_weights),each=n_l2)))
        }

        norm_llh_l3 <- -0.5*rowsum(as.matrix(new_l3^2),group=rep(1:n_top,each=q3))
        llh_l3_point <- rowsum(w2c*llh_l2_inner_sum,
                               group=rep(1:n_top,times=group_sizes)) + norm_llh_l3
        
        llh_l3 <- w3*rowLogSumExps(llh_l3_point + rep(log(l3_weights),each=n_top))
        
        llh <- sum(llh_l3) - Matrix::determinant((R/sqrt(PsiVec))^PsiVec)$modulus
      }else{
        l2_weight_mat <- matrix(rep(log(l2_weights),each=n_l2),nrow=n_l2)
        llh_l2 <- w2c*rowLogSumExps(llh_l2_point + l2_weight_mat)

        llh <- sum(llh_l2) - Matrix::determinant((R/sqrt(PsiVec))^PsiVec)$modulus
      }
      
    }
    llh
  }
}

lnl_by_group <- function(q2,q3,X,y,Z_mat,lambda,Whalf,mu_eta,n_l2,n_top,
                         l2_group_sizes,group_sizes,levels,weights,weights_cond,
                         nQuad,family,
                         parameterization="cholesky") {
  w1 <- weights[[1]]$w
  w1c <- weights_cond[[1]]$w
  w2 <- weights[[2]]$w
  w2c <- weights_cond[[2]]$w
  PsiVec <- rep(weights[[2]]$w,each=q2)
  if(levels == 3){
    PsiVec <- c(PsiVec, rep(weights[[3]]$w,each=q3))
  }
  Psi <- Diagonal(x=PsiVec)
  q_tot <- q3*n_top + q2*n_l2
  if(levels == 2){
    nRE <- list(q2)
  }else{
    nRE <- list(q2,q3)
  }
  old_u <- rep(0,q_tot)
  function(thetabeta,top=TRUE){
    beta <- thetabeta[1:ncol(X)]
    vcov <- thetabeta[-(1:ncol(X))]
    if(parameterization != "cholesky"){
      sigma_mat <- as(nearPD(Matrix(constructSigma(vcov,nRE,n_l2,n_top)))$mat,"sparseMatrix")
      lambda <- as(Cholesky(sigma_mat),"sparseMatrix")
    }else{
      lambda <- updateLambda(lambda,vcov,q2,q3,n_l2,n_top)
    }

    postModeVar <- pirls_u(y,X,Z_mat,lambda,u=old_u,beta,Whalf,mu_eta,family,w1,PsiVec,Psi,
                           l2_group_sizes,group_sizes,n_l2,n_top,levels,q2,q3,by_group=TRUE)
    #old_u <<- 
    muHat <- postModeVar$modes
    
    L <- postModeVar$hess
    R <- Matrix::t(as(L,"sparseMatrix"))
    
    if(nQuad == 1){
      llh <- postModeVar$logLik
    } else {
      # get the adapted quadrature locations for levels 2 and 3
      adapted_points <- makeAdaptedPts(muHat,R,nQuad,levels,q2,q3,n_l2,n_top)
      new_l2 <- adapted_points$adapt_l2
      l2_weights <- adapted_points$l2_weights
      Z_lambda <- Z_mat%*%lambda
      
      combined <- adapted_points$adapt_comb
      eta <- as.matrix(as.vector(X%*%beta) + Z_lambda%*%combined)
      mu <- family$linkinv(eta)
      # weights refers to the unit weights, which default to 1
      glm_llh <- family$lnl(y,mu,w1c,sd=NULL)
      norm_llh_l2 <- -0.5*rowsum(as.matrix(new_l2^2),group=rep(1:n_l2,each=q2))
      llh_l2_point <- rowsum(glm_llh,group=rep(1:n_l2,times=l2_group_sizes)) + norm_llh_l2
      
      # TO-DO: figure out more elegant way to do this - idea is summing up the point
      # values within a single level 3 group (i.e. sum up the 'l' level 2 points at
      # level 3 point 'k')
      
      if(levels == 3){
        w3 <- weights[[3]]$w
        new_l3 <- adapted_points$adapt_l3
        l3_weights <- adapted_points$l3_weights
        
        llh_l2_inner_sum <- vector()
        
        for(i in 1:nQuad^q3){
          llh_l2_inner_sum <- cbind(llh_l2_inner_sum, 
                                    rowLogSumExps(llh_l2_point[,((i-1)*(nQuad^q2) + 1):((nQuad^q2)*i)] + 
                                                    rep(log(l2_weights),each=n_l2)))
        }
        
        norm_llh_l3 <- -0.5*rowsum(as.matrix(new_l3^2),group=rep(1:n_top,each=q3))
        llh_l3_point <- rowsum(w2c*llh_l2_inner_sum,group=rep(1:n_top,times=group_sizes)) + norm_llh_l3
        
        llh_l3 <- w3*rowLogSumExps(llh_l3_point + rep(log(l3_weights),each=n_top))
        
        llh <- llh_l3 - postModeVar$det_by_group
      } else {
        llh_l2 <- w2*rowLogSumExps(llh_l2_point + rep(log(l2_weights),each=n_l2))
        llh <- llh_l2 - postModeVar$det_by_group
      }
      
    }
    if(top){
      sum(llh)
    }else{
      llh
    }
  } #ends  funciton to be return function(par, long=FALSE) 
} 


dropNonPositiveWeights <- function(data, weights) {
  #this removes any zero weight cases if they exist 
  dataW <- data[, weights, drop=FALSE]
  dataW[apply(dataW <= 0, 1, any), ] <- NA
  if(any(!complete.cases(dataW))) {
    warning(paste0("There were ", sum(complete.cases(dataW)==FALSE), " rows with non-positive weights. These have been removed."))
    data <- data[complete.cases(dataW), ]
  }
  return(data)
}

# return unconditional weights
getWgts0 <- function(data, weights, cWeights) {
  wgts0 <- data[ , weights]
  if(cWeights) {
    for(i in (ncol(wgts0)-1):1) {
      wgts0[ , i] <- wgts0[ , i] * wgts0[ , i+1]
    }
  }
  return(wgts0)
}

# get conditional weights
getWgtsC <- function(data, weights, cWeights) {
  wgtsC <- data[ , weights]
  if(!cWeights) {
    for(i in (ncol(wgtsC)-1):1) {
      wgtsC[ , i] <- wgtsC[ , i]/wgtsC[ , i+1]
    }
  }
  return(wgtsC)
}


sumsq <- function(x) {
  sum(x^2)
}

reweight <- function(data, level, method, data_m1=NULL) {
  if(length(method) > 1) {
    message("Length of method larger than one, using first method: ", dQuote(method), ".")
    method <- method[1]
  }
  if(!method %in% c("sample size", "effective sample size")) {
    stop("Unknown re-weighting method: ", dQuote(method), ".")
  }
  data$w0 <- data$w
  if("indexp1" %in% colnames(data)) {
    data$sumn <- ave(data$w, data$indexp1, FUN=length)
    data$sumw <- ave(data$w, data$indexp1, FUN=sum)
    data$sumw2 <- ave(data$w, data$indexp1, FUN=sumsq)
  } else {
    data$sumn <- ave(data$w, data$index, FUN=length)
    data$sumw <- ave(data$w, data$index, FUN=sum)
    data$sumw2 <- ave(data$w, data$index, FUN=sumsq)
  }
  if(method %in% "sample size") {
    data$w <- data$w0 * data$sumn / data$sumw
  }
  if(method %in% "effective sample size") {
    data$w <- data$w0 * data$sumw / data$sumw2
  }
  return(data)
}

do_center_group <- function(center_group, data, groupNames, weights0) {
  if(!is.null(center_group)) {
    #first add nested variables to data set if they exist, this is to handle / and : in group vars
    if (any(grep(":|/", names(center_group)))) {
      nested_groups <- names(center_group)[grep(":|/", names(center_group))]
      for (var in nested_groups){
        vars <- unlist(strsplit(var , ":|/"))
        data[,var] <- paste0(data[ , vars[1]], ":", data[ , vars[2]])
      }
    } #end of if there are : and / 
    if(!all(names(center_group) %in% names(data))){
      unfound_names <- names(center_group)[!names(center_group) %in% names(data)]
      stop("Not all centering group variables are found in the data set. Names include:", paste(unfound_names, collapse=", "))
    }
    if(length(center_group) != length(names(center_group))) {
      stop(paste0("The argument ", dQuote("center_group"), " must be a list where each element is the name of levels."))
    }
    for(name in names(center_group)) {
      # loop is included here for centering at >2 levels 
      # we need to get variable names from model matrix becasue of factor transformations
      # remove the first element becasue it is the intercept

      # identify level--it is the minimum group to account for : and / specified groups
      if(!inherits(center_group[[name]], "formula")) {
        stop("the ", dQuote("center_group"), " argument must be a list of formulas.")
      }
      X <- sparse.model.matrix(center_group[[name]], data=data)
      vars <- colnames(X)
      if(attr(terms(center_group[[name]]), "intercept") %in% 1) {
        vars <- vars[-1]
      }
      if(!all(vars %in% colnames(data))) {
        stop("could not find variable(s) in the ", dQuote("center_group"), " list element ", dQuote(name), " ", paste(vars[!vars %in% colnames(data)], collapse=", "))
      }
      if(!all(unlist(strsplit(name,":|/")) %in% groupNames)) {
        stop("the ", dQuote("center_group"), " argument must be a list with names that are group level names. Could not find level of ", dQuote(name), " in list of level grouping variables: ", paste(dQuote(unique(groupNames)), collapse=", "),".")
      }
      lev <- 1 # use level 1 weights
      # subtract group average from value and put back into data 
      # including scaling factor that accoutns for the fact weights might not sum to 0 l
      data[ , vars] <- sapply(vars, function(var){
                        # X - weighted group average X
                        data[ , var] - 
                          ave(X[ , var] * data[ , weights0[lev]], data[ , name], FUN=sum) /
                          ave(data[ , weights0[lev]], data[ , name], FUN=sum)
                      })
      rm(X)
    } # end for(name in names(center_group))
  } # end if(!is.null(center_group))
  return(data)
}

do_center_grand <- function(center_grand, data, weights0) {
  if(!is.null(center_grand)){
    X <- sparse.model.matrix(center_grand, data=data)
    vars <- colnames(X)
    if(attr(terms(center_grand), "intercept") %in% 1) {
      vars <- vars[-1]
    }
    #subtract overall average from value and put back into X 
    for(var in vars) {
      data[ , var] <- data[ , var] - sum(X[ , var] * data[ , weights0[1]]) / sum(data[ , weights0[1]])
    }
    rm(X)
  }
  return(data)
}

setup_family <- function(family) {
  # setup family
  # if family is set, use it
  if(!is.null(family)) {
    if(inherits(family, "character")) {
      family <- do.call(family, args=list())
    }
    if(!inherits(family, "family")) {
      stop(paste0("The family argument must be of class ", dQuote("family"), "."))
    }
    if(!family$family %in% c("binomial", "poisson", "gaussian")) {
      stop("Unimplemented family, ", dQuote(family$family))
    }
    family$lnl <- switch(family$family,
      binomial = function(y, mu, w, sd) {
                   w * dbinom(x=y, size=rep(1,length(y)), prob=mu, log=TRUE)
                 },
      poisson = function(y, mu, w, sd) {
                   w * dpois(x=y, lambda=mu, log=TRUE)
                 },
      gaussian = function(y, mu, w, sd) {
                   w * dnorm(x=y, mean=mu, sd=sd, log=TRUE)
                 }
    )
  }
  return(family)
}

fit_unweighted_model <- function(formula, data, family, verbose) {
  if(is.null(family)) {
    if(verbose) {
      cat("Using lmer to get an approximate (unweighted) estimate and model structure.\n")
    }
    # warnings happen on e.g. near-singular solve and are not a concern
    suppressWarnings(lme <- lmer(formula, data, REML=FALSE))
  } else {
    if(verbose) {
      cat("Using glmer to get an approximate (unweighted) estimate and model structure.\n")
    }
    lme <- glmer(formula, data, family=family)
  }
  return(lme)
}

# return 0 var if there is only one unit
#' @importFrom stats var
rvar <- function(x) {
  if(length(x) <=1) {
    return(0)
  } else {
    return(var(x))
  }
}
