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
#' @param acc0 integer, the precision of \code{mpfr}, default 120. Only  applies to non-linear models. 
#' @param start optional numeric vector representing the point at which the model should start optimization; takes the shape of c(coef, vars) 
#' from results (see help). 
#' @param family the family; optionally used to specify generalized linear mixed models. Currently only \code{binomial(link="logit")} is supported.
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
#' @importFrom lme4 getME lmer glmer lFormula
#' @importFrom stats dnorm aggregate terms dpois dgamma dbinom ave model.matrix terms.formula as.formula sigma complete.cases
#' @importFrom numDeriv genD hessian grad
#' @importFrom minqa bobyqa 
#' @importFrom Matrix nearPD sparse.model.matrix
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
#' @example \man\examples\mix.R
#' @author Paul Bailey, Claire Kelley, and Trang Nguyen 
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
  if(nQuad <= 0) stop(paste0("The argument ", sQuote("nQuad"), " must be a positive integer."))
  if(!inherits(run, "logical")) stop(paste0("The argument ", sQuote("run"), " must be a logical."))
  if(!inherits(verbose, "logical")) stop(paste0("The argument ", sQuote("verbose"), " must be a logical."))
  if(!inherits(weights, "character")) stop(paste0("The argument ", sQuote("weights"), " must be a character vector of weight column names in ", sQuote("data"), "."))
  if(any(!weights %in% colnames(data))) stop(paste0("The argument ", sQuote("weights"), " must specify valid columns in ", sQuote("data"), "."))
  if(acc0 <= 0) stop(paste0("The argument ", sQuote("acc0"), " must be a positive integer."))
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
  if(length(weights) == 1) {
    # if the length of weights is 1 then below subsets on weights do not work.
    stop(paste0("The argument ", sQuote("weights"), " must be a list of column names with length equal to levels."))
  }
  #this removes any zero weight cases if they exist 
  data[apply(data[ , weights] <= 0, 1, any), weights] <- NA
  if(any(is.na(data[ , weights]))) {
    warning(paste0("There were ", sum(complete.cases(data)==FALSE), " rows with non-positive weights. These have been removed."))
    data <- data[complete.cases(data), ]
  }
  
  # setup family
  # if family is set, use it
  if(!is.null(family)) {
    if(inherits(family, "character")) {
      family <- do.call(family, args=list())
    }
    if(!inherits(family, "family")) {
      stop(paste0("The family argument must be of class ", dQuote("family"), "."))
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
                 },
      Gamma = function(y, mu, w, sd) {
              stop("The gamma family is not implemented.")
                 },
      inverse.gaussian = function(y, mu, w, sd) {
                   stop("The inverse Gaussian family is not implemented.")
                 },
      function(y, mu, w, sd) {
        stop(paste0("Unknown family."))
      }
    )
  }

  # set up initial values 
  adapter <- "MAP" # initial function evaluation through MAP, BLUE estimator can be used post hoc
  weights0 <- weights # keep track of incomming weights column names
  # correctly format acc0
  acc0 <- round(acc0)
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
      stop("Not all centering group variables are found in the data set. ")
    } else {
      for(name in names(center_group)) {
        # loop is included here for centering at >2 levels 
        # we need to get variable names from model matrix becasue of factor transformations
        # remove the first element becasue it is the intercept

        # identify level--it is the minimum group to account for : and / specified groups
        lev <- min(which(groupNames %in% unlist(strsplit(name,":|/"))))
        X <- sparse.model.matrix(center_group[[name]],data=data)
        vars <- colnames(X)[-1]
        X <- cbind(X, data[ , weights0[lev]])
        colnames(X)[ncol(X)] <- weights0[lev]
        # subtract group average from value and put back into data 
        # including scaling factor that accoutns for the fact weights might not sum to 0 l
        data[ , vars] <- sapply(vars, function(var){
                           X[ , var] - 
                             ave(X[ , var] * X[ , weights0[lev]], data[ , name])/
                             (nrow(X)/sum(X[ , weights0[lev]]))
                         })
        # only used for centering
        rm(X)
      } # end for(name in names(center_group))
    } # end else for loop mean centering varibles 
  } # end if(!is.null(center_group))
   
  if(!is.null(center_grand)){
    X <- sparse.model.matrix(center_grand, data=data)
    vars <- colnames(X)[-1]
    
    #subtract overall avvrage from value and put back into X 
    data[ , vars] <- sapply(vars, function(var){X[ , var] - ave(X[ , var])})
    # only used for centering
    rm(X)
  }

  # remove row names so that resorted order is used in lme model 
  row.names(data) <- NULL
  
  # run lmer to get a ballpark fit and model structure information
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
  # get y, test for validity
  mf <- model.frame(lme)
  responseCol <- attributes(attributes(mf)$terms)$response
  y <- as.numeric(mf[ , responseCol])
  if(!is.null(family) && family$family == "binomial") {
    if(length(unique(y)) == 2) {
      y <- ifelse(y == min(y), 0, 1)
    }
    if(any(!y %in% c(0,1))) {
      stop("For a binomial model the outcomes must be 0 or 1.")
    }
  }
  # Get the Z (random effects) matrix from LME 
  model_matrix <- getME(lme,"mmList")
  
  z_groups <- names(model_matrix) #get full names random effects whcih include both variable and group level
  
  #now capture interactions
  all_groups <- names(summary(lme)$ngrps)
  groupNames <- all_groups  
  
  # store the full sample weights (or cWeights) in wgts0
  wgts0 <- data[ , weights]
  if(cWeights) {
    for(i in (ncol(wgts0)-1):1) {
      wgts0[ , i] <- wgts0[ , i] * wgts0[ , i+1]
    }
  }
  # check if weights are potentially conditional

  # create columns for any interactions coming from the formula
  # this will be mainly used in 3 level models and is included for forward compatability 
  missingGroupVars <- all_groups[!all_groups %in% names(data)] #find which group names are not in data set (ie composite groups)

  # drop levels of factors in present vars
  presentVars <- all_groups[all_groups %in% names(data)]
  for(i in seq_along(presentVars)) {
    if(inherits(data[, presentVars[i]], "factor")) {
      data[, presentVars[i]] <- droplevels(data[, presentVars[i]])
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
    stop(paste0("The argument ", sQuote("weights"), " must be a list of column names with length equal to levels."))  
  }
  
  # transform weights into a list of  dataframes where each data frame has
  # columns representing unit weight at that level and index of group names 
  weights <- list()
  for(i in 1:length(nz)) {
    df <- data.frame(w=unname(wgts0[ , i]), stringsAsFactors=FALSE)
    # add the index for this level, when there is one
    if(i < length(nz)) {
      df$indexp1 <- data[ , all_groups[i]]
    }
    if(i > 1) {
      # for levels >1 data frame indexed by group name and values represent weights 
      df$index <- data[ , all_groups[i-1]]
      # return 0 var if there is only one unit
      rvar <- function(x) {
        if(length(x) <=1) {
          return(0)
        } else {
          return(var(x))
        }
      }
      agg <- aggregate(w ~ index, data=df, FUN=rvar)
      if(any(agg$w > sqrt(.Machine$double.eps))) {
        stop(paste0("Some level-", i+1, " weights vary within group."))
      }
      df <- df[!duplicated(df$index), ]
    }
    weights[[i]] <- df
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
               wgtStats=ngrpW)
    class(res) <- "WeMixResults"
    return(res)
  }
  
  ##############################################################
  # 2b) Identify integration parameter for non linear models  #
  ##############################################################
  
  if(verbose) {
    cat("Identifying initial integration locations estimates for random effects.\n")
  }

  # setup these methods that, theoretically, can be used to find the 
  # adaptive integration points. For now, only MAP works.
  # Note: MAP finds the maximum of the likelihood surface.
  # it is not, properly, a MAP

  MAP0 <- MAP(groups=data[ , all_groups, drop=FALSE], y=y, X=X, levels=levels,
              Z=Z, ZFull=ZFull, weights=weights, k=k,
              qp=gauss.quad(nQuad, "hermite"),
              covariance_constructor=covarianceConstructor, verbose=verbose,
              nlmevar=nrow(lmeVarDF)-1, nz=nz, acc=acc0, family=family)
  # notes: BLUE finds the expected value of the likelihood surface
  # it is not, properly, a BLUE
  BLUE0 <- BLUE(groups=data[,all_groups,drop=FALSE], y=y, X=X, levels=levels,
                Z=Z, ZFull=ZFull, weights=weights, k=k,
                qp=gauss.quad(nQuad, "hermite"),
                covariance_constructor=covarianceConstructor, verbose=verbose,
                nlmevar=nrow(lmeVarDF)-1, nz=nz, acc=acc0, family=family)
  
  # form omega0
  # omega0 is a list of matricies where each list element represents a level
  # each matrix column represents the value of b for all units at the level 
  bvec <- getME(lme, "b") #get b values from lme results
  bvecCuts <- getME(lme, "Gp") #find out indices corresponding to different levels
  blist <- vector("list", levels) # NULL for level 1
  
  #startLoc helps transform the z vector into a z matrix, starting from the beginning ( index of 1)
  startLoc <- 1
  
  # number of rows that the Z matrix has at level (l-1), skipped 1 here
  # this is used exclusively to form bmat
  # the issue is that bvecCuts doesn't have cuts for factor levels
  # the n_rows_z code fixes thatS
  comps <- names(getME(lme,"cnms")) #has same n as bvec and tells you group and var
  n_rows_z <- list()
  for (i in 1:length(comps)){
    n_rows_z[i] <- lmeVarDF[lmeVarDF$grp == comps[i],"ngrp"][1]
  }

  #forming matrixes of b values 
  blist <- vector("list", levels) # NULl for level 1
  for(cuti in 2:length(bvecCuts)) {
    bmat <- matrix(bvec[startLoc:bvecCuts[cuti]], nrow=n_rows_z[[cuti-1]])
    li <- unique(lmeVarDF$level[lmeVarDF$ngrp==nrow(bmat)])
    blist[[li]] <- cbind(blist[[li]], bmat)
    startLoc <- bvecCuts[cuti] + 1
  }
  
  omega0 <- blist

  # make MAP estimate with b values and parameter estimates 
  a0 <- MAP0(omega0=omega0, par0=est0)

  # add determinant to scale z values
  # created the adapted quadrature locaiton by moving original points based on curvature Q
  # follows equation of Hartzel et al, pg 87 (Reference in main vignette)
  zScale <- lapply(a0$Qi0, function(Qi0i) {
                     if(is.null(Qi0i)) {
                       return(NULL)
                     }
                     df <- data.frame(detQ=sapply(Qi0i,det))
                     # add group labels
                     for(i in 1:length(groupNames)) {
                       if(length(unique(data[,groupNames[i]])) == nrow(df)) {
                         df[,groupNames[i]] <- unique(data[,groupNames[i]])
                         attr(df,"groups") <- c(attr(df, "groups"), groupNames[i])
                       }  
                     }
                     df
                   })
    
  index <- data.frame(data[,c(groupNames)])
  names(index) <- groupNames
  #add column with determinant of Q to the Z infomration (join by index which is group name)
  for(wi in 2:length(weights)) {
    Zgrps <- attr(zScale[[wi]], "groups")
    weights[[wi]] <- merge(weights[[wi]],zScale[[wi]][,c(Zgrps, "detQ")],by.x="index", by.y=Zgrps)
  }
  
  # first guess
  est <- est0

  # use lnl function to optimize
  qp <- gauss.quad(nQuad,"hermite")
  # the covariance constructor maps the variance terms that are less than 1
  # to exp(x-1) to avoid negative variance terms.
  # fn0 does not use this map and can be easily used by the user
  # fn0R does use that map and is intended for use by the optimizer
  # these calls are otherwise identical.
  fn0 <- param.lnl.quad(y=y,
                        X=X,
                        levels=levels,
                        Z=Z,
                        ZFull=ZFull,
                        Qi=a0$Qi,
                        QiFull=a0$QiFull,
                        omega=a0$omega,
                        omegaFull=a0$omegaFull,
                        W=weights,
                        k=k,
                        qp=qp,
                        cConstructor=covarianceConstructor,
                        acc0=acc0,
                        mappedDefault=FALSE,
                        family=family)
  
  fn0R <- param.lnl.quad(y=y,
                         X=X,
                         levels=levels,
                         Z=Z,
                         ZFull=ZFull,
                         Qi=a0$Qi,
                         QiFull=a0$QiFull,
                         omega=a0$omega,
                         omegaFull=a0$omegaFull,
                         W=weights,
                         k=k,
                         qp=qp,
                         cConstructor=covarianceConstructor,
                         acc0=acc0,
                         mappedDefault=TRUE, 
                         family=family)

  if(!run) {
    # if run option is false the function terminates here, returning the components used in the optimization
    # without performing the optimization
    return(list(lnlf=fn0R, parlme=parlme, omega0=a0$omega0, lme=lme, adapt=a0, weights=weights))
  }
  
  #############################################
  #     3) Maximum Likelihood estimation      #
  #############################################
  
  d1 <- rep(Inf, length(est)) #d1 represent first derivatives of estimates 
  oldlnl <- fn0(est, varFloor=-3.59) # evaluate the likelihood with the estimated coefficients before adaptive quadrature 
  
  a00 <- a0 # start from previous estimates (MAP estimates) 
    
  if(verbose) {
    cat("Starting Newton steps.\n")
  }
  
  #select just variances that are > -3 to aviod continued adaptation below 0 
  covs_and_vars <- est[-(1:k)]
  vars <- covs_and_vars[which(is.na(lmeVarDF$var2))]
  not_0_vars <- which(vars>-3)+k
  est[-(1:k)]  <- ifelse(est[-(1:k)]< -4.6,-4.6,est[-(1:k)])  #reset variance that get below threshold
  # v is the full Newton step vector. Here set to a stub value
  v <- d1

  # Newton steps continue until step size is under threshold, excluding variances < -3
  # the denominator prevents relative changes that are indistinguishable from 0 from being
  # included.
  skipNextHessian <- FALSE
  # every step is of all parts, this works best.
  # but leave this in here since we may later enounter a situation where it helps
  defStepsInds <- list(1:length(est0))
  stepIndQueue <- list()
  dd1 <- d1
  dd2 <- outer(dd1,dd1) 
  iteration <- 0
  varFloorBinding <- FALSE
  oldest <- est

  # first condition, a variance is 0.01, so we look just for fitted values to not change
  # second condition, look for derivatives to be small
  while(all(iteration < max_iteration,
            any(varFloorBinding & max(est - oldest) > 1e-5 ,
               !varFloorBinding & max(abs(dd1[c(1:k, not_0_vars)]/pmax(abs(est[c(1:k, not_0_vars)]), 1e-5))) > 1E-5)
                               )) { 
    iteration <- iteration + 1
    oldest <- est
    if(length(stepIndQueue)==0) {
      stepIndQueue <- defStepsInds
    }
    thisStepInds <- stepIndQueue[[1]]
    # calls helper function to get both first and second derivative 
    d1 <- getGrad(fn0, est, thisStepInds)
    # update full length grad
    dd1[thisStepInds] <- d1
    if(!skipNextHessian) {
      d2 <- getHessian(fn0, est, thisStepInds)
      # update full size Hessian
      dd2[thisStepInds, thisStepInds] <- d2
    }
    # use most recent Hessian
    d2 <- dd2[thisStepInds, thisStepInds]
    fact <- 1 # scaling factor by which this step is multiplied 
    v <- rep(0, length(est0))
    v[thisStepInds] <- solve(d2) %*% d1 # the Newton step
    if(verbose) {
      cat("step:", iteration, "/", max_iteration, "\n")
      cat("lnl:", oldlnl, " max (relative) derivative=", max(abs(dd1[c(1:k,not_0_vars)]/pmax(abs(est[c(1:k,not_0_vars)]), 1e-5))), " ")
      cat("\nCurrent solution, gradient, and Newton step:\n")
      prnt <- cbind(oldEstimate=est, firstDeriv=dd1, proposedNewtonEstimate=est - v)
      rownames(prnt) <- c(names(est0)[1:k], paste0("ln var ", names(est0)[-(1:k)], ""))
      colnames(prnt) <- c("previous Est", "firstDeriv", "Newton Step")
      print(prnt)
    }
    
    # create new estimate by stepping in the direction idicated by gradient and scaled by fact
    newest <- est - fact * v
    newlnl <- fn0(newest, varFloor=-3.59)
    stp <- 0
    
    # make sure Newton step is a step improves lnl
    while(newlnl < oldlnl) {
      stp <- stp + 1
      if(verbose) {
        cat("Halving step size.\n")
      }
      fact <- fact/2
      if(stp > 5 & fact > 0) {
        if(verbose) {
          cat("Reversing step direction.\n")
        }
        ##reverse direction if more than 5 steps have been taken, avoid newtons method getting stuck 
        fact <- -1
        stp <- 0
      }
      if (stp>10) {
        # bad search direction, do nothing
        fact <- 0
        oldlnl <- oldlnl - 1
      }
      newest <- est - fact * v
      newlnl <- fn0(newest, varFloor=-3.59)      
    } # end while(newlnl < oldlnl)
    if(verbose) {
      cat("\n")
    }
    # update current estimate with the step we decided on
    est <- est - fact * v
    # update the lnl
    oldlnl <- newlnl

    # now, worry about threshold
    if(any(est[-(1:k)] < -3.59)) {
      est[-(1:k)]  <- ifelse(est[-(1:k)] < -3.59, -3.59, est[-(1:k)])  #reset variance that get below threshold
      varFloorBinding <- TRUE
      # no need up update oldlnl, this varFloor is always enforced
    }
    #adapts new quadrature points if parameter keepAdapting is true  
    if(keepAdapting) {
      if(verbose) {
        cat("Adapting random effect estimates.\n")
      }
      # adapter is never BLUE now
      if(adapter == "BLUE") {
        a0 <- BLUE0(omega0=a0$omega0, par0=est0, Qi0=a0$Qi0)
      } else {
        # update the b and Q parameters
        a0 <- MAP0(omega0=a0$omega0,
                   par0=est,
                   verb=FALSE)
      }
      
      #recacluate scaled Z values based on new Q (Curvature)
      zScale <- lapply(a0$Qi0, function(Qi0i) {
                   if(is.null(Qi0i)) {
                     return(NULL)
                   }
                   df <- data.frame(detQ=sapply(Qi0i,det))
                   # add group labels
                   for(i in 1:length(groupNames)) {
                     if(length(unique(data[,groupNames[i]])) == nrow(df)) {
                       df[,groupNames[i]] <- unique(data[,groupNames[i]])
                       attr(df,"groups") <- c(attr(df, "groups"), groupNames[i])
                     }  
                   }
                   df
                 })
      
      # move new values of detQ on to dataframes in weights list 
      for(wi in 2:length(weights)) {
        weights[[wi]]$detQ <- NULL
        Zgrps <- attr(zScale[[wi]], "groups")
        weights[[wi]] <- merge(weights[[wi]],zScale[[wi]][,c(Zgrps, "detQ")],by.x="index", by.y=Zgrps)
      }

      # update function for unser and function for optimization
      fn0 <- param.lnl.quad(y=y,
                            X=X,
                            levels=levels,
                            Z=Z,
                            ZFull=ZFull,
                            Qi=a0$Qi,
                            QiFull=a0$QiFull,
                            omega=a0$omega,
                            omegaFull=a0$omegaFull,
                            W=weights,
                            k=k,
                            qp=qp,
                            cConstructor=covarianceConstructor,
                            acc0=acc0,
                            mappedDefault=FALSE,
                            family=family)
      
      fn0R <- param.lnl.quad(y=y,
                             X=X,
                             levels=levels,
                             Z=Z,
                             ZFull=ZFull,
                             Qi=a0$Qi,
                             QiFull=a0$QiFull,
                             omega=a0$omega,
                             omegaFull=a0$omegaFull,
                             W=weights,
                             k=k,
                             qp=qp,
                             cConstructor=covarianceConstructor,
                             acc0=acc0,
                             mappedDefault=TRUE,
                             family=family)
      
      # Process of adapting stops when there the relative difference between chosen points is below threshold 
      if(max(abs(a00$omega0[[2]] - a0$omega0[[2]])/pmax(abs(a0$omega0[[2]]),1E-10)) < 1E-2) {
        if(verbose) {
          cat("Done adapting; the mode is not changing sufficiently.\n")
        }
        keepAdapting <- FALSE
      }
      if(keepAdapting & max(abs(d1)) <= 1E-3) {
        if(verbose) {
          cat("Done adapting: close to a solution.\n")
        }
        keepAdapting <- FALSE
      }
      a00 <- a0
      # the lnl function may have shifted enough that we need to update the current lnl
      # so update it based on these new quad points
      oldlnl <- fn0(est, varFloor=-3.59)
    } # end if(keepAdapting)
    #end of the Newton's method while loop
    
    #re updates the index of vars that are greater than -3 based on new est values 
    covs_and_vars <- est[-(1:k)]
    vars <- covs_and_vars[which(is.na(lmeVarDF$var2))] #select only vars not covars 
    not_0_vars <- which(vars > -3) + k
    if((!skipNextHessian & max(sum(abs(fact*v)/abs(est))) < (.Machine$double.eps)^0.25) & max(abs(dd2)) < Inf) {
      skipNextHessian <- TRUE
    } else {
      skipNextHessian <- FALSE
    }
  } # end of while(max(abs(d1/pmax(abs(est),1e-5))) > 1E-5)
  if(verbose) {
    message("Itterations complete.")
  }
  if (iteration >= max_iteration){
    #warning("Model exceeded maximum number of iterations and may not have converged.")
  }
  #############################################
  #            4) Post Estimation             #
  #############################################
 
  hessian <- dd2
  
  # For each "random effect" calculate the maximum of the likelihood surface (MAP)
  # and the expected value of the likelihood surface (BLUE)
  MAP <- MAP0(omega0=a0$omega0, par0=est, verb=FALSE)$omega0
  BLUE <- BLUE0(omega0=a0$omega0, par0=est, Qi0=a0$Qi0, adapt=FALSE, verb=FALSE)
  
  # make est numeric
  est <- as.numeric(est)
  
  # make sure the names agree with lmer
  names(est) <- names(parlme)
  
  # fix vars that are less than one to be mapped
  covs_and_vars <- est[-(1:k)]
  vars <- covs_and_vars[which(is.na(lmeVarDF$var2))] #select only vars not covars 
  need_fix_vars <- which(vars < 1)
  covs_and_vars[need_fix_vars] <- exp(covs_and_vars[need_fix_vars] - 1)
  
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
    hessian <- getHessian(fn0R, c(est[1:k], covs_and_vars+0.0002*need_fix_vars))
  }

  #calculate interclass corelation (ICC) 
  var_between <- sum(vars[which(!is.na(lmeVarDF$var1) & is.na(lmeVarDF$var2))])
  var_within <- vars[which(lmeVarDF$grp=="Residual")]
  ICC <- var_between/(var_between+var_within)
  
  # it is possible for the lnl to have changed slightly, so update it to avoid confusion
  nobs <- nrow(X)
  names(nobs) <- "Number of obs"
  ngroups <- c(nobs, ngrp)
  
  #set up the variance covariance matrix 
  varDF <- lmeVarDF[,c("grp", "var1", "var2", "vcov", "ngrp", "level")]
  
  varDF$vcov <- 0
  varDF$fullGroup <- paste0(varDF$grp,ifelse(!is.na(varDF$var1),paste0(".",varDF$var1),""))
  
  varDF$vcov <- vars #re assign in variance from mix.  This works without formatting because only 2 level models are posisble
  
  res <- list(lnlf=fn0R, lnl=fn0(est, varFloor=-3.59), coef=est[1:k], vars=vars,
              call=call, levels=levels, ICC=ICC, CMODE=BLUE,
              invHessian=hessian, is_adaptive=TRUE, ngroups=ngroups, varDF=varDF,
              wgtStats=ngrpW)
  class(res) <- "WeMixResults"
  return(res)
}


# finds expected value of the "random effects" likelihood surfaces
# this function cannot find these without input estimates; it cannot be used for adapting.
# @author Paul Bailey
BLUE <- function(groups, y, X, levels, Z, ZFull, weights0, k, qp,
                 covariance_constructor, verbose, nlmevar, nz, acc,
                 family) {

  # must define one of Qi or Qi0. Defining Qi is faster
  #this is the function that is returned when BLUE is called 
  function(omega0, par0, Qi=NULL, Qi0=NULL, verb=verbose, adapt=TRUE) {
    
    # form the Q matrixes from weights
    weights <- weights0 #start with original weights 
    if(is.null(Qi)) {
      Qi <- list(NULL)
      QiFull <- list(NULL)
      for( oi in 2:length(omega0)) {
        map <- groups[,oi-1] # decrement index by one because the first groups data frame doesnt have a column for level 1
        umap <- unique(map)
        nzi <- ncol(Z[[oi]]) # find number of z variables at this level 
        Qi[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(weights[[oi-1]]))
        for(i in 1:nrow(weights[[oi-1]])) {
           #shape Qi, a list with one element per model level
          #each element of Qi has one row for each random effect (Z) and one column for each observation in the level below-random effect combination
          #values are being moved over from Qi0 which is the MAP estimates at the group level
          Qi[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
        }
        QiFull[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(X))
        for(i in 1:nrow(X)) {
          #Also  full version which as one row for each Z and one column for each obs-random effect combination
          QiFull[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
        }
      } # ends  for( oi in 2:length(omega0))
    } # ends if(is.null(Qi)

    # function that creates dataframe to be used to scale Z value by the scaling factor calculated as the determinate of the curvature Q
    zScale <- lapply(Qi0, function(Qi0i) {
                 if(is.null(Qi0i)) {
                   return(NULL)
                 }
                 df <- data.frame(detQ=sapply(Qi0i,det))
                 # add group labels 
                 for(i in 1:ncol(groups)) {
                   if(length(unique(groups[,i])) == nrow(df)) {
                     df[,colnames(groups)[i]] <- unique(groups[,i])
                     attr(df,"groups") <- c(attr(df, "groups"), colnames(groups)[i])
                   }
                 }
                 df
               })

    # for  each data frame in weightes merger on the scaling factor determinat of
    # curvature (detQ)
    for(wi in 2:length(weights)) {
      weights[[wi]]$detQ <- NULL
      Zgrps <- attr(zScale[[wi]], "groups")
      weights[[wi]] <- merge(weights[[wi]], zScale[[wi]][,c(Zgrps, "detQ")],by.x="index", by.y=Zgrps)
    }
     
    # make new omega
    omega <- buildOmega(omega0=omega0, groups=groups, nrowX=nrow(X))
    omegaFull <- buildOmega(omega0=omega0, groups=groups, nrowX=nrow(X), full=TRUE)
    # make stubs
    Qi0_ <- list(NULL)
    Qi_ <- list(NULL)
    tmpomega <- list(NULL)
    for( oi in 2:length(omega0)) {
      omg0 <- omega0[[oi]]
      omg1 <- 2*omg0  #increment up to allow while loop to run
       
      # keep looking for the mean
      while( max(abs( (omg1 - omg0) / pmax(abs(omg0), 1E-5))) > 1E-3) {
        omg1 <- omg0
        tmpomega_ <- c(tmpomega, list(omg0))
        nzi <- ncol(Z[[oi]]) # number of Z columns at olvl
        f <- param.lnl.quad(y, X, oi, Z, ZFull=ZFull, Qi=Qi, QiFull=QiFull,
                            omega, omegaFull=omegaFull, W=weights, k, qp,
                            covariance_constructor, bobyqa=FALSE, verbose=TRUE, acc0=acc,
                            mappedDefault=FALSE, family=family)
        for(ici in 1:ncol(omg0)) {
          f0 <- f(par0, top=FALSE, integralMultiplierExponent=0, integralZColumn=ici)
          f1 <- f(par0, top=FALSE, integralMultiplierExponent=1, integralZColumn=ici)
          # f1 is not logged, f0 is
          omg0[ , ici] <- as.numeric(f1/f0)
        }
        # make omega
        omega0p <- c(tmpomega, list(omg0))
        while( length(omega0p) < length(omega0)) {
          omega0p[[length(omega0p)+1]] <- omega0[[length(omega0p)+1]] 
        }
        omega <- buildOmega(omega0=omega0p, groups=groups, nrowX=nrow(X))
        omegaFull <- buildOmega(omega0=omega0p, groups=groups, nrowX=nrow(X), full=TRUE)
        if(!adapt) {
          # no need to loop, this is the BLUE
          omg1 <- omg0
        }
      }
      if(verb & adapt) {
        cat("BLUE estimates:\n")
        print(omg0)
      }
      # move the integration points to the new BLUE estiamtes and update
      if(adapt) {
        omg0Full <- buildOmega(omega0=tmpomega_, groups=groups, nrowX=nrow(X), full=TRUE)
        derivatives <- genD(adapterLnL(y, X, levels, Z, ZFull, weights, k, qp,
                                       covariance_constructor, omega,
                                       omg0Full,
                                       tmpomega_, par0, verb, Qi, QiFull, oi,
                                       acc, family),
                            rep(0,sum(unlist(nz)[1:oi], na.rm=TRUE)))
        d2 <- derivatives$D[,-(1:nzi),drop=FALSE] # remove first derivatives
        drv <- d2
        # assumes there are exactly two levels
        Qi0_[[oi]] <- lapply(1:nrow(drv), function(i) { 
          scaleQuadPoints(drv[i,], nzi)
        })
        map <- groups[,oi-1]
        umap <- unique(map)
        Qi_[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(X))
        for(i in 1:nrow(X)) {
          Qi_[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0_[[oi]][[(1:length(umap))[map[i]==umap] ]]
        } 
        QiFull[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(X))
        for(i in 1:nrow(X)) {
          QiFull[[oi]][1:nz,(i-1)*nz+1:nz] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
        }
      }
      # now, actually add this to the tmpomega
      tmpomega <- c(tmpomega, list(omg0))
      omg0Full <- buildOmega(omega0=tmpomega, groups=groups, nrowX=nrow(X), full=TRUE)
    }
    if(adapt) {
      return(list(omega0=tmpomega, omega=omega, omegaFull=omg0Full, Qi0=Qi0_, Qi=Qi_, QiFull=QiFull))
    } else {
      return(tmpomega)
    }
  }
}


# finds the maximum of the likelihood surface and (approximate) SD, per "random effect"
# this function can find these without input estimates. This is used for adapting
MAP <- function(groups, y, X, levels, Z, ZFull, weights, k, qp,
                covariance_constructor, verbose, nlmevar, nz, acc, family) {
  ####################################################
  #                        Outline:                  #   
  #        1) Find points for MAP estimation         #
  #        2) MAP estimate                           #
  #         3) Estimate variance                     #
  ####################################################
  
  # when you call adapter, it returns this function that can be called with omega0 and par0
  # and will find a new MAP and SD for each RE
  function(omega0, par0, verb=verbose) {
    #####################################################
    #        1) Find points for MAP estimation         #
    #####################################################
    
    # make omega, which is omega0 at the level one  data length
    omega <- buildOmega(omega0=omega0, groups=groups, nrowX=nrow(X))
    omegaFull <- buildOmega(omega0=omega0, groups=groups, nrowX=nrow(X), full=TRUE)
   
    # f returns the log(likelihood * posterior), per group
    # use Newton's method, per group, to find the maximum of the a posterori surface
    Qi0 <- list(NULL)
    Qi <- list(NULL)
    QiFull <- list(NULL)
    tmpomega <- list(NULL)
    tmpomegaFull <- list(NULL)
    u0 <- 0
    for(oi in 2:length(omega0)) {
      # Over all levels >1 
      omg0 <- omega0[[oi]]
      omg1 <- 1E20*(omg0+1E-15)
      nzi <- nz[[oi]]
      iter <- 0
      while( iter < 25 & max(abs( (omg1 - omg0) / pmax(abs(omg0), 1E-5))) > 1E-3) { # iterate while Omega continues to change more than threshold
        if(iter >= 1) {
          u0 <- max(c(u0, as.vector(abs( (omg1 - omg0) ))))
        }
        iter <- iter + 1
        omg1 <- omg0
        tmpomega_ <- c(tmpomega, list(omg0))
        toF <- buildOmega(omega0=tmpomega_, groups=groups, nrowX=nrow(X), full=TRUE)
        ofn <- adapterLnL(y, X, levels, Z, ZFull, weights, k, qp,
                          covariance_constructor, omega, toF,
                          tmpomega_, par0, verb, Qi, QiFull, oi, acc,
                          family)
        # this is just second derivatives
        d1 <- getJacobian(ofn, rep(0, nz[[oi]], na.rm=TRUE), m=nrow(omg0))
        d2 <- getHessian(ofn, rep(0, nz[[oi]], na.rm=TRUE))
        
        ######################################
        #        2) MAP estimate            #
        ######################################
        
        # this is the Newton step, per group
        omg0 <- lapply(1:length(d2),
                       function(i) {
                         step <- solve(d2[[i]]) %*% d1[[i]]
                         if(iter >= 3) {
                           # first, nudge down to 1/2 step for everything
                           step <- 1/2 * step
                           # contract to no larger than the initial step
                           # slowly moving that max down as iteration number increases
                           ii <- 1
                           while(any(abs(step) > 3*u0/iter)) {
                             ii <- ii + 1
                             step <- 1/2 * step
                             if(ii > 20) {
                              stop("Ridiculous Newton step proposed, MAP not converging.")
                             }
                           }
                         } else {
                           # contract just a bit
                           step <- 1/2 * step
                         }
                         omg0[i,] - step
                       })
        
        # this recudes that to a vector
        omg0 <- t(do.call(cbind, omg0))
        # make omega
        omega0p <- c(tmpomega, list(omg0))
        while( length(omega0p) < length(omega0)) {
          omega0p[[length(omega0p)+1]] <- omega0[[length(omega0p)+1]] 
        }
        omega <- buildOmega(omega0=omega0p, groups=groups, nrowX=nrow(X))
        omegaFull <- buildOmega(omega0=omega0p, groups=groups, nrowX=nrow(X), full=TRUE)
        
      } # end while( max(abs( (omg1 - omg0) / pmax(abs(omg0),1E-5))) > 1E-3)
      if(verb) {
        cat("Estimates:\n")
        print(omg0)
      }
      #####################################################
      #         3) Estimate variance (Fishers I)          #
      #####################################################

      # add this to the tmpomega
      tmpomega <- c(tmpomega, list(omg0))
      tmpomegaFull <- omegaFull
      # get Qi
      drv <- d2
      # assumes there are exactly two levels
      Qi0[[oi]] <- lapply(1:length(drv), function(i) {
        ss <- scaleQuadPoints(drv[[i]], nzi)
        for(j in 1:nrow(ss)) {
          if(ss[j,j] > abs(omg0[i,j])) {
            ss[j,j] <- sqrt(ss[j,j]^2 + omg0[i,j]^2)
            omg0[i,j] <<- 0
          }
        }
        ss
      })
      map <- groups[,oi-1]
      umap <- unique(map)
      nzi <- ncol(Z[[oi]]) # number of Z columns at olvl
      Qi[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(weights[[oi-1]]))
      for(i in 1:nrow(weights[[oi-1]])) {
        Qi[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
      }
      QiFull[[oi]] <- matrix(0, nrow=nzi, ncol=nzi*nrow(X))
      for(i in 1:nrow(X)) {
        QiFull[[oi]][1:nzi,(i-1)*nzi+1:nzi] <- Qi0[[oi]][[(1:length(umap))[map[i]==umap] ]]
      }
      # add to weights for further MAP calculations only
      if(oi < length(omega0)) {
        # add detQ
        df <- data.frame(detQ=sapply(Qi0[[oi]],det))
        # add group labels
        groupNames <- colnames(groups)[oi-1]
        for(i in 1:length(groupNames)) {
          df[,groupNames[i]] <- unique(groups[,groupNames[i]])
        }
        # add detQ to weights for later iterations to use
        weights[[oi]] <- merge(weights[[oi]],df[,c(groupNames, "detQ")],by.x="index", by.y=groupNames)
      }
    } # ends for( oi in 2:length(omega0))
    
    list(omega0=tmpomega, omega=omega, omegaFull=tmpomegaFull, Qi0=Qi0, Qi=Qi, QiFull=QiFull)
  }
}


# helper functions for adapter:

# turn the genD output into the Cholesky of the inverse Fisher information
#' @importFrom Matrix nearPD
scaleQuadPoints <- function(d2, nz){
  solved <- solve(-1*d2)
  res <- NULL
  tryCatch(res <- chol(solved),
           error= function(e) {
             # if solved matrix is not positive definite (PD), use nearPD
             # to find the nearest positive definite matrix.
             # If that fails, because matrix is near negative semi-definite, use
             # the absolute value of the diagonal as an approximation
             tryCatch(solved <- nearPD(solved)$mat,
                      error=function(e){
                        solved <<- diag(abs(diag(solved)))
                      })
             res <<- chol(solved)
           })
  res
}

# turn omega0 (which has the same number of rows as there are units at each level
# into omega (which has the same number of rows as the X data) 
buildOmega <- function(omega0, groups, nrowX, full=FALSE) {
  omega <- list(NULL)
  oind <- 1
  for(o0i in 2:length(omega0)) {
    omega0i <- as.matrix(omega0[[o0i]])
    res <- matrix(0, nrow=nrowX, ncol=ncol(omega0i))
    noind <- ncol(omega0i) # number of columns to copy over
    map <- groups[,o0i-1]
    umap <- unique(map)
    for(i in 1:length(umap)) {
      # for each unique level
      for(oindi in 1:noind) {
        # copy over data by column (oindi) on omega0i
        # effectively duplicated each value of omega0  a number of times coresponding to how many obs there are in that group
        res[which(map==umap[i]),oindi] <- omega0i[i,oindi]
      }
    }
    if(o0i > 2 & !full) {
      res <- res[!duplicated(groups[,o0i-2]),]
    }
    omega <- c(omega, list(res))
    oind <- oind + noind
  } # closes main for loop for(o0i in 2:length(omega0))
  omega
}

# function that returns the likelihood, by group, evaulated at a point
adapterLnL <- function(y, X, levels, Z, ZFull, weights, k, qp,
                       covariance_constructor, omega, omegaFull, omega0, par0,
                       verbose, Qi, QiFull, olvl, acc, family) {
  # olvl is the level we are optimizing now
  function(par, long=FALSE) {
    # notice par here is a new displacement in the random effect location
    # not par0, which is a set of mixed model parameters
    yadj <- 0
    o0 <- omega0
    nzi <- 0 # number of Z columns at olvl
 
    for(i in 1:olvl) {
      if(!is.null(Z[[i]])) {
        ki <- ncol(Z[[i]])
        if(i == olvl) {
          nzi <- ki
        }
        if(ki >= 1) {
          # keep track of update in omega for calculation of prior prob
          # yadj is moved by omega, but also to account for an addition change 
          # controlled by par
          zAdjust <- apply(ZFull[[i]] * omegaFull[[i]],1,sum)
          if(olvl == i) {
            zAdjust <- zAdjust + ZFull[[i]] %*% par[1:ki]
            for(kii in 1:ki) {
              o0[[i]][,kii] <- o0[[i]][,kii] + par[kii]
            }
            par <- par[-(1:ki)]
          }
          # for olvl > 2, zAdjust has the wrong number of rows and need to be mapped
          # back to yp using the weights
          yadj <- yadj + zAdjust
        } # ends if  if(ki >= 1) 
      } # ends  if(!is.null(Z[[i]]))
    } #ends for(i in 1:olvl)
    # evaluate the likelihood
    beta <- par0[1:k]
    parC <- covariance_constructor(par0[-(1:k)])

    # Set curavtures Q for lower level to matrix of 0 
    # This is becasue the curvature is 0 at the maximum 
    # used for calculated likelihood by groups. 
    Qi_ <- matrix(0, nrow=nzi, ncol=nzi*nrow(weights[[olvl-1]])) # this needs to be the Qi for lower levels
    Qi__ <- c(Qi, list(Qi_))
    QiFull_ <- matrix(0, nrow=nzi, ncol=nzi*nrow(X)) # this needs to be the Qi for lower levels
    QiFull__ <- c(Qi, list(QiFull_))

    #use helper funciton to return the likelihood contribution of teach group 
   
    loglikelihoodByGroup <- calc.lin.lnl.quad(y=y, yhat=X %*% beta + yadj, level=olvl,
                                              Z, Qi=Qi__,
                                              omega=lapply(omega, function(omegai) {0*omegai}),
                                              W=weights, C=parC, qp, top=FALSE,
                                              atPoint=TRUE, verbose=verbose,
                                              acc=acc, ZFull=ZFull,
                                              omegaFull=omegaFull, QiFull=QiFull__,
                                              family=family)
    Cl <- parC[[olvl]]
    # apply multivatiate normal pdf to get posterior for each gorup 
    posteriorByGroup <- apply(o0[[olvl]], MARGIN=1, function(p) {
      mvnpdfC(as.matrix(p), rep(0, length = length(p)), varcovM=Cl%*%t(Cl), Log=TRUE)
                        })
    if(long) {
      return(list(res=loglikelihoodByGroup + posteriorByGroup,
                  loglikelihoodByGroup=loglikelihoodByGroup,
                  posteriorByGroup=posteriorByGroup))
    }
    #final return 
    loglikelihoodByGroup + posteriorByGroup
  } #ends  funciton to be return function(par, long=FALSE) 
} 
