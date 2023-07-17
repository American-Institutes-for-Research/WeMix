

# The main function which calculates the analytic solution to the linear mixed effects model. 
# @param weights level-1 weights
# @param y outcome measure. 
# @param X the X matrix.
# @param Zlist, a list of matrixes with the Z values for each level. 
# @param Zlevels, the level corresponding to each matrix in Zlist. 
# @param weights a list of unconditional weights for each model level. 
# @param weightsC a list of conditional weights for each model level. 
# @param groupID a matrix containing the group ids for each level in the model. 
# @param lmeVardf a dataframe containing the variances and covariance of the random effects, in the same format as returned from lme. 
#' @importFrom Matrix Diagonal Matrix bdiag diag
#' @importFrom methods as .hasSlot rbind2 
analyticSolve <- function(y, X, Zlist, Zlevels, weights, weightsC=weights, groupID, lmeVarDF, analyticJacobian=FALSE, forcebaseQR=FALSE, qr_=qr_0, v0, lndetLz0=NULL) {

  if(! inherits(groupID, c("matrix", "data.frame")) ) {
    stop(paste0("Variable", dQuote("groupID")," must be a matrix  or data.frame with IDs at a level per column."))
  }
  # this has the indexes, but we do not need those. Save them here
  groupID0 <- groupID
  # now winnow down to just the expected number of columns, other cols added for conditional weighting
  groupID <- as.matrix(groupID[ , 1:(length(weights)-1), drop=FALSE])
  if(!inherits(y, "numeric")) {
    stop(paste0(dQuote("y"), " must be a numeric vector."))
  }
  ny <- length(y)
  X <- as.matrix(X)
  nX <- nrow(X)
  if(nX != ny) {
    stop(paste0("Length of the ", dQuote("y"), " vector and the ", dQuote("X"), " matrix must agree."))
  }
  if(length(Zlist) != length(Zlevels)) {
    stop(paste0("Length of the ", dQuote("Zlist"), " and ", dQuote("Zlevels"), " must agree."))
  }
  nW1 <- length(weights[[1]])
  if(nW1 != nX) {
    stop(paste0("Number of rows in ", dQuote("weights[[1]]"), " must agree with number of rows in ", dQuote("X"), "."))
  }
  nWc1 <- length(weightsC[[1]])
  if(nWc1 != nX) {
    stop(paste0("Number of rows in ", dQuote("weightsC[[1]]"), " must agree with number of rows in ", dQuote("X"), "."))
  }
  ngid <- nrow(groupID)
  if(ngid != nX) {
    stop(paste0("Number of rows in ", dQuote("groupID"), " must agree with number of rows in ", dQuote("X"), "."))
  }
  # groupID needs to be one indexed on the top level, so enforce that (this is used for robustSE)
  for(gLevel in 1:ncol(groupID)) {
    groupID[ , gLevel] <- as.numeric(as.factor(groupID[ , gLevel]))
  }
  # get the number of groups per level
  nc <- apply(groupID, 2, function(x) { length(unique(x)) } )
  
  # for level 1 cast Z for lowest level as a matrix and transpose to get Zt
  Zt <- Matrix::t(as(Zlist[[1]], "Matrix"))
  
  # for level, build PsiVec, the diagonal of the Psi  (weights) matrix
  PsiVec <- rep(weights[[Zlevels[1]]], ncol(Zlist[[1]])/nc[1])
  
  # Assemble the Zt matrix and psi Vector for levels >1 
  ZColLevels <- rep(Zlevels[1], ncol(Zlist[[1]]))
  if(length(Zlist) > 1) {
    for(zi in 2:length(Zlist)) {
      Zt <- rbind(Zt, Matrix::t(as(Zlist[[zi]], "Matrix")))
      ZColLevels <- c(ZColLevels, rep(Zlevels[zi], ncol(Zlist[[zi]])))
      PsiVec <- c(PsiVec, rep(weights[[Zlevels[zi]]], ncol(Zlist[[zi]])/nc[Zlevels[zi]-1]))
    }
  }
 
  # get unit weights
  W0 <- weights[[1]]
  Psi <- Diagonal(n=nrow(Zt), x=PsiVec) #matrix of weights at level one
  Psi12 <- Diagonal(n=nrow(Zt), x=sqrt(PsiVec)) #matrix of square root of level one weights
  W <- Diagonal(x=(W0))
  W12 <- Diagonal(x=sqrt(W0))
  
  # calculate  conditional weights matrix where level one weights are scaled by level two wieghts
  W12cDiag <- W0
  L1groups <- unique(groupID[ , 1])
  for(gi in 1:length(L1groups)) {
    rowsGi <- groupID[,1] == L1groups[gi]
    W12cDiag[rowsGi] <- sqrt(W12cDiag[rowsGi] / weights[[2]][gi])
  }
  W12c <- Diagonal(x=W12cDiag)
  Z <- Matrix::t(Zt)
  options(Matrix.quiet.qr.R=TRUE) #supress extra print outs 

  M0 <- Matrix(data=0, ncol=ncol(X), nrow=nrow(Zt))

  # to build the Z we need iDelta, which requires a theta.
  # but theta changes by evaluation. So build a pseudo-Delta (not actual theta)
  # only to allow Z to be cached. iDelta is needed only for dimensions,
  # no values are used.
  iDelta <- Delta <- list()
  levels <- max(lmeVarDF$level)
  for (l in 2:levels) {
    #set up matrix of zeros with rows and columns for each random effect
    #number of random effect is the number of variance terms at that level (not including covariance terms )
    n_ref_lev  <- nrow(lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),] ) 
    lambda_i  <- matrix(0, nrow=n_ref_lev, ncol=n_ref_lev)
    row.names(lambda_i) <- lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2), "var1"]
    colnames(lambda_i) <- lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2), "var1"]
    #get only v for this level
    group <- lmeVarDF[lmeVarDF$level==l, "grp"][1]
    v_lev <- v0[grep(paste0("^", group, "."), names(v0))]
    #fill in lambda_i from theta using names
    for (vi in 1:length(v_lev)){
      row_index <- strsplit(names(v_lev[vi]), ".", fixed=TRUE)[[1]][2]
      col_index <- ifelse(length(strsplit(names(v_lev[vi]), ".", fixed=TRUE)[[1]]) > 2, strsplit(names(v_lev[vi]), ".", fixed=TRUE)[[1]][3], row_index)
      lambda_i[row_index, col_index] <- v_lev[vi]
    }
    fixed <- solveFix(lambda_i) # correct non-invertable matrixes to nearby matrixes
    iDelta[[l]] <- fixed$lambda
    Delta[[l]] <- fixed$lambdaInv
  } # end for (l in 2:levels)
  ZiAl <- list()
  l0 <- list()
  for(i in 1:ncol(groupID)) {
    l0 <- c(l0, list(NULL))
  }
  uvl <- lapply(unique(groupID[ , ncol(groupID)]), function(tid) { l0 } )
  uv0 <- 1:ncol(Z)
  # build Z to cache it
  for(level in 1:ncol(groupID)) {
    groups <- unique(groupID[,level])
    Deltai <- Delta[[level+1]]
    topLevel <- level == ncol(groupID)
    Zl <- Z[ , ZColLevels==(level+1), drop=FALSE] # level l Z
    uv0l <- uv0[ZColLevels==(level+1)] # base of possible RE for every group at this level
    goalN <- length(uv0l)/length(groups)
    for(gi in 1:length(groups)) {
      # get Z rows
      Zrows <- groupID[ , level] == groups[gi]
      if(level == 1 || level < ncol(groupID)) {
        if(!topLevel) {
          lp1 <- unique(groupID[Zrows, level+1])
          if(length(lp1) > 1) {
            stop("Not a nested model; WeMix only fits nested models.")
          }
        }
      }
      # filter to the rows for this group, also filter to just columns
      # for this level of the model
      # Zi, used for forming ZiA (stored in ZiAl)
      # Zil, used for finding pointers (stored in uvl)
      if(level == 1 || level < ncol(groupID)) {
        Zi <- Z[Zrows, , drop=FALSE]
        Zil <- Zl[Zrows, , drop=FALSE] 
      } else {
        Zil <- Zi <- Zl[Zrows, , drop=FALSE] # Zi associated with the random effect (~ Z_g in the specs)
      }
      # within columns for this level, only the non-zero columns
      # regard this unit, so filter it thusly
      # the below code uses the pointers from the sparse matrix to identify non zero rows
      # it is equivalent to
      # Zcols <- apply(Zi,2,nozero)
      
      if(.hasSlot(Zi,"p")){
        z_pointers <- Zi@p
        len_pointers <- length(z_pointers)
        Zcols <- which(z_pointers[1:(len_pointers-1)]  - z_pointers[2:len_pointers] !=  0) 
        z_pointers <- Zil@p
        len_pointers <- length(z_pointers)
        Zcolsl <- which(z_pointers[1:(len_pointers-1)]  - z_pointers[2:len_pointers] !=  0) 
      } else {
        Zcols <- apply(Zi, 2, nozero)
        Zcolsl <- apply(Zil, 2, nozero)
      }
      if(length(Zcolsl) < goalN) {
        # if there is a Z column that was identified
        if(length(Zcolsl) > 0) {
          # grab the name of this group
          colname <- colnames(Zil)[Zcolsl[1]]
          # and select all columns with that name
          Zcols <- (1:ncol(Z))[colnames(Z) == colname]
          Zcolsl <- (1:ncol(Z))[colnames(Z) == colname & ZColLevels==(level+1)]
          if(length(Zcols) != goalN) {
            stop("Unable to construct Z matrix.")
          }
        } else {
          # Z is entirely 0, so we can just take the first goalN columns
          # because we know they are all 0
          Zi <- Zil <- Zi[ , 1:goalN, drop=FALSE]
          Zcols <- Zcolsl <- NULL
        }
      } # end if(length(Zcolsl) < goalN)
      if(!is.null(Zcols)) {
        # subset Zi to just columns for this level
        Zi <- Zi[ , Zcols, drop=FALSE]
      }
      # uvi is the column indexes for Zi for this group, at this level
      uvi <- uv0l[Zcolsl]
      
      if(level == 1 || level < ncol(groupID)) {
        ZiA <- W12cDiag[Zrows] * Zi
        ZiAl <- c(ZiAl, list(ZiA))
      }
      tgi <- groupID[groupID[ , level] == groups[gi], ncol(groupID)][1]
      # store the Z columns for REs at this level for this group in uvl
      # for non-top level, this should be many REs.
      uvl[[tgi]][[level]] <- c(uvl[[tgi]][[level]], list(uvi))
    } # end for(gi in 1:length(groups))
  } # end for(level in 1:ncol(groupID))

  preKR <- lapply(2:levels, function(l) {
    diag(1, unique(lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),"ngrp"]))
  })

  # v is the omega vector
  function(v, verbose=FALSE, beta=NULL, sigma=NULL, robustSE=FALSE, returnBetaf=FALSE, getGrad=FALSE, lndetLz0u=TRUE) {
    beta0 <- beta
    sigma0 <- sigma
    
    #Create list delta and lambda Matrixes for each level
    # iDelta is the inverse Delta (lambda pre kronecker) for forming VC estimates in the end
    iDelta <- Delta <- list()
    lambda_by_level <- list()
    
    levels <- max(lmeVarDF$level)
    
    LambdaL <- list()

    for (l in 2:levels){
      #set up matrix of zeros with rows and columns for each random effect
      #number of random effect is the number of variance terms at that level (not including covariance terms )
      n_ref_lev <- nrow(lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2),] ) 
      lambda_i <- matrix(0, nrow=n_ref_lev, ncol=n_ref_lev)
      row.names(lambda_i) <- lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2), "var1"]
      colnames(lambda_i) <- lmeVarDF[lmeVarDF$level==l & is.na(lmeVarDF$var2), "var1"]
      
      #get only v for this level
      group <- lmeVarDF[lmeVarDF$level==l, "grp"][1]
      v_lev <- v[grep(paste0("^",group,"."), names(v))]
      
      #fill in lambda_i from theta using names
      for (vi in 1:length(v_lev)){
        row_index <- strsplit(names(v_lev[vi]), ".", fixed=TRUE)[[1]][2]
        col_index <- ifelse(length(strsplit(names(v_lev[vi]), ".", fixed=TRUE)[[1]])>2, strsplit(names(v_lev[vi]), ".", fixed=TRUE)[[1]][3], row_index)
        lambda_i[row_index,col_index] <- v_lev[vi]
      }
      fixed <- solveFix(lambda_i) # correct non-invertable matrixes to nearby matrixes
      iDelta[[l]] <- fixed$lambda
      Delta[[l]] <- fixed$lambdaInv
      #assemble blockwise lambda by level matrix 
      lambda_by_level[[l]] <- kronecker(lambda_i, preKR[[l-1]])
      LambdaL[[l-1]] <- fixed$lambda
    }

    #Create lambda matrix from diagonal of all lamdas 
    lambda <- bdiag(lambda_by_level[2:length(lambda_by_level)]) # vignette equation 21

    yaugmented <- c(as.numeric(W12 %*% y), rep(0, nrow(Zt))) # first matrix in vingette equation 25 and 60 
    
    A <- rbind(cbind(W12 %*% Z %*% lambda, W12 %*% X),
               cbind(Psi12, M0)) # second matrix in equation (60)/(25). 
    qr_a <- qr_(A, allowPermute=TRUE)
    # equation (26)/ stack vector (u: top part for fixed effects coef estimate, b: bottom part for random effects vector (i.e. of length n_group * n_random_effects)
    suppressWarnings(ub <- Matrix::qr.coef(qr_a, yaugmented))

    # get the gradients
    if(getGrad) {
      b <- ub[-(1:ncol(Z))]
      u <- ub[1:ncol(Z)]
      R <- chr_(A, qr_a)
      Rrows <- nrow(R) - length(b):1+1
      Rrows <- Rrows[Rrows > length(u)]
      RX <- R[Rrows, ncol(R) - length(b):1 + 1, drop=FALSE]
      db <- beta0 - b
      db <- db[colnames(RX)]
      db[is.na(db)] <- beta0[is.na(db)]
      res0 <- res <- 0 * db
      for(bi in 1:length(b)) {
        ei <- 0 * res
        ei[bi] <- 1
        #(RX*RX) is elementwise multiplication
        res[bi] <- -1*(db[bi]) * sum((RX*RX) %*% ei)/sigma^2
      }
      return(res)
    }

    # get ln of determinant of Lz matrix (Bates et al. 2015 notation) / calculation of alpha
    # or the final product in Pinheiro and Bates 2.13 (formula used for weights)
    lndetLz <- 0 # lndetLz is alpha in the specs
    lndetLzg <- list() # by top level group lnldetLz
    # used to find columns with no-non-zero elements
    # grab and the random effect vector
    u <- ub[1:ncol(Z)]
    names(u) <- colnames(Z)
    IMatCols <- ncol(ZiAl[[1]])
    groupsTop <- unique(groupID[ , ncol(groupID)])
    if(lndetLz0u) {
      # only evaluate all of this if we do not already have an lnldetZ value
      for(level in 1:ncol(groupID)) { # really level -1
        groups <- unique(groupID[ , level])
        Deltai <- Delta[[level + 1]]
        # R11 notation from Bates and Pinheiro, 1998
        R11 <- list()
        topLevel <- level == ncol(groupID)
        qi <- nrow(Deltai)
        if(level > 1) {
          # lambdaLev already used at level 1
          # but still postpend new IMat
          IMat <- matrix(0, nrow=qi, ncol=IMatCols)
          IMatCols <- IMatCols - qi
        } else {
          LambdaLev <- bdiag(LambdaL)
          IMat <- matrix(0, nrow=qi, ncol=IMatCols)
          IMatCols <- IMatCols - qi
        } 
        diag(IMat) <- 1
        for(gi in 1:length(groups)) {
          if(level == 1) {
            # get Zi
            ZiA <- ZiAl[[gi]]
            Zrows <- groupID[ , level] == groups[gi]
            ltop <- unique(groupID[Zrows, ncol(groupID)])
            if(!topLevel) {
              # unit ID for level above this
              lp1 <- unique(groupID[Zrows, level+1])
              qp1 <- nrow(Delta[[level+2]])
            }
            if(ncol(ZiA) < qi) {
              ZiA <- cbind(ZiA, matrix(0, nrow=nrow(ZiA), ncol=qi - ncol(ZiA)))
              stop("malformed ZiA")
            }
            ZiA <- ZiA %*% LambdaLev
            # rbind manually, used to do this:
            ZiA <- rbind2(ZiA, IMat) # rbind(Z lambda, sqrt(Psi))
            # in this section we decompose Zi * A using qr decomposition,
            # following equation 43 in the vingette.
            ZiAR <- qr_qrr(ZiA) # this is faster than getchr(ZiA)
            R22 <- ZiAR[1:qi, 1:qi, drop=FALSE]
            lndetLzi <- weights[[level+1]][gi] * sum(log(abs(Matrix::diag(R22))))
            # save R11 for next level up, if there is a next level up
            if (!topLevel) {
              R11i <- sqrt(weightsC[[level+1]][gi])*ZiAR[-(1:qi), -(1:qi), drop=FALSE]
              # load in R11 (a list defined above) to the correct level
              # this if/else allows lp1 to be out of order
              if(length(R11) < lp1 || is.null(R11[[lp1]])) {
                # there was no R11, so start it off with this R11
                R11[[lp1]] <- R11i
              } else {
                R11[[lp1]] <- rbind(R11[[lp1]], R11i)
              }
              if(length(lndetLzg) < ltop || is.null(lndetLzg[[ltop]])) {
                lndetLzg[[ltop]] <- lndetLzi
              } else {
                lndetLzg[[ltop]] <- lndetLzg[[ltop]] + lndetLzi
              }
            } else {
              if(length(lndetLzg) < ltop || is.null(lndetLzg[[ltop]])) {
                lndetLzg[[ltop]] <- lndetLzi
              } else {
                lndetLzg[[ltop]] <- lndetLzg[[ltop]] + lndetLzi
              }
            }
          } else { # end if(level == 1)
            # level >= 2
            ZiA <- rbind(pR11[[groups[gi]]], IMat)
            R <- qr_qrr(ZiA) # probably faster than getchr(ZiA)
            # weight the results
            # groups goes back to this unit
            if (!topLevel) {
              Zrows <- groupID[ , level] == groups[gi]
              # unit ID for level above this
              lp1 <- unique(groupID[Zrows, level+1])
              ltop <- unique(groupID[Zrows, ncol(groupID)])
              R11i <- sqrt(weightsC[[level+1]][gi])*R[-(1:qi), -(1:qi), drop=FALSE]
              # load in R11 (a list defined above) to the correct level
              # this if/else allows lp1 to be out of order
              if(length(R11) < lp1 || is.null(R11[[lp1]])) {
                # there was no R11, so start it off with this R11
                R11[[lp1]] <- R11i
              } else {
                R11[[lp1]] <- rbind(R11[[lp1]], R11i)
              }
              lndetLzi <- weights[[level+1]][groups[gi]] * sum(log(abs(Matrix::diag(R)[1:qi])))
              lndetLzg[[ltop]] <- lndetLzg[[ltop]] + lndetLzi
            } else {
              # top level
              lndetLzi <- weights[[level+1]][groups[gi]] * sum(log(abs(Matrix::diag(R)[1:qi])))
              lndetLzg[[groups[gi]]] <- lndetLzg[[groups[gi]]] + lndetLzi
            }
          } # end else for if(level == 1)
          lndetLz <- lndetLz + lndetLzi
        } # end for(gi in 1:length(groups))
        if(level < ncol(groupID)) {
          pR11 <- R11
        } 
      } #end for(level in 1:ncol(groupID))
    } # end if(lndetLz0)
    # this is the beta vector
    b <- ub[-(1:ncol(Z))]
    if(!is.null(beta0)) {
      if(length(b) != length(beta0)) {
        stop(paste0("The argument ", dQuote("beta"), " must be a vector of length ", length(b), "."))
      }
    }
  
    # wrap the random effect vector to matrix format so that each row regards one group
    bb <- list()
    vc <- matrix(nrow=1+length(v), ncol=1+length(v))
    u0 <- 0
    for(li in 2:length(lambda_by_level)) {
      lli <- lambda_by_level[[li]]
      ni <- nrow(lli)
      bb[[li]] <- lli %*% u[u0 + 1:ni]
      rownames(bb[[li]]) <- names(u[u0 + 1:ni])
      u0 <- u0 + ni
      vc[li,li] <- 1/(ni) * sum((bb[[li]])^2) # mean already 0
    }
    # the discrepancy ||W12(y-Xb-Zu)||^2 + ||Psi(u)||^2
    discrep <- discf(y, Zt, X, lambda, u, Psi12, W12, b)
    # the R22 matrix, bottom right of the big R, conforms with b
    R <- chr_(A, qr_a)
    Rrows <- nrow(R) - length(b):1 +1
    Rrows <- Rrows[Rrows > length(u)]
    RX <- R[Rrows, ncol(R) - length(b):1 +1, drop=FALSE]
    nx <- sum(W0)
    # residual
    sigma <- ifelse(is.null(sigma0), sqrt(discrep / nx), sigma0)
    dev <- 0
    if(lndetLz0u) {
      # lndetLz may be a Matrix, convert it
      dev <- 2*as.numeric(lndetLz)
    } else{
      dev <- 2*lndetLz0
      lndetLz <- lndetLz0
    }
    dev <- dev + nx*log(2*pi*(sigma^2)) + discrep/sigma^2
    # fix rank-deficent RX by setting values to 0
    Rprob <- !rkfnS(R) # do the diagonals add to rank
    # to be on the diagonal of RX, it needs to both be a row and a column
    RXprob <- Rprob[intersect(Rrows, ncol(R)-length(b):1+1)] # for RX
    if(any(RXprob)) {
      Matrix::diag(RX)[RXprob] <- 0
    }
    # add R22 term if beta is not beta-hat
    if(returnBetaf) {
      # a function that returns the lnl varying beta only
      f <- function(beta0) {
        db <- beta0 - b
        # for an illconditioned (or uninvertable) A-matrix there will be NAs
        # if there is an NA, then use the beta value--that is, assume betahat is 0
        db[is.na(db)] <- beta0[is.na(db)]
        # RX will be reordered if X has all 0 columns, db needs to be rearanged similarly
        db <- db[colnames(RX)]
        # Matrix::qr.R makes bad colnames, fix that
        db[is.na(db)] <- 0
        # deviance
        dev <- dev + sum( (RX %*% db)^2 )/sigma^2
        res <- return(dev/-2)
      }
      return(f)
    }
    # add in term for beta not beta-hat
    db <- 0 * b
    dbTerm <- 0
    if(!is.null(beta0)) {
      # difference in beta and betahat (b)
      db <- beta0 - b
      # if there is an NA, then use the beta value--that is, assume betahat is 0
      db[is.na(db)] <- beta0[is.na(db)]
      # RX will be reordered if X has all 0 columns, db needs to be rearanged similarly
      db <- db[colnames(RX)]
      # deviance
      dbTerm <- sum( (RX %*% db)^2 )/sigma^2
      dev <- dev + dbTerm
    }
    #Variance of Beta, from Bates et.al. 2015
    # qr matrix rank, tolerance value and method based on Matrix::rankMatrix API
    if(rkfn(R) == ncol(A)) {
      RXi <- Matrix::solve(RX)
      cov_mat <- (sigma^2) * RXi %*% Matrix::t(RXi)
      varBeta <- Matrix::diag(cov_mat)
    } else {
      cov_mat <- NULL
      varBeta <- rep(NA, length(b))
    }

    sigma <- ifelse(is.null(sigma0), sqrt(discrep / nx), sigma0)
    res <- list(b=b, u=u, ranef=bb, discrep=discrep, sigma=sigma, dev=dev,
                lnl=dev/-2, theta=v, varBeta=varBeta, vars=NULL, cov_mat=cov_mat,
                lndetLz=lndetLz, iDelta=iDelta, db=db, dbTerm=dbTerm, nx=nx, RX=RX, R=R)
    if(robustSE) {
      #Calculate robust standardize effect
      # based on the sandwich estimator of Rabe-Hesketh and Skrondal 2006
      # this whole calculation focuses on calculating the liklihood at the top group level 
      fgroupID <- groupID[ , ncol(groupID)] 
      uf <- unique(fgroupID) # unique final groupIDs
      # store lnl to check later
      lnli2 <- lnli <- vector(length=length(uf))
      Jacobian <- matrix(NA, nrow=length(uf), ncol=length(b)+length(v)+1)
      bwiL <- list()
      #this code block seperates wieghts into the set of weights belonging to each top level group 
      # components of discrep
      wres <- W12 %*% (y - Matrix::t(Zt) %*% lambda %*% u - X %*% b) # residual
      ures <- Psi12 %*% u # augmented ehat
      for(gi in 1:length(uf)) {
        giuf <- uf[gi] # index for this unit
        sgi <- fgroupID == giuf # subset for group gi
        weightsPrime <- list(weights[[1]][sgi])
        weightsCPrime <- list(weightsC[[1]][sgi])
        for(i in 2:length(weights)) {
          theseGroups <- unique(groupID[sgi, i-1])
          weightsPrime[[i]] <- weights[[i]][theseGroups]
          weightsCPrime[[i]] <- weightsC[[i]][theseGroups]
        }
        condenseZ <- function(z, uvi) {
          res <- z[sgi, uvi, drop=FALSE]
        }
        groupIDi <- groupID[sgi,,drop=FALSE]
        lmeVarDFi <- lmeVarDF
        lmeVarDFi$ngrp[lmeVarDFi$level==1] <- sum(sgi)
        for(i in 2:length(weights)) {
          lmeVarDFi$ngrp[lmeVarDFi$level==i] <- length(unique(groupIDi[,i-1]))
        }
        #calculate the group level likelihood by applying the function to only the x or y within the group indexed by sgi
        # overall dev (lnl=dev/-2) 
        # dev <- 2*lndetLz + nx*log(2*pi*(sigma^2)) + discrep/sigma^2
        # VERIFIED: uvl is mapped to the index ID, so use giuf
        lnli2[gi] <- 2*lndetLzg[[giuf]]/-2 + sum(sum(W0[sgi]))*log(2*pi*sigma^2)/-2 + sum(wres[sgi]^2)/sigma^2/-2 + sum(ures[unlist(uvl[[giuf]])]^2)/sigma^2/-2
        # zero index for this level, used in Zlist argument a few lines down
        ind0 <- c(0,cumsum(lapply(Zlist, ncol)))
        bwi <- analyticSolve(y=y[sgi], X=X[sgi,,drop=FALSE],
                  Zlist=lapply(1:length(Zlist), function(lvl) {
                      cols <- sort(unlist(uvl[[giuf]][[lvl]])) - ind0[lvl]
                      tryCatch(Zlist[[lvl]][sgi, cols, drop=FALSE], error=function(e){stop("Could not build Zlist.")})
                    }),
                  Zlevels=Zlevels,
                  weights=weightsPrime,
                  weightsC=weightsCPrime,
                  groupID=groupIDi,
                  lmeVarDF=lmeVarDFi,
                  v0=v,
                  lndetLz0=lndetLzg[[giuf]])
        tryCatch(lnli[gi] <- bwi(v=v, verbose=verbose, beta=b, sigma=sigma, robustSE=FALSE)$lnl,
                 error= function(e) {
                   lnli[gi] <<- NA
                 })
        if( abs(lnli2[gi] - lnli[gi]) > 0.001) {
          # sometimes Matrix::qr tries to solve a singular system and fails
          # normally base::qr works in these cases, so use that instead
          bwi <- analyticSolve(y=y[sgi], X[sgi,,drop=FALSE],
                               Zlist=lapply(1:length(Zlist), function(lvl) {
                                       cols <- sort(unlist(uvl[[giuf]][[lvl]])) - ind0[lvl]
                                       tryCatch(Zlist[[lvl]][sgi, cols, drop=FALSE], error=function(e){stop("Could not build Zlist.")})
                                      }),
                               Zlevels=Zlevels,
                               weights=weightsPrime,
                               weightsC=weightsCPrime,
                               groupID=groupIDi,
                               lmeVarDF=lmeVarDFi,
                               qr_=qr_s,
                               v0=v,
                               lndetLz0=lndetLzg[[giuf]]) 
           bwiL <- c(bwiL, list(bwi))
           tryCatch(lnli[gi] <- bwi(v=v, verbose=verbose, beta=b, sigma=sigma, robustSE=FALSE)$lnl,
                    error= function(e) {
                      lnli[gi] <<- NA
                    })
        } else {
           bwiL <- c(bwiL, list(bwi))
        }

        # if the function can be evaluated for this group.
        if(!is.na(lnli[gi])) {
          # function for partials with respect to beta
          bwiW <- function(bwi, v, sigma, b, ind) {
            bwip <- bwi(v=v, verbose=FALSE,
                        b=b, sigma=sigma,
                        robustSE=FALSE, returnBetaf=TRUE, lndetLz0u=FALSE)
            function(ind) {
              function(par) {
                b[ind] <- par
                bwip(b)
              }
            }
          }
          # function for partials with respect to theta
          bwiT <- function(bwi, v, sigma, t, ind) {
            function(par) {
              if(ind > length(v)) {
                sigma <- par
              } else {
                v[ind] <- par
              }
              bwi(v=v, verbose=FALSE,
                  b=b, sigma=sigma, robustSE=FALSE, lndetLz0u=TRUE)$lnl
            }
          }
          if(analyticJacobian) {
            Jacobian[gi,1:length(b)] <- bwi(v=v, sigma=sigma, b=b, robustSE=FALSE, getGrad=TRUE)
          } else {
            bw <- bwiW(bwi, v=v, sigma=sigma, b=b)
            for(j in 1:length(b)){
              Jacobian[gi,j] <- d(bw(j), par=b[j])
            }
          }
          for(j in 1:(1+length(v))) {
            v0 <- c(v, sigma)[j]
            Jacobian[gi, length(b) + j] <- d(bwiT(bwi, v=v, sigma=sigma, b=b, ind=j), par=v0)
          }
        } else { # end if(!is.na(lnli[gi]))
          warning("Some top level units variance component cannot be computed. Try combining top level units.")
        }
      }
      lnl1 <- dev/-2
      lnliT <- sum(lnli)
      if(!is.na(lnliT)) {
        if( abs(lnl1 - lnliT) > 0.5) {
          # this would be a large difference
          warning("Likelihood estimated at the top group level and summed disagrees with overall likelihood. Standard errrors may not be accurate.")
        }
      }
      J <- matrix(0,ncol=ncol(Jacobian), nrow=ncol(Jacobian))
      nr <- 0
      # a cross product that removes rows (groups) where lnl cannot be evaluated
      for (i in 1:nrow(Jacobian)){
        if(!is.na(lnli[i])) {
          J  <- J  + Jacobian[i,]  %*% t(Jacobian[i,])
          nr <- nr + 1
        }
      }
      J <- (nr/(nr-1))*J
      # just beta part of J
      varBetaRobust <- as(cov_mat %*% J[1:length(b), 1:length(b)] %*% cov_mat , "matrix")
      colnames(J) <- rownames(J) <- c(names(b), names(v), "sigma")
      resid <- y - Matrix::t(Zt) %*% lambda %*% u - X %*% b
      res <- c(res, list(varBetaRobust=varBetaRobust,
                         seBetaRobust=sqrt(diag(varBetaRobust)), iDelta=iDelta,
                         Jacobian=J,
                         resid=resid))
    }
    return(res)
  }
}

# helpers just for analyticSolve

#discrepency function for calculating least squares solution. Follows from Bates et al. 2015 eq 14 and 15
discf <- function(y, Zt, X, lambda, u, Psi12, W12, b) { # return ehatT * ehat
  if(any(is.na(b))) {
    X <- X[,!is.na(b),drop=FALSE]
    b <- b[!is.na(b)]
  }
  wres <- W12 %*% (y - Matrix::t(Zt) %*% lambda %*% u - X %*% b) # residual
  ures <- Psi12 %*% u # augmented ehat
  as.numeric(sum(wres^2) + sum(ures^2))
}

# the sparse QR permultes rows, which causes problems for the method of solving
# used in WeMix. qr_0 is robust to that, and wehn allowPermulte=TRUE it will
# sense that Matrix::qr did that and instead use base::qr.
qr_0 <- function(X, allowPermute=FALSE) { 
  # Matrix::qr does not do well in this case, use base::qr
  if(nrow(X) < ncol(X)) {
    # here base QR elegantly removes coefficients and returns NAs in their place
    X <- as(X, "matrix")
    return(base::qr(X))
  }
  # try to use Matrix
  qr1 <- Matrix::qr(X)
  if(allowPermute || inherits(qr1,"qr")) {
    # if allowPermute (allow Matrix::qr to use permutations) or this is a base::qr, just return
    return(qr1)
  }
  # do not allow permutations. Test if they happened. If they did, use base::qr instead
  if( any(range( order(qr1@q) - 0:max(qr1@q)) > 0) ) {
    X <- as(X, "matrix")
    return(base::qr(X))
  }
  # base case, return the non-permuted Matrix::qr result
  return(qr1)
}

# just use base:qr when qr_s is selected
qr_s <- function(X, allowPermute=FALSE) { 
  X <- as(X, "matrix")
  return(base::qr(X))
}

# run QR and return R
qr_qrr <- function(X) { 
  if(nrow(X) < ncol(X)) {
    X <- as(X, "matrix")
    return(base::qr.R(base::qr(X)))
  }
  qr1 <- Matrix::qr(X)
  if(inherits(qr1,"qr")) {
    return(base::qr.R(qr1))
  }
  # check for permutations and use base if there are some
  if( any(range( order(qr1@q) - 0:max(qr1@q)) > 0) ) {
    return(base::qr.R(base::qr(X)))
  }
  # we can use Matrix:qr
  return(Matrix::qrR(qr1, backPermute=FALSE))
}

# run a chol of AtA and return R
getchr <- function(A) {
  #sometimes we pass an illconditioned matrix and don't care, so eat the warnings
  suppressWarnings(m0 <- as(Matrix::chol(Matrix::t(A) %*% A, pivot=FALSE), "Matrix"))
  colnames(m0) <- colnames(A) # allowable because pivot=FALSE
  return(m0)
}

# get R robust to uninvertable matrixes, falling back on the qr when needed
chr_ <- function(A, qr) {
  m <- tryCatch(getchr(A),
                error=function(e) {
                  if(inherits(qr,"qr")) {
                    return(base::qr.R(qr))
                  } else {
                    A <- as(A, "matrix")
                    return(base::qr.R(base::qr(A)))
                  }
                })
  return(m)
}

# rank of a qr
rkfn <- function(qr_R) {
  sum(rkfnS(qr_R))
}

# is this element of the diag > the tolerance (does it add a value to rank)
rkfnS <- function(qr_R) {
  d <- Matrix::diag(qr_R)
  # add 1 to tolerance to account for imprecision in .Machine$double.eps
  tol <- (length(d)+1) * .Machine$double.eps
  proposed <- abs(d) > max(abs(d)) * tol
  # the condition is a bit different for a "matrix" than a "Matrix"
  if(is.matrix(qr_R)) {
    # increase until condition number of remaining matrix agrees with rank from base package
    while(sum(proposed) > 0 & rcond(qr_R[proposed,proposed]) <= .Machine$double.eps) {
      proposed[proposed][which.min(abs(d[proposed]))] <- FALSE
    }
  }
  return(proposed)
}

# a wrapper for analyticSolve that allows it to be called from an optimizer. Takes the same arguments as analyticSolve. 
# @param weights level-1 weights
# @param y outcome measure. 
# @param X the X matrix.
# @param Zlist, a list of matrixes with the Z values for each level. 
# @param Zlevels, the level corresponding to each matrix in Zlist. 
# @param weights a list of unconditional weights for each model level. 
# @param weightsC a list of conditional weights for each model level. 
# @param  groupID a matrix containing the group ids for each level in the model. 
# @param  lmeVardf a dataframe containing the variances and covariance of the random effects, in the same format as returned from lme. 
#' @importFrom Matrix Diagonal Matrix
#' @importFrom methods as
devG <- function(y, X, Zlist, Zlevels, weights, weightsC=weights, groupID, lmeVarDF, v0) {
  bs <- analyticSolve(y=y, X=X, Zlist=Zlist, Zlevels=Zlevels, weights=weights, weightsC=weightsC, lmeVarDF=lmeVarDF,
           groupID=groupID, v0=v0)
  function(v, getBS=FALSE) {
    if(getBS) {
      return(bs)
    }
    bhat <- bs(v)
    return(bhat$dev)
  }
}

solveFix <- function(lambda_i) {
  # try QR based fix, which may fail
  res <- try(qrFix(lambda_i), silent=TRUE)
  if(inherits(res, "try-error")) {
    res <- try(eigenFix(lambda_i), silent=TRUE)
    if(inherits(res, "try-error")) {
      return(svdFix(lambda_i))
    }
    return(res)
  } else{
    return(res)
  }
  # then try eigenvalue based fix
}

# return a matrix and its inverse, being robust to 0s
qrFix <- function(lambda_i) {
  #set any 0s to smallest non zero number to enable solve in next step
  diag(lambda_i)[abs(diag(lambda_i)) < .Machine$double.eps] <- .Machine$double.eps
  return(list(lambda=lambda_i, lambdaInv=solve(lambda_i)))
}

# return a matrix and its inverse, being robust to 0s
eigenFix <- function(lambda_i) {
  # qrfix was singular, so use eigen
  eig <- eigen(lambda_i, symmetric=FALSE)
  eigV <- eig$values
  eigV[abs(eigV) < 2 * .Machine$double.eps * max(abs(eigV))] <- 2 * .Machine$double.eps * max(abs(eigV))
  eigLambdai <- diag(1/eigV, nrow=length(eigV))
  eigQ <- eig$vectors
  eigQinv <- solve(eigQ)
  eigInv <- eigQ %*% eigLambdai %*% eigQinv
  eig <- eigQ %*% diag(eigV, nrow=length(eigV)) %*% eigQinv
  dimnames(eig) <- dimnames(eigInv) <- dimnames(lambda_i)
  return(list(lambda=eig, lambdaInv=eigInv))
}

svdFix <- function(lambda_i) {
  svdi <- svd(lambda_i)
  s <- svdi$d
  s[s < 2*.Machine$double.eps * max(s)]  <- 2 * .Machine$double.eps * max(s)
  lambda <- svdi$u %*% diag(s) %*% svdi$v
  lambdaInv <- t(svdi$v) %*% diag(1/s) %*% t(svdi$u)
  return(list(lambda=lambda, lambdaInv=lambdaInv))
}
