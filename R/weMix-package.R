#' @title Estimate Weighted Mixed-Effects Models 
#' 
#' @description The WeMix package estimates mixed-effects models (also called multilevel models, 
#' mixed models, or HLMs) with survey weights.
#' 
#' @section Details:
#' 
#' This package is unique in allowing users to analyze data that may have unequal selection
#' probability at both  the individual and group 
#' levels.  For linear models, the model  is evaluated with a weighted version of the estimating equations
#'  used by Bates, Maechler, Bolker, and Walker (2015) in \code{lme4}. In the non-linear case,  WeMix uses numerical 
#'  integration (Gauss-Hermite and adaptive Gauss-Hermite  quadrature) to estimate mixed-effects models with 
#'  survey weights at all levels of the model. 
#'  Note that \code{lme4} is the preferred way to estimate such 
#' models when there are no survey weights or weights only at the lowest level, and our 
#' estimation starts with parameters estimated in lme4. WeMix is intended for use in cases 
#'  where there are weights at all levels and is only for use with fully nested data. 
#' To start using WeMix, see the vignettes covering
#' the mathematical background of mixed-effects model estimation and use the
#' \code{mix} function to estimate models. Use 
#' \code{browseVignettes(package="WeMix")} to see the vignettes.
#' 
#' @references Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting Linear Mixed-Effects
#' Models Using lme4. Journal of Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01
#' 
#' Rabe-Hesketh, S., & Skrondal, A. (2006) Multilevel Modelling of Complex Survey Data. Journal
#' of the Royal Statistical Society: Series A (Statistics in Society), 169, 805-827.
#' https://doi.org/10.1111/j.1467-985X.2006.00426.x
#' 
#' Bates, D. & Pinheiro, J. C. (1998). Computational Methods for Multilevel Modelling. Bell labs working paper.
#' @name WeMix-package
"_PACKAGE"
