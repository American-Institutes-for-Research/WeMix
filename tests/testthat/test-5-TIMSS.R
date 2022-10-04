require(testthat)
skip_on_cran()

options(width = 500)
options(useFancyQuotes = FALSE)

if(!exists("edsurveyHome")) {
  if (Sys.info()[['sysname']] == "Windows") {
    edsurveyHome <- "C:/EdSurveyData/"
  } else {
    edsurveyHome <- "~/EdSurveyData/"
  }
}

context("TIMSS tests")
test_that("TIMSS tests", {
  skip_on_cran()
  # original version by Christian Kjeldsen
  try(downloadTIMSS(root=edsurveyHome, years=2015, cache=FALSE, verbose=FALSE), silent=TRUE)
  dnk15 <- readTIMSS(file.path(edsurveyHome,"/TIMSS/2015"), countries="dnk", gradeLvl=4, verbose=FALSE)
  dnk15dat <- EdSurvey::getData(data=dnk15, varnames=c("atbg01", "asbgsb", "mmat", "asbghrl", "matwgt", "idschool","schwgt"))
  dnk15dat <- subset(dnk15dat, matwgt > 0 & schwgt > 0)
  dnk15dat$cwt2_math <- dnk15dat$schwgt
  dnk15dat$cwt1_math <- dnk15dat$matwgt/dnk15dat$schwgt
  # variance estimation requires a matrix singular by base standards but not Matrix standards
  mm2 <- mix(asmmat01 ~ atbg01 + asbghrl + (1|idschool), data = dnk15dat, weights=c("matwgt","schwgt"))
  
  mm2ref <- structure(c(375.385865901528, -0.553178832770401, 14.9395466422176, 
                        13.8062259868586, 0.378593770436568, 0.994526615735166,
                        27.1896075190162, -1.46114087437972, 15.0217665428432),
                      .Dim = c(3L, 3L), .Dimnames = list(c("(Intercept)", "atbg01", "asbghrl"),
                                                         c("Estimate", "Std. Error", "t value")))
  expect_equal(summary(mm2)$coef, mm2ref)
  # glm Wald test
  mm1 <- mix(I(asmmat01>450) ~ atbg01 + (1|idschool), data=dnk15dat, weights=c("cwt1_math","cwt2_math"),
             cWeights=TRUE, family=binomial(link="logit"))
  mm1s <- summary(mm1)
  mm1w <- WeMix::waldTest(fittedModel=mm1, type="beta", coefs="atbg01")
  
  expect_equal(mm1s$coef[2,3]^2, mm1w$Wald)
})