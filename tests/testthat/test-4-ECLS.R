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

context("Model with top level groups that have entirely 0 columns in Z")
test_that("Model with top level groups that have entirely 0 columns in Z", {
  skip_on_cran()
  require(EdSurvey)
  downloadECLS_K(root=edsurveyHome, years=2011, verbose=FALSE)
  ee <- readECLS_K2011(paste0(edsurveyHome, "ECLS_K/2011/"), verbose=FALSE)
  gg <- EdSurvey::getData(c("x2rscalk5", "childid", "s2_id", "w1_2p0", "x3sumsh", "p1chldbk","p2freerd"), data=ee, dropOmittedLevels=FALSE, returnJKreplicates=FALSE)
  gg$frpl <- ifelse(gg$p2freerd %in% c("1: FREE LUNCH", "2: REDUCED PRICE LUNCH"), 1, 0)
  gg$w1 <- gg$w1_2p0
  gg$w2 <- 1
  gg$n <- ave(gg$s2_id,gg$s2_id, FUN=length)
  gg2 <- gg[!is.na(gg$x2rscalk5) & gg$w1>0 & !is.na(gg$p1chldbk) & gg$n > 15 & gg$s2_id < 1500,]
  m3 <- mix(x2rscalk5 ~ p1chldbk + frpl + (1+frpl|s2_id), data=gg2, weights=c("w1", "w2"), verbose=FALSE)
  # regression tests
  expect_equal(m3$lnl, -4076916.913043, tol=1e-5)
  expect_equal(m3$coef, c(`(Intercept)` = 68.9855317534602, p1chldbk = 0.00849589076218699, frpl = -4.09352553895415), tol=1e-5)
  expect_equal(m3$SE, c(`(Intercept)` = 0.5581593716085, p1chldbk = 0.0032088689344584, frpl = 0.648499658236494), tol=1e-5)
  varDF0 <- structure(list(grp = c("s2_id", "s2_id", "s2_id", "Residual"),
                           var1 = c("(Intercept)", "frpl", "(Intercept)", NA),
                           var2 = c(NA, NA, "frpl", NA),
                           vcov = c(60.7787852874039, 92.4431165557077, -41.5862322447801, 146.8701875209),
                           ngrp = c(267, 267, 267, 4157),
                           level = c(2, 2, 2, 1),
                           SEvcov = c(6.9518, 11.6894, 6.6054, 6.470529),
                           fullGroup = c("s2_id.(Intercept)", "s2_id.frpl", "s2_id.(Intercept)", "Residual")),
                      row.names = c(NA, -4L), class = "data.frame")
  # large error in SEvcov only
  expect_equal(m3$varDF, varDF0, tol=1e-2)
  m3$varDF$SEvcov <- NULL
  varDF0$SEvcov <- NULL
  expect_equal(m3$varDF, varDF0, tol=1e-4)
})

context("ECLSK three level unordered model")
test_that("ECLSK three level unordered model", {
  skip_on_cran()
  skip_if_not_installed("EdSurvey")
  skip_if_not_installed("tidyr")
  require(EdSurvey)
  require(tidyr)
  
  eclsk11 <- readECLS_K2011(path = paste0(edsurveyHome, "ECLS_K/2011"))
  
  myDataWide <- EdSurvey::getData(eclsk11, c("childid", "x1mscalk5", "x2mscalk5",
                                             "x3mscalk5","x4mscalk5", "x5mscalk5", 
                                             "x6mscalk5", "x7mscalk5", "x8mscalk5",
                                             "x9mscalk5", "w8cf8p_80", "s2_id", "w2sch0"),
                                  returnJKreplicates=FALSE)
  myDataWide <- subset(myDataWide, w8cf8p_80 > 0)
  myDataWide <- subset(myDataWide, w2sch0 > 0)
  
  myDataWide$nsch <- ave(rep(1, nrow(myDataWide)), myDataWide$s2_id, FUN=sum)
  
  # require n students per school; filter by school id for speed
  myDataWide <- subset(myDataWide, nsch == 10 & s2_id < 1800)
  
  myDataTall <- gather(data=myDataWide, key="scorevar", value="score",
                       c("x1mscalk5", "x2mscalk5", "x3mscalk5", 
                         "x4mscalk5", "x5mscalk5", "x6mscalk5", 
                         "x7mscalk5", "x8mscalk5", "x9mscalk5") )
  
  myDataTall$wave <- substr(myDataTall$scorevar, 2, 2)
  
  myDataTall$calYear <- ifelse(myDataTall$wave==1, 0, NA) # fall 2010 (October)
  myDataTall$calYear <- ifelse(myDataTall$wave==2, 6, myDataTall$calYear) # Spring 2011 (April)
  myDataTall$calYear <- ifelse(myDataTall$wave==3, 12, myDataTall$calYear) # Fall 2011 
  myDataTall$calYear <- ifelse(myDataTall$wave==4, 18, myDataTall$calYear) # Spring 2012
  myDataTall$calYear <- ifelse(myDataTall$wave==5, 24, myDataTall$calYear) # Fall 2012
  myDataTall$calYear <- ifelse(myDataTall$wave==6, 30, myDataTall$calYear) # Spring 2013
  myDataTall$calYear <- ifelse(myDataTall$wave==7, 42, myDataTall$calYear) # Spring 2014
  myDataTall$calYear <- ifelse(myDataTall$wave==8, 54, myDataTall$calYear) # Spring 2015
  myDataTall$calYear <- ifelse(myDataTall$wave==9, 66, myDataTall$calYear) # Spring 2016???
  
  table(myDataTall$wave, myDataTall$calYear,dnn=c("Wave","Calendar Year"), useNA="ifany")
  
  myDataTall$acaYear <- ifelse(myDataTall$wave==1, 0, NA) # fall 2010 (October)
  myDataTall$acaYear <- ifelse(myDataTall$wave==2, 6, myDataTall$acaYear) # Spring 2011 (April)
  myDataTall$acaYear <- ifelse(myDataTall$wave==3, 12-3, myDataTall$acaYear) # Fall 2011 
  myDataTall$acaYear <- ifelse(myDataTall$wave==4, 18-3, myDataTall$acaYear) # Spring 2012
  myDataTall$acaYear <- ifelse(myDataTall$wave==5, 24-6, myDataTall$acaYear) # Fall 2012
  myDataTall$acaYear <- ifelse(myDataTall$wave==6, 30-6, myDataTall$acaYear) # Spring 2013
  myDataTall$acaYear <- ifelse(myDataTall$wave==7, 42-9, myDataTall$acaYear) # Spring 2014
  myDataTall$acaYear <- ifelse(myDataTall$wave==8, 54-12, myDataTall$acaYear) # Spring 2015
  myDataTall$acaYear <- ifelse(myDataTall$wave==9, 66-12, myDataTall$acaYear) # Spring 2016???
  
  myDataTall$w1c <- 1
  myDataTall$w2 <- myDataTall$w8cf8p_80
  myDataTall$w3 <- myDataTall$w2sch0
  myDataTall$w1 <- myDataTall$w1c * myDataTall$w2
  
  #write.csv(myDataTall, "myDataTall.csv")
  
  # check ranef orgering on complicated data
  suppressWarnings(lmu <- lmer(score ~ calYear + (1+calYear|childid) + (1|s2_id), data=myDataTall, verbose=FALSE, REML=FALSE))
  myDataTall$one <- 1
  suppressWarnings(mu <- mix(score ~ calYear + (1+calYear|childid) + (1|s2_id), data=myDataTall, verbose=FALSE, weights=c("one", "one", "one")))
  lmeRanef <- ranef(lmu)$childid
  # not in WeMix output
  attr(lmeRanef, "postVar") <- NULL
  expect_equal(lmeRanef, mu$ranefMat$childid, 100*(.Machine$double.eps)^0.25)
  lmeRanef <- ranef(lmu)$s2_id
  # not in WeMix output
  attr(lmeRanef, "postVar") <- NULL
  expect_equal(lmeRanef, mu$ranefMat$s2_id, 10* (.Machine$double.eps)^0.25)
  
  # sometimes this gives a warning, but not always. The important part is the results.
  suppressWarnings(m1 <- mix(score ~ calYear + (1+calYear|childid) + (1|s2_id), data=myDataTall, verbose=FALSE, weights=c("w1", "w2", "w3")))
  expect_equal(m1$coef, c(`(Intercept)` = 49.6594791871134, calYear = 1.23043346166935), tolerance=200*sqrt(.Machine$double.eps))  
  expect_equal(m1$lnl, -6258739.70379471)
  
  sumRef <- structure(c(49.6594791871134, 1.23043346166935, 2.7913905771771, 
                        0.020424766920315, 17.790229569856, 60.242227804594),
                      .Dim = 2:3,
                      .Dimnames = list(c("(Intercept)", "calYear"),
                                       c("Estimate", "Std. Error", "t value")))
  
  expect_equal(summary(m1)$coef, sumRef, tolerance = (.Machine$double.eps)^0.25)
  
  sumVarDF <- structure(list(grp = c("childid", "childid", "childid", "s2_id", "Residual"),
                             var1 = c("(Intercept)", "calYear", "(Intercept)", "(Intercept)", NA),
                             var2 = c(NA, NA, "calYear", NA, NA),
                             vcov = c(87.9540891469932, 0.0203169310526866, 0.682227179260158, 52.321204750664, 79.1772622305314),
                             ngrp = c(120, 120, 120, 12, 1080),
                             level = c(2, 2, 2, 3, 1),
                             SEvcov = c(20.6124219047353, 0.00775092348656896, 0.234090937293996, 19.9393913552464, 8.61080027941434),
                             fullGroup = c("childid.(Intercept)", "childid.calYear", "childid.(Intercept)", "s2_id.(Intercept)", "Residual")),
                        row.names = c(NA, -5L),
                        class = "data.frame")
  expect_equal(m1$varDF, sumVarDF, tolerance = 2 * (.Machine$double.eps)^0.25) 
})