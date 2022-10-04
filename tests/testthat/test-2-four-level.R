require(testthat)
skip_on_cran()

# for sleepstudy, calls in here to lmer
library(lme4)
# for grad function

options(width = 500)
options(useFancyQuotes = FALSE)

### Data Used: sleepstudy
tolerance <- 1E-3
data("sleepstudy")

context("four level model")
test_that("unweighted four level model v lmer", {
  skip_on_cran()
  sleepstudy2 <- sleepstudy
  sleepstudy2$Group <- 1
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("310", "309", "349", "335"), 2, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("330", "352", "337", "369"), 3, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("331", "332", "334", "372"), 4, sleepstudy2$Group)
  sleepstudy2$w1 <- 1 
  sleepstudy2$w2 <- 1
  sleepstudy2$w3 <- 1
  sleepstudy2$w4 <- 1
  sleepstudy2$Subject <- as.character(sleepstudy2$Subject)
  set.seed(2)
  for(i in 1:20) {
    sleepstudyTmp <- sleepstudy2
    sleepstudyTmp$superGroup <- LETTERS[i]
    sleepstudyTmp$Subject <- paste0(LETTERS[i],"_",sleepstudy2$Subject)
    sleepstudyTmp$Group <- paste0(LETTERS[i],"_",sleepstudy2$Group)
    sleepstudyTmp$Reaction <- sleepstudyTmp$Reaction + 36 * rnorm(1) + 31*rnorm(nrow(sleepstudyTmp))
    if(i == 1) {
      sleepstudy3 <- sleepstudyTmp
    } else {
      sleepstudy3 <- rbind(sleepstudy3, sleepstudyTmp)
    }
  }
  
  lmr <- lmer(Reaction ~ Days + (1|Subject) + (1|Group) + (1|superGroup), data=sleepstudy3, REML=FALSE)
  wm0 <- mix(Reaction ~ Days + (1|Subject) + (1|Group) + (1|superGroup), data=sleepstudy3, weights=c("w1", "w2","w3", "w4"))
  expect_equal(wm0$lnl, as.numeric(logLik(lmr)), tol=1e-3)
  expect_equal(coef(wm0), fixef(lmr), tol=1e-4)
  expect_equal(unname(wm0$vars[length(wm0$vars)]), unname(lmr@devcomp$cmp["sigmaML"]^2), tol=1e-4)
})


test_that("Weighted v unweighted replicated four level model", {
  skip_on_cran()
  sleepstudy2 <- sleepstudy
  sleepstudy2$Group <- 1
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("310", "309", "349", "335"), 2, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("330", "352", "337", "369"), 3, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("331", "332", "334", "372"), 4, sleepstudy2$Group)
  sleepstudy2$w1 <- 1
  sleepstudy2$w2 <- 1
  sleepstudy2$w3 <- 1
  sleepstudy2$w4 <- 1
  sleepstudy2$Subject <- as.character(sleepstudy2$Subject)
  set.seed(2)
  for(i in 1:20) {
    sleepstudyTmp <- sleepstudy2
    sleepstudyTmp$superGroup <- LETTERS[i]
    sleepstudyTmp$Subject <- paste0(LETTERS[i],"_",sleepstudy2$Subject)
    sleepstudyTmp$Group <- paste0(LETTERS[i],"_",sleepstudy2$Group)
    sleepstudyTmp$Reaction <- sleepstudyTmp$Reaction + 36 * rnorm(1) + 31*rnorm(nrow(sleepstudyTmp))
    if(i == 1) {
      sleepstudy3 <- sleepstudyTmp
    } else {
      sleepstudy3 <- rbind(sleepstudy3, sleepstudyTmp)
    }
  }
  # all weights must be 1 or 2
  w10 <- sample(c(rep(2,600), rep(1,3600-600)), 3600)
  sleepstudyrep <- sleepstudy3
  sleepstudy3$w1 <- w10
  for(i in 1:3600) {
    if(w10[i] > 1) {
      sst <- sleepstudy3[i,]
      sleepstudyrep <- rbind(sleepstudyrep, sst)
    }
  }
  ug <- unique(sleepstudy3$Subject)
  w20 <- sample(c(rep(2,120), rep(1,length(ug)-120)), 360)
  for(i in 1:length(ug)) {
    if(w20[i] > 1) {
      sst <- sleepstudyrep[sleepstudyrep$Subject == ug[i],]
      sst$Subject <- paste0("r2_", sst$Subject)
      sleepstudyrep <- rbind(sleepstudyrep, sst)
      sleepstudy3[sleepstudy3$Subject == ug[i],"w2"] <- w20[i]
    }
  }
  ug <- unique(sleepstudy3$Group)
  w30 <- sample(c(rep(2,20), rep(1,length(ug)-20)), length(ug))
  for(i in 1:length(ug)) {
    if(w30[i] > 1) {
      sst <- sleepstudyrep[sleepstudyrep$Group == ug[i],]
      sst$Subject <- paste0("r3_", sst$Subject)
      sst$Group <- paste0("r3_", sst$Group)
      sleepstudyrep <- rbind(sleepstudyrep, sst)
      sleepstudy3[sleepstudy3$Group == ug[i],"w3"] <- w30[i]
    }
  }
  ug <- unique(sleepstudy3$superGroup)
  w40 <- sample(c(rep(2,5), rep(1,length(ug)-5)), length(ug))
  for(i in 1:length(ug)) {
    if(w40[i] > 1) {
      sst <- sleepstudyrep[sleepstudyrep$superGroup == ug[i],]
      sst$Subject <- paste0("r4_", sst$Subject)
      sst$Group <- paste0("r4_", sst$Group)
      sst$superGroup <- paste0("r4_", sst$superGroup)
      sleepstudyrep <- rbind(sleepstudyrep, sst)
      sleepstudy3[sleepstudy3$superGroup == ug[i], "w4"] <- w40[i]
    }
  }
  sleepstudyrep$w1 <- sleepstudyrep$w2 <- sleepstudyrep$w3 <- sleepstudyrep$w4 <- 1
  sleepstudy3$w3u <- sleepstudy3$w3 * sleepstudy3$w4
  sleepstudy3$w2u <- sleepstudy3$w2 * sleepstudy3$w3u
  sleepstudy3$w1u <- sleepstudy3$w1 * sleepstudy3$w2u
  
  suppressWarnings(lmr <- lmer(Reaction ~ Days + (1|Subject) + (1|Group) + (1|superGroup), data=sleepstudyrep, REML=FALSE))
  wm0 <- mix(Reaction ~ Days + (1|Subject) + (1|Group) + (1|superGroup), data=sleepstudy3, weights=c("w1", "w2","w3", "w4"), cWeights=TRUE)
  wmr <- mix(Reaction ~ Days + (1|Subject) + (1|Group) + (1|superGroup), data=sleepstudyrep, weights=c("w1", "w2","w3", "w4"))
  summary(lmr)
  summary(wm0)
  summary(wmr)
  expect_equal(wm0$lnl, wmr$lnl, tol=1e-3)
  expect_equal(wm0$lnl, as.numeric(logLik(lmr)), tol=1e-3)
  expect_equal(coef(wm0), fixef(lmr), tol=1e-4)
  expect_equal(coef(wm0), coef(wmr), tol=1e-4)
  expect_equal(unname(wm0$vars[length(wm0$vars)]), unname(lmr@devcomp$cmp["sigmaML"]^2), tol=1e-4)
  expect_equal(wm0$vars, wmr$vars, tol=1e-4)
})