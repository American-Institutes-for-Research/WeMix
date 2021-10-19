require(testthat)

# for sleepstudy, calls in here to lmer
library(lme4)
# for grad function

options(width = 500)
options(useFancyQuotes = FALSE)

### Data Used: sleepstudy
tolerance <- 2E-5
data("sleepstudy")
### Unweigted =====================================
sleepstudyU <- sleepstudy
sleepstudyU$weight1L1 <- 1
sleepstudyU$weight1L2 <- 1

context("The model runs")
test_that("The model runs", {
  system.time(wm0 <- mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2")))
  lme1 <- lmer(Reaction ~ Days + (1|Subject), data=sleepstudyU, REML=FALSE)
  expect_equal(wm0$lnl, -897.039321502613, tolerance=tolerance*897) # value from  lmer(Reaction ~ Days + (1 | Subject), data=sleepstudy, REML=FALSE)
  expect_equal(unname(c(wm0$coef)), unname(fixef(lme1)), tolerance=tolerance)
  expect_equal(unname(wm0$theta), unname(lme1@theta), tolerance=tolerance)
  expect_equal(coef(wm0),
               expected = getME(lme1, "fixef"),
               tolerance = tolerance)
  lmevars1 <- data.frame(summary(lme1)$varcor)$sdcor
  expect_equal(unname(wm0$vars),
               expected = unname(lmevars1)^2,
               tolerance = tolerance)
  # agrees with GLLAMM
  #source: logit_command.do
  gllamm_model1 <- c("(Intercept)" = 251.4051,
                     "Days"         = 10.46729,
                     "Residual"     = 1296.8768,
                     "Subject"      = 954.52789,
                     "lnl"          = -897.03932)
  expect_equal(unname(coef(wm0)),
               expected = unname(gllamm_model1[1:2]),
               tolerance = abs(tolerance))
  expect_equal(unname(wm0$vars),
               expected = unname(gllamm_model1[3:4]),
               tolerance = tolerance)
  expect_equal(unname(wm0$lnl),
               expected=unname(gllamm_model1[5]),
               tolerance=abs(tolerance))
  
  sleepstudy2 <- sleepstudy
  sleepstudy2$block <- rep(1:6,each=30)
  sleepstudy2$weight1L1 <- sleepstudy2$weight1L2 <- sleepstudy2$weight1L3 <- 1
  # change Reaction by the artifical block we just added
  sleepstudy2$Reaction <- sleepstudy2$Reaction + sleepstudy2$block*20
  names(sleepstudy2) <- gsub("weight1L","pwt",names(sleepstudy2))
  lmr <- lmer(Reaction ~ Days + (1|Subject) + (1|block), data=sleepstudy2, REML=FALSE)
  mm <- mix(Reaction ~ Days + (1|Subject) + (1|block), data=sleepstudy2, weights=c("pwt1","pwt2","pwt3"))
  expect_equal(unname(mm$lnl),
               expected=unname(logLik(lmr)[[1]]),
               tolerance=tolerance)
})


context("Factor binomial")
test_that("Factor binomial", {
  sleepstudyM <- sleepstudyU
  sleepstudyM$highR <- factor(ifelse(sleepstudyM$Reaction > 340, 2, 1), levels=c(1,2), labels=c("L","H"))
  sleepstudyM$Sub <- factor(sleepstudyM$Subject, levels=300:400) 
  m1 <- mix(highR ~ Days + (1|Subject), data=sleepstudyM, weights=c("weight1L1", "weight1L2"), family="binomial")
  expect_equal(coef(m1), c(`(Intercept)` = -8.69336570860094, Days = 1.17236933934874), tol=1e-6)
  summaryREF <- c("Call:",
                  "mix(formula = highR ~ Days + (1 | Subject), data = sleepstudyM, ", 
                  "    weights = c(\"weight1L1\", \"weight1L2\"), family = \"binomial\")", 
                  "",
                  "Variance terms:",
                  " Level   Group        Name Variance Std. Error Std.Dev.", 
                  "     2 Subject (Intercept)     7.99       5.13     2.83", "Groups:", 
                  " Level   Group n size mean wgt sum wgt",
                  "     2 Subject     18        1      18",
                  "     1     Obs    180        1     180",
                  "",
                  "Fixed Effects:",
                  "            Estimate Std. Error t value",
                  "(Intercept)   -8.693      2.072   -4.20", 
                  "Days           1.172      0.278    4.21",
                  "",
                  "lnl= -53.90 ")
  withr::with_options(list(digits=2),
                       co <- capture.output(summary(m1))
                     )
  expect_equal(co, summaryREF)
})

context("Unweighted model with 2 random effects")
test_that("Agrees with lme4 3,handles missing data", {
  sleepstudyM <- sleepstudyU
  #introduce a missing value 
  sleepstudyM$Days[3] <- NA
  lme2 <- lmer(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), data=sleepstudyM, REML=FALSE)
  # the dropped row should cause a warning
  expect_warning(wm2 <- mix(Reaction ~ Days + (1 | Subject) + (0 + Days | Subject), data=sleepstudyM, weights=c("weight1L1", "weight1L2")),"with missing data")
  expect_equal(wm2$lnl, as.numeric(logLik(lme2)), tol=1E-7*abs(as.numeric(logLik(lme2))))

    # check coef
  expect_equal(coef(wm2),
               expected = getME(lme2, "fixef"),
               tolerance = tolerance)
  # check vars
  lmewm2vars <- data.frame(summary(lme2)$varcor)$sdcor
  expect_equal(unname(wm2$vars),
               expected = unname(lmewm2vars)^2,
               tolerance = tolerance)
})
 
context("Mean Centering Matches HLM results")
test_that("Mean Centering Matches HLM results", {
 wm1 <- mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), nQuad=13,center_group=list("Subject"= as.formula(~Days)), verbose=FALSE)
 expect_equal(unname(wm1$coef), c(298.507892, 10.467286), tolerance=1E-1)
 expect_equal(wm1$lnl, -897.0393, tolerance=tolerance*897) # value from  lmer(Reaction ~ Days + (1 | Subject), data=sleepstudy, REML=FALSE)
})

ss <- sleepstudy
ss1 <- ss
ss2 <- ss
doubles <- c(308, 309, 310) # subject with double obs

ss2 <- rbind(ss, subset(ss, Subject %in% doubles))

ss$w <- ifelse(ss$Subject %in% doubles, 2, 1)
contrasts(ss2$Subject) <- "contr.sum"
ss1$W1 <- ifelse(ss1$Subject %in% doubles, 2, 1)
ss1$W2 <- 1
ss1$bin <- ifelse(sleepstudy$Reaction<300,0,1) #for the binomial test

ss2$W2 <- ss2$W1 <- 1

doubles <- c(308, 309, 310) # subject with double obs
ss30 <- subset(ss, Subject %in% doubles)
ss30$Subject <- as.numeric(as.character(ss30$Subject)) + 1000
ss0 <- ss
ss0$Subject <- as.numeric(as.character(ss$Subject))
ss3 <- rbind(ss0, ss30)
ss3$Subject <- as.factor(ss3$Subject)

ss3$W2 <- 1
ss3$W1 <- 1

ss4 <- ss
ss4$W2 <- ifelse(ss4$Subject %in% doubles, 2, 1)
ss4$W1 <- ifelse(ss4$Subject %in% doubles, 2, 1) # make unconditional weight

context("GLM works: Binomial")
test_that("GLM works: Binomial", {
  skip_on_cran()
  #full test for binomial 
  bi_1 <- mix(bin~Days + (1|Subject), data=ss1, family=binomial(link="logit"), verbose=FALSE,
              weights=c("W1", "W2"), nQuad=13)
  expect_equal(unname(bi_1$coef), c(-3.3448,.5928), tolerance=1E-3)
  expect_equal(bi_1$lnl, -93.751679, tolerance=1E-5)
  sum_bi <-  summary(bi_1)
  expect_is(summary(bi_1), "summaryWeMixResults")
})

context("Repeating is the same as weighting: L1 replicate vs weighting")
test_that("Repeating is the same as weighting: L1 replicate vs weighting", {
  # mix for L1, weighted
  wmeL1W <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss1,
                weights=c("W1", "W2"))

  # mix for L1, duplicated
  system.time(wmeL1D <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss2,
                            weights=c("W1", "W2")))

  # check weighted agrees with duplicated lme4 results
  expect_equal(wmeL1W$lnl, -1048.34318418762, tolerance=1050*tolerance)

  # check duplicated agrees with duplicated lme4 results
  expect_equal(wmeL1D$lnl, -1048.34318418762, tolerance=1050*tolerance)
  
 })

context("grouping factor not sorted")
test_that("grouping factor not sorted", {
  skip_on_cran()
  ss1_mixed <- ss1[c(125:180,1,100,2:99,101:124),]
  row.names(ss1_mixed) <- NULL
  # mix for L1, weighted
  wmeL1W <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss1_mixed,
                weights=c("W1", "W2"), nQuad=13, run=TRUE,  verbose=FALSE)

  # mix for L1, duplicated
  system.time(wmeL1D <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss2,
                            weights=c("W1", "W2"), nQuad=13, run=FALSE,  verbose=FALSE))
  #statares <- c(251.4619, 10.40726, 1000.7466, 1338.0865) # not used

  # check weighted agrees with duplicated lme4 results
  expect_equal(wmeL1W$lnl, -1048.34318418762, tolerance=tolerance)
  expect_equal(wmeL1D$lnl, -1048.34318418762, tolerance=tolerance)

  # check final results
  suppressWarnings(mix1 <-  mix(formula=Reaction ~ Days + (1 | Subject), data=ss1_mixed,
                                weights=c("W1", "W2")))
  suppressWarnings(mix1REF <-  mix(formula=Reaction ~ Days + (1 | Subject), data=ss1,
                                   weights=c("W1", "W2")))
  expect_equal(mix1$coef, mix1REF$coef, tolerance=1e3)
  expect_equal(mix1$vars, mix1REF$vars, tolerance=1e-3)
  expect_equal(mix1$lnl, mix1REF$lnl, tolerance=1e-3)
  #check  weighted fixed effects variances
  expect_equal(unname(sqrt(diag(mix1$cov_mat))), unname(sqrt(diag(mix1REF$cov_mat))),tolerance = tolerance)
  
  #
  
})

context("Weighted three level model unsorted")
test_that("Weighted three level model unsorted", {
  skip_on_cran()
  sleepstudy2 <- sleepstudy
  sleepstudy2$Group <- 1
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("310", "309", "349", "335"), 2, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("330", "352", "337", "369"), 3, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("331", "332", "334", "372"), 4, sleepstudy2$Group)
  sleepstudy2$Group <- factor(sleepstudy2$Group)
  ss2 <- sleepstudy2
  w1c <- w1 <- rep(1,180)
  w2c <- w2 <- rep(1,18)
  w3c <- w3 <- rep(1,4)
  # unbalanced (non-identical within a Subject), level-1 (obs level) weights
  w1c[1:4] <- w1[1:4] <- 2
  uR <- sleepstudy2[1:4,]
  sleepstudy2 <- rbind(sleepstudy2, uR)
  # level-2 weights
  w2c[1] <- w2[1] <- 2
  w1[ss2$Subject=="308"] <- 2*w1[ss2$Subject=="308"]
  sR <- subset(sleepstudy2, Subject == "308")
  sR$Subject <- "S308"
  sleepstudy2 <- rbind(sleepstudy2, sR)
  # level-3 weights
  w3c[1] <- w3[1] <- 2
  w2[c(1,7,12,13,16,17)] <- 2*w2[c(1,7,12,13,16,17)]
  w1[ss2$Group==1] <- 2*w1[ss2$Group==1]
  
  gR <- subset(sleepstudy2, Group %in% 1)
  gR$Subject <- paste0("G", gR$Subject)
  gR$Group <- "11"
  sleepstudy2 <- rbind(sleepstudy2, gR)

  # lmr for reference
  lmr <- lmer(Reaction ~ Days + (1|Subject) + (1|Group), 
              data=sleepstudy2,  control=lmerControl(optimizer="bobyqa"),
              REML=FALSE)
  
  ss2$w1 <- w1
  ss2$w2 <- rep(w2,each=10)
  ss2$w3 <- ifelse(ss2$Subject %in% c("308", "333", "350", "351", "370", "371"),2,1)

  ss3 <- ss2[sample(row.names(ss2), size=nrow(ss2)), ]
  wm0 <- mix(Reaction ~ Days + (1|Subject) + (1|Group), data=ss3, 
             weights=c("w1", "w2", "w3"))

  # check coef
  expect_equal(coef(wm0),
               expected = getME(lmr, "fixef"),
               tolerance = tolerance)
  # check vars
  lmewm2vars <- data.frame(summary(lmr)$varcor)$sdcor
  expect_equal(unname(wm0$vars),
               expected = unname(lmewm2vars)^2,
               tolerance = tolerance)
  
  # var beta not expected to be equal  
})

context("repeating is the same as weighting: L2 replicate vs weighting")
test_that("L2 replicate vs weighting", {
  # mix for L2, duplicated
  system.time(wmeL2D <- mix(formula=Reaction ~ Days + (1 | Subject),
                            data=ss3, weights=c("W1", "W2")))

  # mix for L2, weighted
  system.time(wmeL2W <- mix(formula=Reaction ~ Days + (1 | Subject), data=ss4,
                            weights=c("W1", "W2")))

  expect_equal(wmeL2W$lnl, -1055.34690957995, tolerance=2E-7)
  expect_equal(wmeL2D$lnl, -1055.34690957995, tolerance=2E-7)
})

context("Repeating is the same as weighting: L1 replicate vs weighting, 2 REs")
test_that("Repeating is the same as weighting: L1 replicate vs weighting, 2 REs", {
  # mix for L1, weighted, 2 REs
  wmeL1WRE2 <- mix(formula=Reaction ~ Days + (1 | Subject) + (0+Days|Subject),
                   data=ss1, weights=c("W1", "W2"), nQuad=13, run=FALSE, verbose=FALSE)

  # mix for L1, duplicated, 2 REs
  wmeL1DRE2 <- mix(formula=Reaction ~ Days + (1 | Subject) + (0+Days|Subject),
                   data=ss2, weights=c("W1", "W2"),nQuad=13, run=FALSE, verbose=FALSE)

  expect_equal(wmeL1WRE2$lnl, -1018.29298875158, tolerance=1050*2E-7)

  expect_equal(wmeL1DRE2$lnl, -1018.29298875158, tolerance=1050*2E-7)
})

ssB <- sleepstudy
set.seed(2)
ssB$Reaction <- ssB$Days * 3.141 + rnorm(nrow(ssB))
ssB$W2 <- 1
ssB$W1 <- 1

context("Zero variance estimate")
test_that("Zero variance estimate", {
  skip_on_cran()
  # this has 0 variance estimate in lmer
  lmeB <- lmer(Reaction ~ Days + (1|Subject), data=ssB, REML=FALSE)
  mixB <- mix(formula=Reaction ~ Days + (1 | Subject), data=ssB,
              weights=c("W1", "W2"),  nQuad=13, run=TRUE,   verbose=FALSE)
  expect_equal(mixB$lnl, as.numeric(logLik(lmeB)), tol=1e-5)
  expect_equal(coef(mixB), fixef(lmeB), tol=1e-5)
  expect_equal(unname(mixB$vars[length(mixB$vars)]), unname(lmeB@devcomp$cmp["sigmaML"]^2), tol=1e-5)

  ss1 <- sleepstudy
  #add group variables for 3 level model 
  ss1$Group <- 1
  ss1$Group <- ifelse(ss1$Subject %in% c(349,335,330, 352, 337, 369), 2, ss1$Group)

  # Create weights
  ss1$W1 <- ifelse(ss1$Subject %in% c(308, 309, 310), 2, 1)
  ss1$W2 <- 1
  ss1$W3 <- ifelse(ss1$Group == 2, 2, 1)

  #Run three level model with random slope and intercept. 
  three_level <- mix(Reaction ~ Days + (1|Subject) + (1+Days|Group), data=ss1, weights = c("W1","W2","W3"))
  expect_is(three_level, "WeMixResults")
})

context("Unweighted three level model")
test_that("Unweighted three level model", {
  sleepstudy2 <- sleepstudy
  sleepstudy2$Group <- 1
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("310", "309", "349", "335"), 2, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("330", "352", "337", "369"), 3, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("331", "332", "334", "372"), 4, sleepstudy2$Group)
  sleepstudy2$Group <- factor(sleepstudy2$Group)
  sleepstudy2$w1 <- 1 
  sleepstudy2$w2 <- 1
  sleepstudy2$w3 <- 1
  wm0 <- mix(Reaction ~ Days + (1|Subject) + (0+Days|Subject) + (1 | Group), data=sleepstudy2, weights=c("w1", "w2","w3"), verbose=FALSE, run=TRUE)
  lm0 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject) + (1 | Group), data=sleepstudy2,REML=FALSE)
  # check vars
  lmevars1 <- data.frame(summary(lm0)$varcor)$sdcor
  expect_equal(unname(wm0$vars),
               expected = unname(lmevars1)^2,
               tolerance = 1e-3)
  expect_equal(wm0$lnl, as.numeric(logLik(lm0)), tol=1e-3)
  expect_equal(coef(wm0), fixef(lm0), tol=1e-4)
  expect_equal(unname(wm0$vars[length(wm0$vars)]), unname(lm0@devcomp$cmp["sigmaML"]^2), tol=1e-4)
})

context("Three level model slash and colon")
test_that("Three level model slash and colon", {
  sleepstudy2 <- sleepstudy
  sleepstudy2$Group <- 1
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("310", "309", "349", "335"), 2, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("330", "352", "337", "369"), 3, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("331", "332", "334", "372"), 4, sleepstudy2$Group)
  sleepstudy2$Group <- factor(sleepstudy2$Group)
  # make subj that restarts at 1 per group, so there are four Subjects with subj=1, one per group.
  sleepstudy2 <- sleepstudy2[order(sleepstudy2$Group, sleepstudy2$Subject),]
  sleepstudy2$subj <- factor(rep(c(1:6,1:4,1:4,1:4),each=10))
  # table(sleepstudy2$subj, sleepstudy2$Group) # shows each group has a subject 1:4 and group 1 as a 5 and 6 too
  ss2 <- sleepstudy2
  ss2$w1 <- w1c <- w1 <- rep(1,180)
  w2c <- w2 <- rep(1,18)
  w3c <- w3 <- rep(1,4)
  ss2$w2 <- rep(1, 180)
  ss2$w3 <- rep(1, 180)

  lmr <- lmer(Reaction ~ Days + (1|Group) + (1|Subject), data=sleepstudy2, REML=FALSE)
  
  # next line is a bad specification: confounds group=1, subj=1 person with group=2, subj=2 person
  # lmr <- lmer(Reaction ~ Days + (1|Group) + (1|subj), data=sleepstudy2, REML=FALSE)
  wm0 <- mix(Reaction ~ Days + (1|Group:subj) + (1|Group) , data=ss2, weights=c("w1", "w2","w3"))
  expect_equal(wm0$lnl, as.numeric(logLik(lmr)), tol=1e-6)
  expect_equal(coef(wm0), fixef(lmr), tol=1e-6)
  expect_equal(unname(wm0$vars[length(wm0$vars)]), unname(lmr@devcomp$cmp["sigmaML"]^2), tol=1e-4)
  
  #group mean center Days to test group mean centering as well
  # we know the mean is 4.5 for all groups because all individuals were tested on all days
  sleepstudy2$gmc_days <- sleepstudy2$Days-4.5
  
  lmr <- lmer(Reaction ~ gmc_days + (1|Group/subj), data=sleepstudy2, REML=FALSE)
  # above is same as:
  # lmr <- lmer(Reaction ~ Days + (1|Group) + (1|Group:subj), data=sleepstudy2, REML=FALSE)
  wm0 <- mix(Reaction ~ Days + (1|Group/subj), data=ss2, weights=c("w1", "w2","w3"),center_group = list("Group/subj" = ~Days))
  expect_equal(wm0$lnl, as.numeric(logLik(lmr)), tol=1e-6)
  expect_equal(unname(coef(wm0)), unname(fixef(lmr)), tol=1e-6)
  expect_equal(unname(wm0$vars[length(wm0$vars)]), unname(lmr@devcomp$cmp["sigmaML"]^2), tol=1e-4)

  # non-nested model
  #lmr <- lmer(Reaction ~ Days + (1|Subject/Group), data=sleepstudy2, REML=FALSE)
  expect_error(wm0 <- mix(Reaction ~ Days + (1|Subject/Group), data=ss2, weights=c("w1", "w2", "w3")))
})

# check the format of summary output
context("Summary output format")
test_that("summary output format", {
  skip_on_cran()
  sleepstudy2 <- sleepstudy
  sleepstudy2$Group <- 1
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("310", "309", "349", "335"), 2, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("330", "352", "337", "369"), 3, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("331", "332", "334", "372"), 4, sleepstudy2$Group)
  sleepstudy2$Group <- factor(sleepstudy2$Group)
  # make subj that restarts at 1 per group, so there are four Subjects with subj=1, one per group.
  sleepstudy2 <- sleepstudy2[order(sleepstudy2$Group, sleepstudy2$Subject),]
  sleepstudy2$subj <- factor(rep(c(1:6,1:4,1:4,1:4),each=10))
  sleepstudy2$carPr <- pnorm(sleepstudy2$Reaction -285,sd=50)
  set.seed(2) 
  sleepstudy2$car <- runif(180) < sleepstudy2$carPr
  # table(sleepstudy2$subj, sleepstudy2$Group) # shows each group has a subject 1:4 and group 1 as a 5 and 6 too
  ss2 <- sleepstudy2
  ss2$w1 <- w1c <- w1 <- rep(1,180)
  w2c <- w2 <- rep(1,18)
  w3c <- w3 <- rep(1,4)
  ss2$w2 <- rep(1, 180)
  ss2$w3 <- rep(1, 180)

  co0 <- c("Call:",
          "mix(formula = Reaction ~ Days + (Days + car | Subject), data = ss2, ", 
          "    weights = c(\"w1\", \"w2\"))", "",
          "Variance terms:",
          " Level    Group        Name Variance Std. Error Std.Dev. Corr1 Corr2", 
          "     2  Subject (Intercept)    629.4      254.1    25.09            ", 
          "     2  Subject        Days     36.7       12.2     6.06  0.17      ", 
          "     2  Subject     carTRUE    648.5      309.3    25.46 -0.42 -0.39", 
          "     1 Residual                528.8      145.7    23.00            ", 
          "Groups:",
          " Level   Group n size mean wgt sum wgt",
          "     2 Subject     18        1      18",
          "     1     Obs    180        1     180",
          "",
          "Fixed Effects:",
          "            Estimate Std. Error t value", 
          "(Intercept)   252.72       6.53   38.70",
          "Days           11.19       1.50    7.44",
          "",
          "lnl= -868.13 ",
          "Intraclass Correlation= 0.713 ")

  wm0 <- mix(Reaction ~ Days + (Days+car|Subject), data=ss2, weights=c("w1", "w2"))
  withr::with_options(list(digits=2),
                       co <- capture.output(summary(wm0))
                     )
  expect_equal(co, co0)

  co1 <- c("Call:",
           "mix(formula = Reaction ~ Days + (car || Subject), data = ss2, ",
           "    weights = c(\"w1\", \"w2\"))", "", "Variance terms:",
           " Level    Group        Name Variance Std. Error Std.Dev. Corr1 Corr2", 
           "     2  Subject (Intercept)      698        461     26.4            ", 
           "     2  Subject    carFALSE      214        117     14.6     0      ", 
           "     2  Subject     carTRUE     1619        737     40.2     0  0.92", 
           "     1 Residual                  808        162     28.4            ", 
           "Groups:",
           " Level   Group n size mean wgt sum wgt",
           "     2 Subject     18        1      18",
           "     1     Obs    180        1     180",
           "",
           "Fixed Effects:",
           "            Estimate Std. Error t value",
           "(Intercept)   236.84       6.71   35.32",
           "Days            9.59       1.56    6.15",
           "",
           "lnl= -890.68 ",
           "Intraclass Correlation= 0.758 ")
  wm0 <- mix(Reaction ~ Days + (car||Subject), data=ss2, weights=c("w1", "w2"))
  withr::with_options(list(digits=2),
                       co <- capture.output(summary(wm0))
                     )
  expect_equal(co, co1)
})

context("Weighted three level model")
test_that("Weighted three level model", {
  sleepstudy2 <- sleepstudy
  sleepstudy2$Group <- 1
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("310", "309", "349", "335"), 2, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("330", "352", "337", "369"), 3, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("331", "332", "334", "372"), 4, sleepstudy2$Group)
  sleepstudy2$Group <- factor(sleepstudy2$Group)

  #sleepstudy2$NewVar <- sleepstudy2$Reaction/100 + rnorm(180,10,5)
  sleepstudy2$NewVar <- c(9.43457771946133, 12.7183652886321, 11.242490364576, 8.32049769489214, 21.4769136458094, 12.801691804759, 11.2944480319612, 23.2877967612343, 17.6554807720438, 15.6839576150421,
                          12.0415732799454, 12.6326929342708, 10.7497191412678, 8.67408880103278, 16.7872963953428, 14.6248001125594, 14.3553882674867, 11.8275049579491, 18.7525279773655, 5.57650843995549, 
                          19.7623973348143, 19.1858339464439, 16.9845555399241, 5.79214586016569, 13.6760775160985, 8.20870365747105, 10.469622698121, 9.5294364080553, 0.474272552642729, 9.90998530696626,
                          7.85598926755497, 13.6908324775899, 8.363578606421, 12.4585008913942, 20.8371544197279, 6.95572306644593, 4.49104017816579, 13.529348727477, 8.31255622061886, 8.78452997885101, 
                          17.216378291436, 13.6427335980228, 19.498592420836, 7.77906722601152, 15.1266474871889, 12.8030124499139, 6.62054433117513, 8.66974559938289, 11.9540398510289, 14.4271543736958,
                          10.0496588773239, 15.5177485124749, 9.12750122498998, 10.1801282848309, 14.2562515194416, 19.3249322019648, 10.7021090784807, 19.5336266451282, 7.58133671725254, 19.889529248361, 
                          5.32630146259616, 13.6044477039579, 4.98570236262749, 23.8801945918399, 14.447634009947, 15.6043025989833, 22.411260279495, 16.66695588071, 13.2723030644539, 17.6752877679585,
                          13.4050771535649, 16.8476633007641, 14.4754165647234, 13.2597805356783, 12.2720489082021, 16.3515751846435, 18.2589186962967, 13.5608781488381, 22.128982372665, 14.0340540905538,
                          2.08454544532827, 10.8863647522786, 21.0103204678449, 18.4595691535641, 11.9914499297614, 9.96791953740012, 20.0728999616675, 10.7237225063556, 18.6117287544885, 15.6290844263721,
                          13.584186380461, 12.8188973592142, 8.73588740878364, 22.1049445716224, 8.82190889109436, 16.7501795316447, 11.7894955061621, 8.67751729426691, 8.28560204932531, 25.3725616103908,
                          14.0620375346182, 14.9882805901214, 17.181272042065, 19.6581190406994, 11.6757977643058, 16.638959090019, 9.68436466184703, 23.1282879805077, 20.7690736926973, 14.3401954001056,
                          16.2747349705209, 8.61220295948687, 19.869302293393, 5.78478905807136, 6.59988762289215, 20.489844013367, 27.239468499143, 20.5797203990527, 19.4544929129352, 12.5320006698243, 
                          18.2664525619468, 16.6022500953422, 13.0216995414837, 12.4573456587016, 18.6922622366627, 7.37338081469625, 18.2016743123929, 12.9216600044986, 17.2543801572655, 20.1850819959469,
                          9.46342758905614, 16.1505649461445, 11.644647125078, 17.6992809329738, 13.5726089392145, 11.9006738563255, 9.11529530502925, 11.9473408244809, 13.7037965361481, 21.052192466607, 
                          8.94718178506257, 15.165337112518, 8.23261412105177, 12.2188157585207, 14.1455895553682, 4.94479286361234, 16.9113105275149, 12.7165809453375, 5.58774136619933, 8.31670629685991,
                          6.67728085523989, 11.7741431392146, 7.91590113917014, 9.13428547474213, 18.5689544076891, 18.1995195536489, 13.3510327364336, 21.4047919488063, 17.7821335119331, 9.19515401988151, 
                          10.7918240635605, 6.63952336871624, 15.3733302621749, 15.5262403960572, 13.4421915705244, 11.5980148647886, 17.3920393426159, 6.6883267354037, 10.6865661488425, 15.1104585129107,
                          9.99190979270527, 16.54240123493, 13.5687018760124, 13.5592897858269, 19.4019579232211, 15.0561754718535, 9.36036662501835, 4.27173388144395, 19.2040511963462, 6.36268714321505)

  ss2 <- sleepstudy2
  w1c <- w1 <- rep(1,180)
  w2c <- w2 <- rep(1,18)
  w3c <- w3 <- rep(1,4)
  # unbalanced (non-identical within a Subject), level-1 (obs level) weights
  w1c[1:4] <- w1[1:4] <- 2
  uR <- sleepstudy2[1:4,]
  sleepstudy2 <- rbind(sleepstudy2, uR)
  # level-2 weights
  w2c[1] <- w2[1] <- 2
  w1[ss2$Subject=="308"] <- 2*w1[ss2$Subject=="308"]
  sR <- subset(sleepstudy2, Subject == "308")
  sR$Subject <- "S308"
  sleepstudy2 <- rbind(sleepstudy2, sR)
  # level-3 weights
  w3c[1] <- w3[1] <- 2
  w2[c(1,7,12,13,16,17)] <- 2*w2[c(1,7,12,13,16,17)]
  w1[ss2$Group==1] <- 2*w1[ss2$Group==1]

  gR <- subset(sleepstudy2, Group %in% 1)
  gR$Subject <- paste0("G", gR$Subject)
  gR$Group <- "11"
  sleepstudy2 <- rbind(sleepstudy2, gR)

  ss2$w1 <- w1
  ss2$w2 <- rep(w2,each=10)
  ss2$w3 <- ifelse(ss2$Subject %in% c("308", "333", "350", "351", "370", "371"),2,1)
 
  # lmr for reference
  lmr <- lmer(Reaction ~ Days + (1|Subject) + (1|Group), data=sleepstudy2, REML=FALSE)

  #these are the conditional weights, used in the stata comparison
  ss2$c1 <- w1c
  ss2$c2 <- rep(w2c,each=10)
  ss2$c3 <- ifelse(ss2$Subject %in% c("308", "333", "350", "351", "370", "371"),2,1)

  #compare against mixed stata resutls
  wm0 <- mix(Reaction ~ Days + (1|Subject) + (1|Group), data=ss2, weights=c("w1", "w2","w3"))
  expect_equal(wm0$lnl, -1389.0983, tolerance=1e-3)
  expect_equal(unname(wm0$coef), c(243.6831,12.67954),tolerance = 1e-3)
  expect_equal(unname(wm0$vars), c(360.97 ,756.5857, 1153.704 ),tolerance = 1e-3)
  expect_equal(unname( wm0$SE), c( 11.68881 , 2.468474  ),tolerance = 1e-3)

  wm1 <- mix(Reaction ~ Days + (1 |Subject) + (1|Group)+ (0+Days|Group), data=ss2, weights=c("w1", "w2","w3"))
  expect_equal(wm1$lnl,  -1377.6876  , tolerance=1e-3)
  expect_equal(unname(wm1$coef), c(247.6234 , 11.71211 ),tolerance = 1e-3)
  expect_equal(unname(wm1$vars), c(397.6301 , 220.4182  , 16.86835 , 1022.88  ),tolerance = 1e-3)
  expect_equal(unname( wm1$SE), c(9.97378  , 2.60348  ),tolerance = 1e-3)
})

context("Wald Test")
test_that("Test for Wald Tests using  ", {
  wm1<- mix(Reaction ~ Days + (0+Days|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"))
  
  beta_1 <- waldTest(wm1, type="beta", coefs="Days", hypothesis=10.4073- 1.5458 *1.96)
  #wald test and T test are the same for individual betas 
  expect_equal(beta_1$p,.05,.01)
  
  beta_2<- waldTest(wm1, type="beta", coefs="(Intercept)", hypothesis=251.4051- 6.8246 *1.96)
  #wald test and T test are the same for individual betas 
  expect_equal(beta_2$p,.05,.01)
  
  #  Ensure that tests for Lambda run
  wm2<- mix(Reaction ~ Days + (1+Days|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"))
  waldTest(wm2,type="Lambda")
  waldTest(wm2,type="Lambda",coefs = c("Subject.Days.(Intercept)","Subject.Days" ))
})

context("Complex weighted three level model")
test_that("Complex weighted three level model", {
  skip_on_cran()
  sleepstudy2 <- sleepstudy
  sleepstudy2$Group <- 1
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("310", "309", "349", "335"), 2, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("330", "352", "337", "369"), 3, sleepstudy2$Group)
  sleepstudy2$Group <- ifelse(sleepstudy2$Subject %in% c("331", "332", "334", "372"), 4, sleepstudy2$Group)
  sleepstudy2$Group <- factor(sleepstudy2$Group)
  sleepstudy2$NewVar <- c(9.43457771946133, 12.7183652886321, 11.242490364576, 8.32049769489214, 21.4769136458094, 12.801691804759, 11.2944480319612, 23.2877967612343, 17.6554807720438, 15.6839576150421,
                          12.0415732799454, 12.6326929342708, 10.7497191412678, 8.67408880103278, 16.7872963953428, 14.6248001125594, 14.3553882674867, 11.8275049579491, 18.7525279773655, 5.57650843995549, 
                          19.7623973348143, 19.1858339464439, 16.9845555399241, 5.79214586016569, 13.6760775160985, 8.20870365747105, 10.469622698121, 9.5294364080553, 0.474272552642729, 9.90998530696626,
                          7.85598926755497, 13.6908324775899, 8.363578606421, 12.4585008913942, 20.8371544197279, 6.95572306644593, 4.49104017816579, 13.529348727477, 8.31255622061886, 8.78452997885101, 
                          17.216378291436, 13.6427335980228, 19.498592420836, 7.77906722601152, 15.1266474871889, 12.8030124499139, 6.62054433117513, 8.66974559938289, 11.9540398510289, 14.4271543736958,
                          10.0496588773239, 15.5177485124749, 9.12750122498998, 10.1801282848309, 14.2562515194416, 19.3249322019648, 10.7021090784807, 19.5336266451282, 7.58133671725254, 19.889529248361, 
                          5.32630146259616, 13.6044477039579, 4.98570236262749, 23.8801945918399, 14.447634009947, 15.6043025989833, 22.411260279495, 16.66695588071, 13.2723030644539, 17.6752877679585,
                          13.4050771535649, 16.8476633007641, 14.4754165647234, 13.2597805356783, 12.2720489082021, 16.3515751846435, 18.2589186962967, 13.5608781488381, 22.128982372665, 14.0340540905538,
                          2.08454544532827, 10.8863647522786, 21.0103204678449, 18.4595691535641, 11.9914499297614, 9.96791953740012, 20.0728999616675, 10.7237225063556, 18.6117287544885, 15.6290844263721,
                          13.584186380461, 12.8188973592142, 8.73588740878364, 22.1049445716224, 8.82190889109436, 16.7501795316447, 11.7894955061621, 8.67751729426691, 8.28560204932531, 25.3725616103908,
                          14.0620375346182, 14.9882805901214, 17.181272042065, 19.6581190406994, 11.6757977643058, 16.638959090019, 9.68436466184703, 23.1282879805077, 20.7690736926973, 14.3401954001056,
                          16.2747349705209, 8.61220295948687, 19.869302293393, 5.78478905807136, 6.59988762289215, 20.489844013367, 27.239468499143, 20.5797203990527, 19.4544929129352, 12.5320006698243, 
                          18.2664525619468, 16.6022500953422, 13.0216995414837, 12.4573456587016, 18.6922622366627, 7.37338081469625, 18.2016743123929, 12.9216600044986, 17.2543801572655, 20.1850819959469,
                          9.46342758905614, 16.1505649461445, 11.644647125078, 17.6992809329738, 13.5726089392145, 11.9006738563255, 9.11529530502925, 11.9473408244809, 13.7037965361481, 21.052192466607, 
                          8.94718178506257, 15.165337112518, 8.23261412105177, 12.2188157585207, 14.1455895553682, 4.94479286361234, 16.9113105275149, 12.7165809453375, 5.58774136619933, 8.31670629685991,
                          6.67728085523989, 11.7741431392146, 7.91590113917014, 9.13428547474213, 18.5689544076891, 18.1995195536489, 13.3510327364336, 21.4047919488063, 17.7821335119331, 9.19515401988151, 
                          10.7918240635605, 6.63952336871624, 15.3733302621749, 15.5262403960572, 13.4421915705244, 11.5980148647886, 17.3920393426159, 6.6883267354037, 10.6865661488425, 15.1104585129107,
                          9.99190979270527, 16.54240123493, 13.5687018760124, 13.5592897858269, 19.4019579232211, 15.0561754718535, 9.36036662501835, 4.27173388144395, 19.2040511963462, 6.36268714321505)
  
  ss2 <- sleepstudy2
  w1c <- w1 <- rep(1,180)
  w2c <- w2 <- rep(1,18)
  w3c <- w3 <- rep(1,4)
  # unbalanced (non-identical within a Subject), level-1 (obs level) weights
  w1c[1:4] <- w1[1:4] <- 2
  uR <- sleepstudy2[1:4,]
  sleepstudy2 <- rbind(sleepstudy2, uR)
  # level-2 weights
  w2c[1] <- w2[1] <- 2
  w1[ss2$Subject=="308"] <- 2*w1[ss2$Subject=="308"]
  sR <- subset(sleepstudy2, Subject == "308")
  sR$Subject <- "S308"
  sleepstudy2 <- rbind(sleepstudy2, sR)
  # level-3 weights
  w3c[1] <- w3[1] <- 2
  w2[c(1,7,12,13,16,17)] <- 2*w2[c(1,7,12,13,16,17)]
  w1[ss2$Group==1] <- 2*w1[ss2$Group==1]
  
  gR <- subset(sleepstudy2, Group %in% 1)
  gR$Subject <- paste0("G", gR$Subject)
  gR$Group <- "11"
  sleepstudy2 <- rbind(sleepstudy2, gR)
  
  
  ss2$w1 <- w1
  ss2$w2 <- rep(w2,each=10)
  ss2$w3 <- ifelse(ss2$Subject %in% c("308", "333", "350", "351", "370", "371"),2,1)
  
  ss2$n1 <- ss2$n2 <- ss2$n3 <- 1
  #lme reference with duplicated data
  lmr <- lmer(Reaction ~ Days + (1|Subject) + (1+Days|Group), data=sleepstudy2, REML=FALSE)

  #Does WeMix match to duplicated lmr with one+two random effects
  wmr <- mix(Reaction ~ Days + (1|Subject) + (1+Days|Group), data=ss2,weights=c("w1","w2","w3"))
  expect_equal(wmr$lnl, as.numeric(logLik(lmr)), tol=1e-3)
  expect_equal(coef(wmr), fixef(lmr), tol=1e-4)
  lmevars1 <- data.frame(summary(lmr)$varcor)
  vars <- lmevars1[is.na(lmevars1$var2),"vcov"]
  expect_equal(unname(wmr$vars),
               expected = unname(vars),
               tolerance = 1e-3)
  expect_equal(coef(wmr), fixef(lmr), tol=1e-4)
  expect_equal(unname(wmr$vars[length(wmr$vars)]), unname(lmr@devcomp$cmp["sigmaML"]^2), tol=1e-4)
  
  #lme reference with duplicated data
  lmr2 <- lmer(Reaction ~ Days + (1 + NewVar |Subject) + (1+Days|Group), 
               data=sleepstudy2, control=lmerControl(optimizer="bobyqa"), REML=FALSE)
  #wemix with weights
  wmr2 <- mix(Reaction ~ Days + (1+ NewVar |Subject) + (1+Days|Group), data=ss2,weights=c("w1","w2","w3"))
  #Does WeMix match to duplicated lmr with two correlated random effects at two levesl 
  expect_equal(wmr2$lnl, as.numeric(logLik(lmr2)), tol=1e-3)
  expect_equal(coef(wmr2), fixef(lmr2), tol=1e-4)
  
  lmevars2 <- data.frame(summary(lmr2)$varcor)
  vars <- lmevars2[is.na(lmevars2$var2),"vcov"]
  expect_equal(unname(wmr2$vars),
               expected = unname(vars),
               tolerance = 1e-3)
  
  expect_equal(coef(wmr2), fixef(lmr2), tol=1e-4)
  expect_equal(unname(wmr2$vars[length(wmr2$vars)]), unname(lmr2@devcomp$cmp["sigmaML"]^2), tol=1e-4)
}) 

if(!exists("edsurveyHome")) {
  if (Sys.info()[['sysname']] == "Windows") {
    edsurveyHome <- "C:/EdSurveyData/"
  } else {
    edsurveyHome <- "~/"
  }
}

context("PISA tests")
test_that("PISA tests", {
  skip_on_cran()
  require(EdSurvey)
  #read in data 

  options(timeout=60*60)
  downloadPISA(root=edsurveyHome, years=2012, cache=FALSE, verbose=FALSE)
  cntl <- readPISA(file.path(edsurveyHome, "PISA/2012"), countries = "USA", verbose=FALSE)
  om <- getAttributes(cntl, "omittedLevels")
  data <- getData(cntl,c("schoolid","pv1math","st29q03","sc14q02","st04q01",
                         "escs","w_fschwt","w_fstuwt"), 
                  omittedLevels = FALSE, addAttributes = FALSE)
  
  # Remove NA and omitted Levels
  om <- c("Invalid","N/A","Missing","Miss",NA,"(Missing)")
  for (i in 1:ncol(data)) {
    data <- data[!data[,i] %in% om,] 
  }

  #relevel factors for model 
  data$st29q03 <- relevel(data$st29q03,ref="Strongly agree")
  data$sc14q02 <- relevel(data$sc14q02,ref="Not at all")
  
  # Multivariate model with random intercept
  m1 <- mix(pv1math ~ st29q03 + sc14q02 + st04q01 + escs + (1|schoolid), data=data, weights=c("w_fstuwt", "w_fschwt"))
  m1bref <- matrix(c(486.8037, 7.777978, -11.1083,  5.699849, -19.25533, 5.455594, -41.5422,  6.864339,
                    -21.34052, 17.06059, -11.78236, 12.82083, -26.91253, 7.657342,  9.507693, 2.986006,
                     25.56825, 2.117479), ncol=2, byrow=TRUE)
  expect_equal(unname(summary(m1)$coef[,1:2]), m1bref, tol=1E-5)
  
  #test variance
  m1vref <- c(1413.81, 5264.799)
  expect_equal(unname(m1$vars), m1vref, tol=1E-5)
  
  #test lnl 
  expect_equal(m1$lnl ,-12789991.91, tol=1E-5)
  
  # var of var
  m1s <- summary(m1)
  # this is a regression test. See vignette for range of reasonable results.
  expect_equal(m1s$varsmat[,5], c(327.5015, 152.4944), tol=1e-4)
  
  data$pwt2 <- data$w_fschwt
  data$pwt1 <- data$w_fstuwt / data$w_fschwt
  # check conditional weights
  m1c <- mix(pv1math ~ st29q03 + sc14q02 + st04q01 + escs + (1|schoolid), data=data, weights=c("pwt1", "pwt2"), cWeight=TRUE)
  # call should disagree, so remove that
  m1$call <- NULL
  m1c$call <- NULL
  # the cConstructor disagrees
  m1$lnlf <- NULL
  m1c$lnlf <- NULL
  expect_equal(m1, m1c, tol=1e-5)
  #test complicated model
  m2 <- mix(pv1math ~ st29q03 + sc14q02 + st04q01 + escs + (1|schoolid) + (0+escs|schoolid), data=data, weights=c("w_fstuwt", "w_fschwt"))
  expect_equal(m2$lnl ,-12741522.65, tol=1E-5)
    
  m2bref <- matrix(c(483.8881, 7.405499, -10.62803, 5.490497, -17.30314, 5.372283, -38.70806, 6.768197, -20.06462, 15.58282, -12.59165, 12.66878, -39.94083, 7.239317, 10.29371, 3.028518, 28.0161, 2.563144), ncol=2, byrow=TRUE)
  expect_equal(unname(summary(m2)$coef[,1:2]), m2bref,tol=1E-5)
  
  #test  SE of variance based on mixed
  se_var2 <- summary(m2)$vars[,2] # regression test, ideal unclear
  expect_equal(se_var2 ,c(287.89, 67.99, 137.70), tol=3e-5)
  
  #test variance
  m2vref <- c(1354.71,370.33,4967.36)
  expect_equal(unname(m2$vars), m2vref, tol=2E-5)

  suppressWarnings(m3 <- mix(pv1math ~ st29q03 + sc14q02 + st04q01 + escs+ (1+escs|schoolid), data=data, weights=c("w_fstuwt", "w_fschwt")))
  expect_equal(m3$lnl, -12741285.44518178, tol=1E-5)

  m3bref <- matrix(c(483.5955, 7.494755,
                     -10.62607, 5.486497, -17.26411, 5.371982, -38.6784, 6.757546,
                     -20.55204, 15.84887, -8.744658, 12.42958, -33.18773, 6.897756,
                     10.31113, 3.02734,
                     27.63378, 2.435956), ncol=2, byrow=TRUE)
  expect_equal(unname(summary(m3)$coef[,1:2]), m3bref, tol=1E-5)
})


context("Model Matrix has a hard time with")
test_that("Model Matrix has a hard time with", {
  skip_on_cran()
  require(EdSurvey)
  sdf <- readNAEP(system.file("extdata/data", "M36NT2PM.dat", package = "NAEPprimer"))
  gg <- getData(varnames=c("composite", "dsex", "b017451", "scrpsu", "origwt", "smsrswt"), data=sdf, returnJKreplicates=FALSE)
  gg2 <- gg[gg$origwt > 0 & gg$smsrswt > 0,]
  suppressMessages(m4 <- mix(mrpcm2 ~ dsex + b017451 + (1|scrpsu), data=gg2, weights=c("origwt", "smsrswt")))
  expect_equal(m4$lnl, -81882.3634148408, tol=1e-5)
  expect_equal(m4$coef, c(`(Intercept)` = 270.571817625202, dsexFemale = -2.14907551600309, `b017451Once every few weeks` = 3.84433128452533, `b017451About once a week` = 9.19954009166631, `b0174512 or 3 times a week` = 12.8701977366809, `b017451Every day` = 6.29635843233831))
  expect_equal(m4$SE, c(`(Intercept)` = 1.14859092700564, dsexFemale = 0.636207386283148, `b017451Once every few weeks` = 1.05190060207538, `b017451About once a week` = 0.986062663493941, `b0174512 or 3 times a week` = 1.00366545196335, `b017451Every day` = 1.10858750669702))
  m4varDF <- structure(list(grp = c("scrpsu", "Residual"),
                            var1 = c("(Intercept)", NA),
                            var2 = c(NA_character_, NA_character_),
                            vcov = c(300.695102536901, 969.501918646696),
                            ngrp = c(672, 16331),
                            level = c(2, 1),
                            SEvcov = c(25.0004604746421, 16.2937862202021),
                            fullGroup = c("scrpsu.(Intercept)", "Residual")),
                       row.names = c(NA, -2L),
                       class = "data.frame")
  expect_equal(m4$varDF, m4varDF, tol=1e-5)
})



context("examples run")
test_that("examples run", {
  skip_on_cran()
  ss1 <- sleepstudy

  # Create weights
  ss1$W1 <- ifelse(ss1$Subject %in% c(308, 309, 310), 2, 1)
  ss1$W2 <- 1

  # Run random-intercept 2-level model 
  two_level <- mix(Reaction ~ Days + (1|Subject), data=ss1, weights=c("W1", "W2"))
  expect_is(two_level, "WeMixResults")

  #Run random-intercept 2-level model with group-mean centering
  grp_centered <- mix(Reaction ~ Days + (1|Subject), data=ss1,
                      weights = c("W1", "W2"),
                      center_group = list("Subject" = ~Days))
  expect_is(grp_centered, "WeMixResults")

  #Run three level model with random slope and intercept. 
  #add group variables for 3 level model 
  ss1$Group <- 3
  ss1$Group <- ifelse(as.numeric(ss1$Subject) %% 10 < 7, 2, ss1$Group)
  ss1$Group <- ifelse(as.numeric(ss1$Subject) %% 10 < 4, 1, ss1$Group)
  # level-3 weights
  ss1$W3 <- ifelse(ss1$Group == 2, 2, 1)

  three_level <- mix(Reaction ~ Days + (1|Subject) + (1+Days|Group), data=ss1, 
                     weights=c("W1", "W2", "W3"))
  expect_is(three_level, "WeMixResults")

  # Conditional Weights
  # use vignette example
  library(EdSurvey)

  #read in data 
  cntl <- readPISA(file.path(edsurveyHome,"PISA/2012"), countries="USA", verbose=FALSE)
  data <- getData(cntl,c("schoolid","pv1math","st29q03","sc14q02","st04q01",
                         "escs","w_fschwt","w_fstuwt"), 
                  omittedLevels=FALSE, addAttributes=FALSE)

  # Remove NA and omitted Levels
  om <- c("Invalid", "N/A", "Missing", "Miss", NA, "(Missing)")
  for (i in 1:ncol(data)) {
    data <- data[!data[,i] %in% om,]
  }

  #relevel factors for model 
  data$st29q03 <- relevel(data$st29q03, ref="Strongly agree")
  data$sc14q02 <- relevel(data$sc14q02, ref="Not at all")

  # run with unconditional weights
  m1u <-  mix(pv1math ~ st29q03 + sc14q02 +st04q01+escs+ (1|schoolid), data=data, 
              weights=c("w_fstuwt", "w_fschwt"))
  expect_is(m1u, "WeMixResults")

  # conditional weights
  data$pwt2 <- data$w_fschwt
  data$pwt1 <- data$w_fstuwt / data$w_fschwt

  # run with conditional weights
  m1c <-  mix(pv1math ~ st29q03 + sc14q02 +st04q01+escs+ (1|schoolid), data=data, 
              weights=c("pwt1", "pwt2"), cWeights=TRUE)
  expect_is(m1c, "WeMixResults")
  # the results are, up to rounding, the same in m1u and m1c, only the calls are different

})

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
  downloadECLS_K(root=edsurveyHome, years=2011)
  ee <- readECLS_K2011(paste0(edsurveyHome, "ECLS_K/2011/"), verbose=FALSE)
  gg <- getData(c("x2rscalk5", "childid", "s2_id", "w1_2p0", "x3sumsh", "p1chldbk","p2freerd"), data=ee, omittedLevels=FALSE, returnJKreplicates=FALSE)
  gg$frpl <- ifelse(gg$p2freerd %in% c("1: FREE LUNCH", "2: REDUCED PRICE LUNCH"), 1, 0)
  gg$w1 <- gg$w1_2p0
  gg$w2 <- 1
  gg$n <- ave(gg$s2_id,gg$s2_id, FUN=length)
  gg2 <- gg[!is.na(gg$x2rscalk5) & gg$w1>0 & !is.na(gg$p1chldbk) & gg$n > 15 & gg$s2_id < 1500,]
  m3 <- mix(x2rscalk5 ~ p1chldbk + frpl + (1+frpl|s2_id), data=gg2, weights=c("w1", "w2"), verbose=FALSE)
  # regression tests
  expect_equal(m3$lnl, -1513119.73817376, tol=1e-5)
  expect_equal(m3$coef, c(`(Intercept)` = 67.5434743459653, p1chldbk = 0.0287511027423215, frpl = -4.62934048945893), tol=1e-5)
  expect_equal(m3$SE, c(`(Intercept)` = 1.07449754776874, p1chldbk = 0.00668096358477409, frpl = 1.12769131710247), tol=1e-5)
  varDF0 <- structure(list(grp = c("s2_id", "s2_id", "s2_id", "Residual"), 
                           var1 = c("(Intercept)", "frpl", "(Intercept)", NA),
                           var2 = c(NA, NA, "frpl", NA),
                           vcov = c(82.1115665535057, 103.939072336388, -57.1975297909512, 151.622543466079),
                           ngrp = c(97, 97, 97, 1520),
                           level = c(2, 2, 2, 1),
                           SEvcov =  c(17.082223467814, 20.4018480839976, 15.328078288835, 10.0505684893662),
                           fullGroup = c("s2_id.(Intercept)", "s2_id.frpl", "s2_id.(Intercept)", "Residual")),
                       row.names = c(NA, -4L), class = "data.frame")
  expect_equal(m3$varDF, varDF0, tol=1e-5)
})

context("TIMSS tests")
test_that("TIMSS tests", {
  skip_on_cran()
  require(EdSurvey)

  # original version by Christian Kjeldsen
  downloadTIMSS(root=edsurveyHome, years=2015, cache=FALSE, verbose=FALSE)
  dnk15 <- readTIMSS(file.path(edsurveyHome,"/TIMSS/2015"), countries="dnk", gradeLvl=4, verbose=FALSE)
  dnk15dat <- getData(data=dnk15, varnames=c("atbg01", "asbgsb", "mmat", "asbghrl", "matwgt", "idschool","schwgt"))
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
  mm1w <- waldTest(fittedModel=mm1, type="beta", coefs="atbg01")
 
  expect_equal(mm1s$coef[2,3]^2, mm1w$Wald)

})
