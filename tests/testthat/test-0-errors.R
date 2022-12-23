require(testthat)

options(width = 500)
options(useFancyQuotes = FALSE)

sleepstudyU <- sleepstudy
sleepstudyU$weight1L1 <- 1
sleepstudyU$weight1L2 <- 1
sleepstudyU$weight1L3 <- 1

context("Errors correctly report")
test_that("Errors correctly report", {
  # error: weight not on data frame
  # should mention all the weights that are not on the data frame
  expect_error(
    mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1ZZ", "weight1L2ZZ"))
    , regexp="(?=.*weight1L1ZZ)(?=.*weight1L2ZZ)", perl=TRUE)
  # error: too many weights, should 
  expect_error(
    mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2", "weight1L3"))
    , "(?=.*3)(?=.*expecting 2)", perl=TRUE)
  # different message for only one weight because other tests can't proceed with only one weight
  expect_error(
    mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1"))
    , "(?=.*1)(?=.*expecting 2)", perl=TRUE)
  # unimplemented families
  expect_error(
    mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), family=Gamma())
    , "Gamma", perl=TRUE)
  expect_error(
    mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), family=inverse.gaussian())
    , "inverse", perl=TRUE)
  expect_error(
    mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), family=quasibinomial())
    , "quasibinomial", perl=TRUE)
  # a center_group must be a list named after groups
  # legit call:
  # mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), center_group=list("Subject"= ~Days))
  # not a list, just a formula
  expect_error(
    mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), center_group=Subject ~ novar)
    , "must be a list", perl=TRUE)
  # list but non-group-level name
  expect_error(
    mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), center_group=list("Days"= ~Days))
    , "Could not find level of \"Days\"", perl=TRUE)
  # non-0/1 outcomes; note the binomial will trip over 0 <= y <= 1, so needs to be in that range
  expect_error(
    suppressMessages(suppressWarnings(
      mix(I(Reaction/max(Reaction)) ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2"), family=binomial())
    ))
    , "outcomes must be 0 or 1", perl=TRUE)
  #   
  set.seed(2)
  sleepstudyU$weight1L2Bad <- sample(c(rep(1, nrow(sleepstudyU)-3), runif(3)))
  expect_error(
    capture.output(suppressMessages(
      mix(Reaction ~ Days + (1|Subject), data=sleepstudyU, weights=c("weight1L1", "weight1L2Bad"))
    ))
    , "level-2 weights vary", perl=TRUE)

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
  
  expect_error(
    capture.output(suppressMessages(
      wm0 <- mix(Reaction ~ Days + (1|Subject/Group), data=ss2, weights=c("w1", "w2", "w3"))
    ))
    , "nested model"
    )


})
