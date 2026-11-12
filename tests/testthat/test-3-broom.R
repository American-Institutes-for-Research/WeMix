context("broom helpers work")
require(lme4)
library(broom)
require(dplyr)

slp <- sleepstudy
slp$w1 <- 1
slp$w2 <- 1

m1 <- WeMix::mix(Reaction ~ Days + (Days | Subject), slp, weights=c("w1", "w2"))
m2 <- lme4::lmer(Reaction ~ Days + (Days | Subject), slp, weights=slp$w1, REML=FALSE) # WeMix uses REML


test_that("terms, model.frame, model.matrix work", {
  # should be almost equivalent
  expect_equal(terms(m1), terms(m2, fixed.only=FALSE))
  # ignores attributes:
  expect_equivalent(model.frame(m1) %>% select(-"w1", -"w2"), model.frame(m2) %>% select(-"(weights)"))
  expect_equal(model.matrix(m1), model.matrix(m2))
})

test_that("tidy works", {
  skip_if_not_installed("broom.mixed")
  library(broom.mixed)
  # check the first three columns (the description)
  expect_equal(tidy(m1)[,1:4], tidy(m2)[,1:4], tolerance=1e-5 )
  # the standard errors from WeMix are robust and so cannot be compared to lme4::lmer
})

m3 <- WeMix::mix(I(Reaction>290) ~ Days + (Days | Subject), slp, weights=c("w1", "w2"), family=binomial(), nQuad=1)
m4 <- lme4::glmer(I(Reaction>290) ~ Days + (Days | Subject), slp, weights=slp$w1, family=binomial())
test_that("non-linear models: terms, model.frame, model.matrix, tidy", {
  expect_equal(terms(m3), terms(m4, fixed.only=FALSE))
  expect_equivalent(model.frame(m3) %>% select(-"w1", -"w2"), model.frame(m4) %>% select(-"(weights)"))
  expect_equal(model.matrix(m3), model.matrix(m4))
  expect_equal(tidy(m3)[,1:4],tidy(m4)[,1:4], tolerance=1e-2)
})