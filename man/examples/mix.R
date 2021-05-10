\dontrun{
library(lme4)

data(sleepstudy)
ss1 <- sleepstudy

# Create weights
ss1$W1 <- ifelse(ss1$Subject %in% c(308, 309, 310), 2, 1)
ss1$W2 <- 1

# Run random-intercept 2-level model 
two_level <- mix(Reaction ~ Days + (1|Subject), data=ss1, weights=c("W1", "W2"))

#Run random-intercept 2-level model with group-mean centering
grp_centered <- mix(Reaction ~ Days + (1|Subject), data=ss1,
                    weights = c("W1", "W2"),
                    center_group = list("Subject" = ~Days))

#Run three level model with random slope and intercept. 
#add group variables for 3 level model 
ss1$Group <- 3
ss1$Group <- ifelse(as.numeric(ss1$Subject) %% 10 < 7, 2, ss1$Group)
ss1$Group <- ifelse(as.numeric(ss1$Subject) %% 10 < 4, 1, ss1$Group)
# level-3 weights
ss1$W3 <- ifelse(ss1$Group == 2, 2, 1)

three_level <- mix(Reaction ~ Days + (1|Subject) + (1+Days|Group), data=ss1, 
                   weights=c("W1", "W2", "W3"))

# Conditional Weights
# use vignette example
library(EdSurvey)

#read in data 
downloadPISA("~/", year=2012)
cntl <- readPISA("~/PISA/2012", countries="USA")
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
m1u <- mix(pv1math ~ st29q03 + sc14q02 +st04q01+escs+ (1|schoolid), data=data, 
           weights=c("w_fstuwt", "w_fschwt"))
summary(m1u)

# conditional weights
data$pwt2 <- data$w_fschwt
data$pwt1 <- data$w_fstuwt / data$w_fschwt

# run with conditional weights
m1c <- mix(pv1math ~ st29q03 + sc14q02 +st04q01+escs+ (1|schoolid), data=data, 
            weights=c("pwt1", "pwt2"), cWeights=TRUE)
summary(m1c)
# the results are, up to rounding, the same in m1u and m1c, only the calls are different

}
