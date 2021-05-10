## ----code options, echo = FALSE----------------------------------------------------
options(width = 85)

## ----source package, eval=FALSE----------------------------------------------------
#  install.packages("WeMix")

## ----load package, results="hide", message=FALSE,eval=FALSE------------------------
#  library(WeMix)

## ----wemix, eval=FALSE-------------------------------------------------------------
#  # model with one random effect
#  m1 <-  mix(pv1math ~ st29q03 + sc14q02 +st04q01+escs+ (1|schoolid), data=data,
#      weights=c("w_fstuwt", "w_fschwt"))
#  
#  # model with two random effects assuming zero correlation between the two
#  m2 <- mix(pv1math ~ st29q03 + sc14q02 +st04q01+escs+ (1|schoolid)+ (0+escs|schoolid),
#      data=data, weights=c("w_fstuwt", "w_fschwt"))
#  
#  #Wald tests for beta and Lambda parameters
#  waldTest(m2,type="beta")
#  waldTest(m2,type="Lambda")

## ----eval=FALSE--------------------------------------------------------------------
#  library(lme4)
#  library(WeMix)
#  ss1 <- sleepstudy
#  doubles <- c(308, 309, 310) # subject with double obs
#  
#  # Create weights
#  ss1$W1 <- ifelse(ss1$Subject %in% doubles, 2, 1)
#  ss1$W2 <- 1
#  
#  # Create binary outcome variable called "over300"
#  ss1$over300 <- ifelse(sleepstudy$Reaction<300,0,1)
#  
#  # Run mixed model with random intercept and fixed slope
#  bi_1 <- mix(over300~ Days + (1|Subject),data=ss1,
#              family=binomial(link="logit"),verbose=FALSE,
#              weights=c("W1", "W2"),nQuad=13)

## ----Centering, eval=FALSE---------------------------------------------------------
#  library(lme4)  #to use example sleep study data
#  
#  #create dummy weights
#  sleepstudy$weight1L1 <- 1
#  sleepstudy$weight1L2 <- 1
#  
#  # Group mean centering of the variable Days within group Subject
#  group_center <- mix(Reaction ~ Days + (1|Subject), data=sleepstudy,
#                      center_group=list("Subject"= ~Days),
#                      weights=c("weight1L1", "weight1L2"))
#  
#  # Grand mean centering of the variable Days
#  grand_center <- mix(Reaction ~ Days + (1|Subject), data=sleepstudy,
#                      center_grand=~Days,weights=c("weight1L1", "weight1L2"))

## ----Stata, eval=FALSE-------------------------------------------------------------
#  import delimited "PISA2012_USA.csv"
#  
#  generate intercept = 1
#  eq intercept: intercept
#  eq slope: escs
#  tabulate st29q03, generate (st29q03d)
#  tabulate sc14q02, generate (sc14q02d)
#  tabulate st04q01, generate (st04q01d)
#  
#  //Random intercept model
#  gllamm pv1math st29q03d1 st29q03d2 st29q03d4 sc14q02d1 sc14q02d3 sc14q02d4 st04q01d2
#  escs, i(schoolid) pweight(pwt) l(identity) f(gau) nip(27) nrf(1) eqs(intercept)
#  adapt nocorrel
#  
#  
#  //Random slope and intercept model
#  gllamm pv1math st29q03d1 st29q03d2 st29q03d4 sc14q02d1 sc14q02d3 sc14q02d4 st04q01d2
#  escs, i(schoolid) pweight(pwt) l(identity) f(gau) nip(13) nrf(2) eqs(intercept slope)
#  adapt nocorrel

## ----Statamixed, eval=FALSE--------------------------------------------------------
#  import delimited "PISA2012_USA.csv"
#  tabulate st29q03, generate (st29q03d)
#  tabulate sc14q02, generate (sc14q02d)
#  tabulate st04q01, generate (st04q01d)
#  
#  //Random intercept model
#  mixed pv1math st29q03d1 st29q03d2 st29q03d4 sc14q02d1 sc14q02d3 sc14q02d4 st04q01d2
#  escs [pw = pwt1] || schoolid: , pweight (pwt2)
#  
#  //Random slope and intercept model
#  mixed pv1math st29q03d1 st29q03d2 st29q03d4 sc14q02d1 sc14q02d3 sc14q02d4 st04
#  q01d2 escs [pw = pwt1] || schoolid: escs, pweight (pwt2)

## ----sas, eval=FALSE---------------------------------------------------------------
#  PROC IMIPORT DATAFILE="PISA2012_USA.csv"
#       OUT=pisa_data DBMS=csv REPLACE;
#  RUN;
#  
#  PROC GLIMMIX DATA=pisa_data METHOD=quadrature(qpoints=27) EMPIRICAL=classical NOREML;
#    NLOPTIONS GCONV=1E-10 TECHNIQUE=TRUREG;
#    CLASS  sc14q02(REF='Not at all') st04q01(REF='Female') st29q03(REF='Strongly agree');
#    MODEL pv1math = escs sc14q02 st04q01 st29q03/ OBSWEIGHT=pwt1 SOLUTION;
#    RANDOM INT/subject=schoolid WEIGHT=pwt2;
#  RUN;
#  
#  
#  PROC GLIMMIX DATA=pisa_data METHOD=quadrature(qpoints=13) EMPIRICAL=classical NOREML;
#    NLOPTIONS GCONV=1E-10 TECHNIQUE=TRUREG;
#    CLASS  sc14q02(REF='Not at all') st04q01(REF='Female') st29q03(REF='Strongly agree');
#    MODEL pv1math = escs sc14q02 st04q01 st29q03/ OBSWEIGHT=pwt1 SOLUTION;
#    RANDOM intercept escs/subject=schoolid WEIGHT=pwt2;
#  RUN;

## ----M_plus_one, eval=FALSE--------------------------------------------------------
#  DATA: FILE= pisa_2012_with_dummies.csv;
#  
#  VARIABLE: NAMES= id schoolid pv1math escs pwt1 pwt2 s29q03d1 s29q03d2
#  s29q03d4 int c14q02d1 c14q02d2 c14q02d4 s04q01d2;
#  CLUSTER = schoolid;
#  WITHIN = escs s29q03d1 s29q03d2 s29q03d4 c14q02d1 c14q02d2 c14q02d4 s04q01d2;
#  
#  USEVARIABLES= schoolid pv1math escs  s29q03d1 s29q03d2
#  s29q03d4 c14q02d1 c14q02d2 c14q02d4 s04q01d2 pwt1 pwt2;
#  
#  WEIGHT IS pwt1;
#  BWEIGHT = pwt2;
#  wtscale=unscaled;
#  bwtscale=unscaled;
#  
#  ANALYSIS: TYPE = TWOLEVEL;
#  
#  MODEL:
#  %WITHIN%
#  pv1math ON escs s29q03d1 s29q03d2 s29q03d4 c14q02d1 c14q02d2 c14q02d4 s04q01d2;
#  %BETWEEN%
#  pv1math;

## ----M_plus, eval=FALSE------------------------------------------------------------
#  DATA: FILE=pisa_2012_with_dummies.csv;
#  
#  VARIABLE: NAMES= id schoolid pv1math escs pwt1 pwt2 s29q03d1 s29q03d2
#  s29q03d4 int c14q02d1 c14q02d2 c14q02d4 s04q01d2 escs2;
#  CLUSTER = schoolid;
#  WITHIN = escs s29q03d1 s29q03d2 s29q03d4 c14q02d1 c14q02d2 c14q02d4 s04q01d2;
#  
#  USEVARIABLES= schoolid pv1math escs  s29q03d1 s29q03d2
#  s29q03d4 c14q02d1 c14q02d2 c14q02d4 s04q01d2 pwt1 pwt2;
#  
#  WEIGHT IS pwt1;
#  BWEIGHT = pwt2;
#  wtscale=unscaled;
#  bwtscale=unscaled;
#  
#  ANALYSIS: TYPE = TWOLEVEL RANDOM;
#  
#  MODEL:
#  %WITHIN%
#  pv1math on s29q03d1 s29q03d2 s29q03d4 c14q02d1 c14q02d2 c14q02d4 s04q01d2;
#  slope | pv1math ON escs ;
#  %BETWEEN%
#  pv1math with slope @0

## ----data_prep, eval=FALSE---------------------------------------------------------
#  library(EdSurvey)
#  
#  #read in data
#  downloadPISA("~/", year=2012)
#  cntl <- readPISA("~/PISA/2012", countries="USA")
#  data <- getData(cntl,c("schoolid", "pv1math", "st29q03", "sc14q02", "st04q01",
#                         "escs", "w_fschwt", "w_fstuwt"),
#                  omittedLevels=FALSE, addAttributes=FALSE)
#  
#  # conditional weights
#  data$pwt2 <- data$w_fschwt
#  data$pwt1 <- data$w_fstuwt / data$w_fschwt
#  
#  # Remove NA and omitted Levels
#  om <- c("Invalid","N/A","Missing","Miss",NA,"(Missing)")
#  for (i in 1:ncol(data)) {
#    data <- data[!data[,i] %in% om,]
#  }
#  
#  # relevel factors for model
#  data$st29q03 <- relevel(data$st29q03, ref="Strongly agree")
#  data$sc14q02 <- relevel(data$sc14q02, ref="Not at all")
#  
#  write.csv(data, file="PISA2012_USA.csv")

