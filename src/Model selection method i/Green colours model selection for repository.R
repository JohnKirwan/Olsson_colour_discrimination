rm(list = ls()) #clear the workspace of anything before starting








## Loading packages required to run the analysis and the plotting
library(lattice)
library(graphics)
library(psyphy)
library(lmerTest)
library(MASS)
library(glm2)
library(readr)
#Load data
load("(replace with folder)/Green colours data.RData")
#or
# data <- read.csv('green colours.csv', header = TRUE, sep = ";", row.names = NULL) # change the "Working Directory" to the folder containing your data file. Change the name xxxxx.csv to your filename
data$background <- factor(data$background)
data$ind <- factor(data$ind)
data$batch <- as.factor(data$batch)


#Green colour
CheckRank <- function(mod){
  cf <- coef(mod)
  if(length(cf)>mod$rank){
    cbind(names(cf)[is.na(cf)])
  }else{'FULL RANK'}
}#a function to check if models are rank deficient and return missing coef

NumCoef <- function(mod, coef = 'ind'){
  cf <- coef(mod)#coefficients
  fit <- c(names(cf[!is.na(cf)]))#which ones did fit
  return(length(grep(coef, fit)))#how many are there
}
#Model testing, null model and full model, and models removing one factor
mod.0 <- glm.WH(cbind(corr,incorr)~1, data = data, lambda.init = 0.1, interval = c(0, 1), NumAlt = 2) # Null model, data is explained only by variance. lambda.init is your estimate of the upper asymptote and NumAlt specifies the number of alternatives e.g. 2 alternative forced choice
CheckRank(mod.0)
#[1] "FULL RANK"
NumCoef(mod.0)
#[1] 0

# full model with 3-way interaction
mod.1 <- glm.WH(cbind(corr,incorr)~Colour.difference*background*sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=1000) # Model where data is affected by the stimulus level
CheckRank(mod.1)
                 
#[1,] "indD3"                   
#[2,] "indD4"                   
#[3,] "batchb"                  
#[4,] "batchc"                  
#[5,] "batchd"                  
#[6,] "backgroundorange:sexmale"


# compare null and full model
anova(mod.0,mod.1,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ 1
# 2: cbind(corr, incorr) ~ Colour.difference * background * sex + 
#  ind + batch
#Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
#1        90     804.06                          
#2        71     117.72 19   686.34 < 2.2e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

AIC(mod.0,mod.1)
#df       AIC
#mod.0  1 1090.6942
#mod.1 20  442.3526


### removing 3-way interaction
mod.2 <- glm.WH(cbind(corr,incorr)~Colour.difference*background+Colour.difference*sex+background*sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.2)
#[1,] "indD3"                   
#[2,] "indD4"                   
#[3,] "batchb"                  
#[4,] "batchc"                  
#[5,] "batchd"                  
#[6,] "backgroundorange:sexmale"

#comparison with previous prefered model
anova(mod.1,mod.2,test="Chisq")

#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background * sex + 
#  ind + batch
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + background * sex + ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        71     117.72                     
#2        72     120.09 -1  -2.3746   0.1233

AIC(mod.1,mod.2)
#df      AIC
#mod.1 20 442.3526
#mod.2 19 442.7272



#removing interaction background:sex
mod.3a <- glm.WH(cbind(corr,incorr)~Colour.difference*background+Colour.difference*sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.3a)
#[1,] "indD3" 
#[2,] "indD4" 
#[3,] "batchb"
#[4,] "batchc"
#[5,] "batchd"

#comparison with previous prefered model
anova(mod.2,mod.3a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + ind + batch
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        72     120.09                     
#2        72     120.09  0        0   
AIC(mod.3a,mod.2)
#df      AIC
#mod.3a 19 442.7272
#mod.2 19 442.7272
  
#removing interaction colour difference:background
mod.4a <- glm.WH(cbind(corr,incorr)~Colour.difference+background+Colour.difference*sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level

CheckRank(mod.4a)
#[1,] "indD3" 
#[2,] "indD4" 
#[3,] "batchb"
#[4,] "batchc"
#[5,] "batchd

anova(mod.4a,mod.3a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + ind + batch
#Model 2: cbind(corr, incorr) ~ Colour.difference + background + Colour.difference * 
#  sex + background * sex + ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        72     120.09                     
#2        73     120.95 -1 -0.85909    0.354
AIC(mod.4a,mod.3a)
#df      AIC
#mod.5a 18 441.5863
#mod.4a 19 442.7272



#removing interaction colour difference:sex
mod.5a <- glm.WH(cbind(corr,incorr)~Colour.difference+background+sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.5a)
#[1,] "indD3" 
#[2,] "indD4" 
#[3,] "batchb"
#[4,] "batchc"
#[5,] "batchd"

anova(mod.5a,mod.4a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + background + sex + 
#  ind + batch
#Model 2: cbind(corr, incorr) ~ Colour.difference + background + Colour.difference * 
#  sex + ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        74     121.09                     
#2        73     120.95  1  0.14278   0.7055

AIC(mod.5a,mod.4a)
#df      AIC
#mod.5a 17 439.7290
#mod.4a 18 441.5863


#removing batch
mod.6a <- glm.WH(cbind(corr,incorr)~Colour.difference+background+sex+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.6a)
#[1,] "indD3"
#[2,] "indD4"

anova(mod.6a,mod.5a,test="Chisq")
#Analysis of Deviance Table

# 1: cbind(corr, incorr) ~ Colour.difference + background + sex + 
#  ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + background + sex + 
#  ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        74     121.09                     
#2        74     121.09  0        0      
AIC(mod.6a,mod.5a)
#df     AIC
#mod.6a 17 439.729
#mod.5a 17 439.729

#removing ind
mod.7a <- glm.WH(cbind(corr,incorr)~Colour.difference+background+sex, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.7a)
#[1] "FULL RANK"

anova(mod.7a,mod.6a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + background + sex
#Model 2: cbind(corr, incorr) ~ Colour.difference + background + sex + 
#  ind
#Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
#1        87     193.30                          
#2        74     121.09 13   72.206 3.142e-10 ***
# ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

AIC(mod.7a, mod.6a)
#df      AIC
#mod.7a  4 485.9355
#mod.6a 17 439.7290

#removing sex
mod.7b <- glm.WH(cbind(corr,incorr)~Colour.difference+background+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.7b)
#[1,] "indD4"

anova(mod.7b,mod.6a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + background + ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + background + sex + 
#  ind
#Resid. Df Resid. Dev Df   Deviance Pr(>Chi)
#1        74     121.09                       
#2        74     121.09  0 8.5265e-14 

AIC(mod.7b,mod.6a)
#df     AIC
#mod.7b 17 439.729
#mod.6a 17 439.729

#removing background
mod.8a <- glm.WH(cbind(corr,incorr)~Colour.difference+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.8a)
#[1] "FULL RANK"
anova(mod.7b,mod.8a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + background + ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + ind
#Resid. Df Resid. Dev Df   Deviance Pr(>Chi)
#1        74     121.09                       
#2        74     121.09  0 5.6843e-14 

AIC(mod.8a,mod.7b)
#df     AIC
#mod.8a 17 439.729
#mod.7b 17 439.729

#removing ind
mod.9a <- glm.WH(cbind(corr,incorr)~Colour.difference, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.9a)
#[1] "FULL RANK"
anova(mod.9a,mod.8a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference
#Model 2: cbind(corr, incorr) ~ Colour.difference + ind
#Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
#1        89     227.24                          
#2        74     121.09 15   106.14 8.837e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

AIC(mod.9a,mod.8a)
#df      AIC
#mod.9a  2 515.8722
#mod.8a 17 439.7290

#removing colour diff.
mod.9b <- glm.WH(cbind(corr,incorr)~ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.9b)
#[1] "FULL RANK"
anova(mod.9b,mod.8a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + ind
# Df Resid. Dev Df Deviance  Pr(>Chi)    
#1        75     712.55                          
#2        74     121.09  1   591.46 < 2.2e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

AIC(mod.9b,mod.8a)
#df      AIC
#mod.9b 16 1029.191
#mod.8a 17  439.729

summary(mod.8a)
#Call:
#  glm(formula = formula, family = binomial(probit.lambda(NumAlt, 
#                                                         lam.c)), data = data, maxit = 100)
#
#Deviance Residuals: 
#  Min         1Q     Median         3Q        Max  
#-2.365566  -0.775267   0.001406   0.712578   2.511404  
#
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)        -2.2783     0.5042  -4.519 6.21e-06 ***
#  Colour.difference   4.1467     0.6513   6.367 1.93e-10 ***
#  indaB3             -1.3540     0.5606  -2.415 0.015718 *  
#  indaB4             -0.3866     0.5076  -0.762 0.446259    
#indaB5             -0.4293     0.5063  -0.848 0.396540    
#indbB22            -1.2040     0.5532  -2.176 0.029529 *  
#  indbB32            -1.9958     0.6064  -3.291 0.000997 ***
#  indbB42            -1.2695     0.5567  -2.280 0.022584 *  
#  indbB62            -1.0579     0.5455  -1.939 0.052457 .  
#indC1              -4.8498     0.9642  -5.030 4.91e-07 ***
#  indC2              -1.1198     0.5163  -2.169 0.030086 *  
#  indC3              -1.4399     0.5448  -2.643 0.008213 ** 
#  indC4               0.1441     0.4735   0.304 0.760874    
#indD1              -1.2797     0.5291  -2.418 0.015588 *  
#  indD2              -2.8411     1.3181  -2.155 0.031128 *  
#  indD3              -4.9984     0.9638  -5.186 2.14e-07 ***
#  indD4              -4.5181     0.9769  -4.625 3.75e-06 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 804.06  on 91  degrees of freedom
#Residual deviance: 121.09  on 74  degrees of freedom
#AIC: 439.73

#Number of Fisher Scoring iterations: 100

#lambda	 0.0651




##### Forward Selection below
#adding background
mod10a <- glm.WH(cbind(corr,incorr)~Colour.difference+background+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod10a)
#[1,] "indD4"
anova(mod.8a,mod10a,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + background + ind
#Resid. Df Resid. Dev Df    Deviance Pr(>Chi)
#1        74     121.09                        
#2        74     121.09  0 -5.6843e-14 

#adding sex
moda <- glm.WH(cbind(corr,incorr)~Colour.difference+sex+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod10a)
#[1,] "indD4"
anova(mod.8a,moda,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + sex + ind
#Resid. Df Resid. Dev Df    Deviance Pr(>Chi)
#1        74     121.09                        
#2        74     121.09  0 -5.6843e-14     

#adding batch
modb <- glm.WH(cbind(corr,incorr)~Colour.difference+batch+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(modb)
#[1,] "indbB62"
#[2,] "indC4"  
#[3,] "indD4"  
anova(mod.8a,modb,test="Chisq")

#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + batch + ind
#Resid. Df Resid. Dev Df    Deviance Pr(>Chi)
#1        74     121.09                        
#2        74     121.09  0 -7.1054e-14

#adding colour.difference:sex
modc <- glm.WH(cbind(corr,incorr)~Colour.difference+ind+Colour.difference:sex, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(modc)
#[1] "FULL RANK"
anova(mod.8a,modc,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + ind + Colour.difference:sex
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        74     121.09                     
#2        73     120.95  1  0.14278   0.7055

#adding background:colour diff # DOESN'T WORK
modd <- glm.WH(cbind(corr,incorr)~Colour.difference+background+ind+Colour.difference:background, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
anova(mod.8a,modd,test="Chisq")



#adding background:sex
mode <- glm.WH(cbind(corr,incorr)~Colour.difference+ind+background:sex, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mode)
#[1,] "backgroundgreen:sexfemale" 
#[2,] "backgroundorange:sexfemale"
#[3,] "backgroundgreen:sexmale"   
#[4,] "backgroundorange:sexmale" 

anova(mod.8a,mode,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference + ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + ind + background:sex
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        74     121.09                     
#2        74     121.09  0        0 


#adding colour.diff:background:sex = DOESN'T WORK
modf <- glm.WH(cbind(corr,incorr)~Colour.difference+ind+Colour.difference:background:sex, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
anova(mod.8a,modf,test="Chisq")











