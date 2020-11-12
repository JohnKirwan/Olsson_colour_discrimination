rm(list = ls()) #clear the workspace of anything before starting








## Loading packages required to run the analysis and the plotting
library(lattice)
library(graphics)
library(psyphy)
library(lmerTest)
library(MASS)
library(glm2)
# reading the data

load("(correct folder)/orange colours data.RData")
#or
#data <- read.csv('orange colours.csv', header = TRUE, sep = ";", row.names = NULL) # change the "Working Directory" to the folder containing your data file. Change the name xxxxx.csv to your filename
data$background <- factor(data$background)
data$ind <- factor(data$ind)
data$batch <- as.factor(data$batch)
## Setting up models

#Rank deficiency is a big problem for model selection
#i.e. The design is not balanced, so that at least one of the covariates can be written as an exact linear combination of other covariates.

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

#Orange colours

#Model testing, null model and full model, and models removing one factor
mod.0 <- glm.WH(cbind(corr,incorr)~1, data = data, lambda.init = 0.1, interval = c(0, 1), NumAlt = 2) # Null model, data is explained only by variance. lambda.init is your estimate of the upper asymptote and NumAlt specifies the number of alternatives e.g. 2 alternative forced choice
CheckRank(mod.0)
#[1] "FULL RANK"

# full model with 3-way interaction
mod.1 <- glm.WH(cbind(corr,incorr)~Colour.difference*background*sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=1000) # Model where data is affected by the stimulus level
CheckRank(mod.1)
#[1,] "indD3"                   
#[2,] "indD4"                   
#[3,] "batchb"                  
#[4,] "backgroundorange:sexmale"

anova(mod.0,mod.1,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ 1
#Model 2: cbind(corr, incorr) ~ Colour.difference * background * sex + 
#  ind + batch
#Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
#1        78     725.99                          
#2        59      96.66 19   629.33 < 2.2e-16 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
AIC(mod.0,mod.1)
#df      AIC
#mod.0  1 972.5317
#mod.1 20 381.2061

#removing 3-way interaction
mod.2 <- glm.WH(cbind(corr,incorr)~Colour.difference*background+Colour.difference*sex+background*sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.2)
#[1,] "indD3"                   
#[2,] "indD4"                   
#[3,] "batchb"                  
#[4,] "backgroundorange:sexmale"

anova(mod.1,mod.2, test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background * sex + 
#  ind + batch
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + background * sex + ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        59     96.664                     
#2        60     96.801 -1 -0.13684   0.7114
AIC(mod.1,mod.2)
#df      AIC
#mod.0  1 972.5317
#mod.1 20 381.2061


#removing background:sex
mod.21 <- glm.WH(cbind(corr,incorr)~Colour.difference*background+Colour.difference*sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.21)
  
#[1,] "indD3" 
#[2,] "indD4" 
#[3,] "batchb"

anova(mod.2,mod.21,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + background * sex + ind + batch
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + background + sex + ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        58     93.566                     
#2        58     93.566  0        0 

AIC(mod.2,mod.21)
#df      AIC
#mod.2  19 379.3429
#mod.21 19 379.3429


## Removing colour.difference:sex
mod.22 <- glm.WH(cbind(corr,incorr)~Colour.difference*background+sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.22)
#[1,] "indD3" 
#[2,] "indD4" 
#[3,] "batchb"

anova(mod.2,mod.22,test="Chisq")
#Analysis of Deviance Table
#
#Model 1: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + background * sex + ind + batch
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference + 
#  sex + background * sex + ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        58     93.566                     
#2        59     93.823 -1 -0.25705   0.6122

AIC(mod.2,mod.22)
#df      AIC
#mod.2  19 379.3429
#mod.22 18 377.4202

#removing colour.difference:background
mod.23 <- glm.WH(cbind(corr,incorr)~Colour.difference+background+sex+ind+batch, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.23)
#[1,] "indD3" 
#[2,] "indD4" 
#[3,] "batchb"

anova(mod.22,mod.23,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + background * sex + ind + batch
# 2: cbind(corr, incorr) ~ Colour.difference + background + Colour.difference * 
 # sex + background * sex + ind + batch
#Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
#1        58     93.566                          
#2        59    133.933 -1  -40.368 2.104e-10 ***

AIC(mod.22,mod.23)
#df      AIC
#mod.22 18 377.4202
#mod.23 17 410.2438

#removing batch
mod.3 <- glm.WH(cbind(corr,incorr)~Colour.difference*background+sex+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.3)
#[1,] "indD3"
#[2,] "indD4"

anova(mod.22,mod.3,test="Chisq")

#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference * 
#  sex + background + sex + ind + batch
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + Colour.difference + 
#  sex + background + sex + ind + batch
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        58     93.566                     
#2        59     93.823 -1 -0.25705   0.6122
AIC(mod.22,mod.3)
#df      AIC
#mod.22 18 377.4202
#mod.3  18 377.4202

#removing sex
mod.31 <- glm.WH(cbind(corr,incorr)~Colour.difference*background+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.31)
#[1,] "indD4"

anova(mod.3,mod.31,test="Chisq")

#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + sex + 
#  ind
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + ind
#Resid. Df Resid. Dev Df   Deviance Pr(>Chi)
#1        61     96.878                       
#2        61     96.878  0 1.5632e-13

AIC(mod.3,mod.31)
#df      AIC
#mod.3  18 377.4202
#mod.31 18 377.4202

## Removing ind
mod.4 <- glm.WH(cbind(corr,incorr)~Colour.difference*background, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.4)
#[1] "FULL RANK"

anova(mod.31,mod.4,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + sex + 
#  ind
#Model 2: cbind(corr, incorr) ~ Colour.difference * background
#Resid. Df Resid. Dev  Df Deviance  Pr(>Chi)    
#1        61     96.878                           
#2        75    138.715 -14  -41.837 0.0001312 ***
AIC(mod.31,mod.4)

#df      AIC
#mod.31 18 365.8286
#mod.4  18 365.8286



## removing colour.difference:background interaction
mod.5 <- glm.WH(cbind(corr,incorr)~Colour.difference+background+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.5)
#[1,] "indD4"
anova(mod.31,mod.5,test="Chisq")

#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + ind
#Model 2: cbind(corr, incorr) ~ Colour.difference + background + ind
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)    
#1        61     96.878                         
#2        62    131.701 -1  -34.824 3.61e-09 ***

AIC(mod.31,mod.5)
#df      AIC
#mod.31 18 377.4202
#mod.5  17 410.2438



summary(mod.31)

#Call:
#  glm(formula = formula, family = binomial(probit.lambda(NumAlt, 
 #                                                        lam.c)), data = data, maxit = 100)

#Deviance Residuals: 
#  Min        1Q    Median        3Q       Max  
#-3.23494  -0.60308   0.05277   0.60417   2.62016  

#Coefficients: (1 not defined because of singularities)
#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                         -2.4733     0.4179  -5.918 3.25e-09 ***
#  Colour.difference                    1.1681     0.1544   7.563 3.93e-14 ***
#  backgroundorange                    -2.5515     1.1018  -2.316 0.020572 *  
#  indaA3                              -0.1029     0.6247  -0.165 0.869211    
#indaA4                               1.8000     0.5757   3.127 0.001768 ** 
#  indaA5                               0.1428     0.5969   0.239 0.810888    
#indbA12                             -0.3762     0.6597  -0.570 0.568467    
#indbA2                               0.3379     0.5724   0.590 0.554985    
#indbA52                             -0.1045     0.6214  -0.168 0.866512    
#indbA6                               0.2956     0.5770   0.512 0.608423    
#indC1                                0.6186     0.3743   1.652 0.098445 .  
#indC2                                0.2994     0.3826   0.783 0.433916    
#indC3                                0.6271     0.3742   1.676 0.093793 .  
#indC4                                0.5101     0.3762   1.356 0.175168    
#indD1                                1.4122     0.3818   3.699 0.000216 ***
#  indD2                                0.2510     0.3846   0.653 0.514063    
#indD3                                0.6331     0.3741   1.692 0.090645 .  
#indD4                                    NA         NA      NA       NA    
#Colour.difference:backgroundorange   3.4175     0.9124   3.745 0.000180 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#(Dispersion parameter for binomial family taken to be 1)

#Null deviance: 725.989  on 79  degrees of freedom
#Residual deviance:  96.878  on 61  degrees of freedom
#AIC: 377.42

#Number of Fisher Scoring iterations: 100

#lambda	 0.0407 	

	


## Forward Selection below

#Adding Colour.difference:sex
mod.31a <- glm.WH(cbind(corr,incorr)~Colour.difference*background+ind+Colour.difference:sex, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.31a)
#[1,] "indD4"
anova(mod.31a,mod.31,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + ind + 
#  Colour.difference:sex
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + ind
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        60     96.801                     
#2        61     96.878 -1  -0.0773    0.781

AIC(mod.31,mod.31a)
#df      AIC
#mod.31  18 377.4202
#mod.31a 19 379.3429  

#Adding background:sex
mod.31b <- glm.WH(cbind(corr,incorr)~Colour.difference*background+ind+background:sex, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level
CheckRank(mod.31b)
                   
#[1,] "indD4"                   
#[2,] "backgroundgreen:sexmale" 
#[3,] "backgroundorange:sexmale"

anova(mod.31b,mod.31,test="Chisq")
#Analysis of Deviance Table

#Model 1: cbind(corr, incorr) ~ Colour.difference * background + ind + 
#  background:sex
#Model 2: cbind(corr, incorr) ~ Colour.difference * background + ind
#Resid. Df Resid. Dev Df Deviance Pr(>Chi)
#1        61     96.878                     
#2        61     96.878  0        0         

AIC(mod.31,mod.31b)
#df      AIC
#mod.31  18 377.4202
#mod.31b 18 377.4202


#Adding colour.difference:background:sex = DOESN'T Work
mod.31d <- glm.WH(cbind(corr,incorr)~Colour.difference*background*sex+ind, data = data, lambda.init = 0.05, interval = c(0, 1), NumAlt = 2, maxit=100) # Model where data is affected by the stimulus level











