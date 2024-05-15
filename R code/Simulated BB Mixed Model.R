## Simulated Smoking Cessation Data
## Random Intercept Model/ Random Slope Model 

library(VGAM)    # rbetabinom()
library(PROreg)  # bbmm()
library(lme4)    # Create random effect disgn matrix
library(mvtnorm) # rmvnorm()

rm(list = ls())
#------------------------------------------------------------------------------#

#------------------------------------#
# Generate BB Random Intercept Model #
#------------------------------------#


set.seed(2202023)
ntrial<-7                        # prior week (cluster size)
n<-250                           # total subjects
nis<-sample(1:12,n,T)            # repeated measure per subject (unbalanced)
id<-rep(1:n,nis)                 # id
N<-length(id)                    # total observations

# Fixed effect 
t<-rep(0,N)
for (i in 1:n) t[id==i]<-sort(sample(1:12,nis[i]))  # time variable
X<-cbind(1,t)                       # Design Matrix
beta<-c(-2.5,0.35)                  # Fixed effects

truerho<-rho<-0.25                  # correlation parameter (0<rho<1)

# Random effect
truesigmab<-sigmab<-0.5
trueB<-b<-rnorm(n,sd=sqrt(sigmab))  # random intercept 

# BB Outcomes 
eta<-X%*%beta+rep(b,nis)
mu<-exp(eta)/(1+exp(eta))


y<-rbetabinom(N,size=ntrial,prob=mu,rho=rho) # outcome

# Dataset
dat<-data.frame(id=id,y=y,m=7,t=t)
dat$id<-factor(dat$id)

#------------------------------------------------------------------------------#


#------------------------#
# Random Intercept Model #
#------------------------#

model <- BBmm(fixed.formula = y~t,
              random.formula = ~id,m=7,data=dat)
summary(model)
# Call:	BBmm(fixed.formula = y ~ t, random.formula = ~id, m = 7, data = dat)
# 
# Fixed effects coefficients:
#   
#   Estimate    StdErr t.value   p.value    
# (Intercept) -2.341474  0.079102 -29.601 < 2.2e-16 ***
#   t            0.344399  0.010777  31.956 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# ---------------------------------------------------------------
#   Random effects dispersion parameter(s):
#   
#   Estimate     StdErr
# id 0.768242 0.03868365
# 
# ---------------------------------------------------------------
#   Logarithm of beta-binomial dispersion parameter log(phi):
#   
#   Estimate     StdErr
# 1 -1.187925 0.06270429    # exp(-1.187925)/(1+exp(-1.187925))=0.23 (rho)
# 



#--------------------------------------#
# Random Intercept Model (Alternative) #
#--------------------------------------#

# Generate random effect design matrix
randformula = "y ~ t  + (1 | id)"
tmp <- lFormula(eval(randformula), dat)
Z <- t(as.matrix(tmp$reTrms$Zt))

model2 <-BBmm(X=X,y=y,Z=Z,m=7,nRandComp=n)
summary(model2)  # same as model 1

# Call:	BBmm(X = X, y = y, Z = Z, nRandComp = n, m = 7)
# 
# Fixed effects coefficients:
#   
#   Estimate    StdErr t.value   p.value    
# -2.341474  0.079102 -29.601 < 2.2e-16 ***
#   t  0.344399  0.010777  31.956 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# ---------------------------------------------------------------
#   Random effects dispersion parameter(s):
#   
#   Estimate     StdErr
# 1 0.768242 0.03868365
# 
# ---------------------------------------------------------------
#   Logarithm of beta-binomial dispersion parameter log(phi):
#   
#   Estimate     StdErr
# 1 -1.187925 0.06270429



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
rm(list = ls())
##-------------------------------#
# Generate BB Random Slope Model #
#--------------------------------#
set.seed(2202023)
ntrial<-7                        # prior week (cluster size)
n<-250                           # total subjects
nis<-sample(1:12,n,T)            # repeated measure per subject (unbalanced)
id<-rep(1:n,nis)                 # id
N<-length(id)                    # total observations

# Fixed effect 
t<-rep(0,N)
for (i in 1:n) t[id==i]<-sort(sample(1:12,nis[i]))  # time variable
X<-cbind(1,t)                       # Design Matrix
beta<-c(-2.5,0.35)                  # Fixed effects

truerho<-rho<-0.25                  # correlation parameter (0<rho<1)

# Random effect
truesigmab<-sigmab<-matrix(c(0.5,0.10,0.10,0.25),2,2)
trueB<-B<-rmvnorm(n,sigma=sigmab)  # random intercept/slope

b1<-B[,1]                    # Random intercepts
b2<-B[,2]                    # Random slopes

# BB Outcomes 
eta<-X%*%beta+rep(b1,nis)+rep(b2,nis)*t
mu<-exp(eta)/(1+exp(eta))

y<-rbetabinom(N,size=ntrial,prob=mu,rho=rho) # outcome

# Dataset
dat<-data.frame(id=id,y=y,m=7,t=t)
dat$id<-factor(dat$id)


#--------------------#
# Random Slope Model #
#--------------------#
randformula = "y ~ t + (1 + t | id)"
tmp <- lFormula(eval(randformula), dat)
Z <- t(as.matrix(tmp$reTrms$Zt))
Z1<-Z[,seq(1,dim(Z)[2],by=2)]    # Random Intercpet matrix
Z2<-Z[,seq(2,dim(Z)[2],by=2)]    # Random Slope matrix
Zstar<-cbind(Z1,Z2)

model3 <-BBmm(X=X,y=y,Z=Zstar,m=7,nRandComp=c(n,n))
summary(model3)  

# Call:	BBmm(X = X, y = y, Z = Zstar, nRandComp = c(n, n), m = 7)
# 
# Fixed effects coefficients:
#   
#   Estimate    StdErr t.value   p.value    
# -2.264760  0.092396 -24.512 < 2.2e-16 ***
#   t  0.367935  0.014832  24.808 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# ---------------------------------------------------------------
#   Random effects dispersion parameter(s):
#   
#   Estimate     StdErr
# 1 0.8024727 0.06914888
# 2 0.5489055 0.03017336
# 
# ---------------------------------------------------------------
#   Logarithm of beta-binomial dispersion parameter log(phi):
#   
#   Estimate     StdErr
# 1 -1.134703 0.08193738






