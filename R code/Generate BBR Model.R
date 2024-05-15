library(VGAM) # rbetabinom()
library(aod)  # betabin()

rm(list=ls())

set.seed(1234)

n<-350          # Number of obs
m<-20           # Number of Bernolli trials (can by varied)

# Covariates
x1<-rnorm(n)  
x2<-rbinom(n,1,0.5)
X<-cbind(1,x1,x2)

# Fixed effects
beta<-c(-0.5,1.0,-0.5)

# Correlation parameter (0<rho<1)
rho<-0.25

# Linear predictors
eta<-X%*%beta
# Logit link
pi<-exp(eta)/(1+exp(eta))

# Generate outcomes
y<-rbetabinom(n=n,size=m,prob=pi,rho=rho)

# Create dataset
dat<-data.frame(y=y,m=m,x1=x1,x2=x2)

# Fit BB model from "aod" package
mod.aod<-betabin(cbind(y,m-y)~x1+x2,random=~1,data=dat)
summary(mod.aod)

# Fit BB model from "VGAM" package
mod.vgam <- vglm(cbind(y,m-y)~x1+x2,
                 family=betabinomial,
                 data = dat)
summary(mod.vgam)

exp(-1.24926)/(1+exp(-1.24926)) # = 0.228 # Names of linear predictors: logitlink(rho), default 


mod2.vgam <- vglm(cbind(y,m-y)~x1+x2,
                 family=betabinomial(lrho="loglink"),  # loglink(rho)
                 data = dat)
summary(mod2.vgam)
exp(-1.50136)
