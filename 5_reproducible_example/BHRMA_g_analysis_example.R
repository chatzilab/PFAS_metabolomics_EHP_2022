### BHRMA: Script to perform Bayesian Ridge Regression with g prior and selection for mixtures
# authors: Jingxuan He, Jesse Goodrich, and David Conti
# last updated: Mar 2022

library(R2jags)
library(tidyverse)


# Description:
# input variables
# X: A NxP matrix of exposures for mixture analysis (on the original scale with NA's for individuals with BLD)
# Y: A N-length vector for a continuous outcome
# U: A NxQ matrix of covariates (variables included in the regression model but not included in the g-estimation)
# LOD: A P-length vector of LODs for each exposure. Individuals with missing data will have data imputed below this level of detection  
# profiles: A 2xP matrix of two counterfactual profiles of exposures for which a potential outcomes risk difference is calculated (as the exposures are standardized within the function, these profiles should be on the standard normal scale)

# Create Function ---------------------------------------------
BHRMA.g <- function(X=NULL, Y=NULL, U=NULL, LOD=NULL, profiles=NULL) {
  
  # JAGS model
  ridge.BDL.model <- 
    "model {
  for(i in 1:N) {
    Y[i] ~ dnorm(mu[i], prec.sigma.Y)
    mu[i] <- alpha + inprod(beta[1:P], X.s[i,1:P]) + inprod(delta[1:Q], U[i,1:Q])
    
    # imputation BDL
    for(p in 1:P) {
      X[i,p] ~ dnorm(X.true[i,p],prec.X[p]) 
      X.true[i,p] <- X.notmiss[i,p]*(1-R[i,p]) + X.miss[i,p]*R[i,p]
      X.notmiss[i,p] ~ dnorm(mu.X[p], tau.X[p])T(LOD[p], )
      X.miss[i,p] ~ dnorm(mu.X[p], tau.X[p])T( , LOD[p])
      X.s[i,p] <- (X.true[i,p] - mu.X[p])/sigma.X[p]
    }
  }
  # prior on outcome variance
  prec.sigma.Y <- 1/(sigma.Y*sigma.Y)
  sigma.Y ~ dunif(0,3)
  
  # prior on covariate effects
  for(q in 1:Q) { delta[q] ~ dnorm(0, 1.0E-06) }
  
  # prior on intercept
  alpha ~ dnorm(0, 1.0E-06)
  
  # prior on exposure effects
  beta[1:P] ~ dmnorm(mu.beta[1:P], T[1:P, 1:P])
  for(j in 1:P) {
    mu.beta[j] <- (1-gamma[j])*prop.mu.beta[j]
    b[j] <- beta[j]*gamma[j]
    gamma[j] ~ dbern(pi)
    for(k in 1:P) {
      T[j,k] <- gamma[j]*gamma[k]*XtX[j,k]/(G) + (1-gamma[j]*gamma[k])*equals(j,k)*pow(prop.sd.beta[j],-2)
    }
    tau.X[j] <- 1/(sigma.X[j]*sigma.X[j])
    sigma.X[j] ~ dunif(0,5)
    mu.X[j] ~ dnorm(0, 1.0E-06)
    prec.X[j] <- 10000
  }
  pi ~ dbeta(1,P) 
  #pi ~ dbeta(P, 1)

  # semi-Bayes
  #G <- w/(1-w)
  #w <- 0.9   # w -> 0 shrink to common mean; as w -> inf toward the maximum likelihood estimate

  # Zellner and Siow prior on G
  #b0 <- 0.5*N
  #inv.G ~ dgamma(0.5, b0)
  #G <- 1/inv.G
  #w <- G/(G+1)

  # Hyper-g prior (following Perrakis 2018, note that this is on the G^-1 so the Beta distribution is switchd in terms of a and b from Li and Clyde 2019 equation 34)
  a <- 3
  bw <- a/2 - 1
  w~dbeta(1,bw) T(0.00001,0.99999) # JG: added T(0.0001,0.9999) 
  G <- w/(1-w)

  # Hyper-g/n prior (following Perrakis 2018, note that this is on the G^-1 so the Beta distribution is switchd in terms of a and b from Li and Clyde 2019 equation 34)
  #a <- 3
  #bw <- a/2 - 1
  #w~dbeta(1,bw)
  #G <- N*w/(1-w)

  #beta-prime 
  #G <- w/(1-w)
  #w ~ dbeta(bw, .25)
  #bw <- (N-P_m-1.5)/2
  #P_m <- sum(gamma[1:P])

  # g-estimation
  eta.low <- inprod(b[1:P], profiles[1,1:P])
  eta.high <- inprod(b[1:P], profiles[2,1:P])
  psi <- eta.high-eta.low
  
}"
  
  # Other Stuff
  N <- length(Y)
  P <- ncol(X)
  Q <- ncol(U)
  R <- ifelse(is.na(X), 1,0)
  exposure.Names <- colnames(X)
  
  ### get the univariate result
  univariate.results <- t(sapply(1:P, FUN=function(p) {  # using index p facilitate write
    x <- as.matrix(X[,p])
    reg <- glm(Y~x, family=gaussian)    # perform logistic regression
    s.reg <- summary(reg)                 # get the summary for the regression
    c.reg <- s.reg$coef["x",]             # select the coefficients for the exposure
    return(c.reg)                         # to avoid potential memory issues only return coefficients if small number of exposures
  }, simplify=T))
  univariate.results <- data.frame(exposure.Names,univariate.results)
  
  ### g prior model result
  prop.mu.beta <- rep(0, P)
  prop.sd.beta <- univariate.results$Std..Error
  XtX <- t(as.matrix(X))%*%as.matrix(X) 
  
  # run jags
  jdata <- list(N=N, Y=Y, X=X, R=R, U=U, P=P, Q=Q, 
                profiles=profiles, LOD=LOD,XtX=XtX, 
                prop.mu.beta=prop.mu.beta, 
                prop.sd.beta=prop.sd.beta)
  
  var.s <- c("beta", "gamma", "eta.low", "eta.high",  "psi")
  model.fit <- jags.model(file=textConnection(ridge.BDL.model), data=jdata, n.chains=1, n.adapt=4000, quiet=T)
  update(model.fit, n.iter=1000, progress.bar="none")
  model.fit <- coda.samples(model=model.fit, variable.names=var.s, n.iter=5000, thin=1, progress.bar="none")
  
  # summarize results
  r <- summary(model.fit)
  var.names <- c(paste(exposure.Names, "beta", sep="."),
                 "eta.high",
                 "eta.low",
                 paste(exposure.Names, "gamma", sep="."),
                 "psi")
  ridge.BDL.results <- data.frame(var.names, 
                                  r$statistics[,1:2], 
                                  r$quantiles[,c(1,5)])
  wald = abs(ridge.BDL.results[,"Mean"]/ridge.BDL.results[,"SD"])
  ridge.BDL.results$p.val = (2*(1-pnorm(wald,0,1)))
  return(ridge.BDL.results)
}

# Set up Data  -------------------------------
simulated_data <- read_rds(here::here("5_reproducible_example",
                                      "simulated_data.RDS"))

# X: A NxP matrix of exposures for mixture analysis (on the original scale with NA's for individuals with BLD)
X = simulated_data %>% 
  select(pfos, pfhxs, pfhps, pfoa, pfna, pfda)

# P: the number of exposures in the mixture
P = ncol(X)

# Y: A N-length vector for a continuous outcome
Y = simulated_data$metabolite_1

# U: A NxQ matrix of covariates (variables included in the regression model but not included in the g-estimation)
U = simulated_data %>% 
  select(age, sex.num, bmi, tanner, ses.num, wave.num)

# LOD: A P-length vector of LODs for each exposure. Individuals with missing data will have data imputed below this level of detection  
LOD = c(0.01, 0.01, 0.05, 0.01, 0.01, 0.01, 0.43, 0.01)

# profiles: A 2xP matrix of two counterfactual profiles of exposures for which a potential outcomes risk difference is calculated (as the exposures are standardized within the function, these profiles should be on the standard normal scale)
profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P)

# Names of exposures
exposure.Names = colnames(X)

# Run Function ------------------------------------------------------
# should take ~4 minutes on a standard laptop computer
fit = BHRMA.g(X=X, Y=Y, U=U, LOD=LOD, profiles=profiles)

fit
