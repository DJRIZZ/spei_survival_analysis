# Hierarchical CJS model in JAGS
# code by Cat Bradley

library(jagsUI)
# load data
orig.dat <- read.table("SPEI-92thru15-Adult.txt", header = F)
CH <- as.matrix(orig.dat)
## FUNCTION TO CONVERT CH TO M-ARRAY ##
marray <- function(CH) {
  n.ind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions + 1, nrow = n.occasions)
  # calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  } # t
  for (i in 1:n.ind){
    pos <- which(CH[i,] != 0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z], pos[z+1]] <- m.array[pos[z], pos[z+1]] + 1
    } # z
  } # i
  # calculate the number of individuals never recaptured
  for (t in 1:n.occasions){
    m.array[t, n.occasions + 1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  } # t
  out <- m.array[1:(n.occasions - 1), 2:(n.occasions + 1)]
  return(out)
}
# create the m-array from the capture histories
marr <- marray(CH)
releases <- apply(marr, 1, sum)

# minimum sea ice vector; 39 min.ice days in the winter from 1992 - 1993...
# 93 min.ice days in the winter 2014 - 2015
min.ice <- c(39, 49, 37, 52, 68, 48, 61, 28, 114, 43, 74, 74, 67, 33, 55, 54, 48, 43, 64, 45, 44, 65, 93)
min.ice.scale <- as.vector(scale(min.ice))
min.ice2.scale <- as.vector(scale(min.ice^2))

### random time effect on detection, covariate and residual variance on survival model
sink("cjs-RET-RETcov-mnl.mod")
cat("
  model{
  # Priors and constraints
  for (t in 1:(noccasions - 1)){
    logit(phi[t]) <- mu.phi + gamma*min.ice.scale[t] + gamma2*min.ice2.scale[t] + epsilon.phi[t]
    logit(p[t]) <- mu.p + epsilon.p[t]
  } # t

  for (t in 1:(noccasions - 1)){
    epsilon.phi[t] ~ dnorm(0, tau.phi)
    phi.est[t] <- 1/(1+exp(-mu.phi - gamma*min.ice.scale[t] - gamma2*min.ice2.scale[t] - epsilon.phi[t]))
    epsilon.p[t] ~ dnorm(0, tau.p)
    p.est[t] <- 1/(1+exp(-mu.p - epsilon.p[t]))
  } # t
  
  mean.phi ~ dunif(0,1) # prior for mean survival
  mu.phi <- log(mean.phi/(1-mean.phi)) # logit transformation
  sigma.phi ~ dunif(0, 4)
  tau.phi <- pow(sigma.phi, -2)
  sigma2.phi <- pow(sigma.phi, 2) # temporal variance

  mean.p ~ dunif(0,1) # prior for mean survival
  mu.p <- log(mean.p/(1-mean.p)) # logit transformation
  sigma.p ~ dunif(0, 4)
  tau.p <- pow(sigma.p, -2)
  sigma2.p <- pow(sigma.p, 2) # temporal variance

  gamma ~ dnorm(0, 0.001)I(-10,10)
  gamma2 ~ dnorm(0, 0.001)I(-10,10)

  # Define the multinomial likelihood
  for (t in 1:(noccasions - 1)){
    marr[t,1:noccasions] ~ dmulti(pr[t,], r[t])
  } # t

  # Define cell probabilities of the m-array
  # Main diagonal
  for (t in 1:(noccasions - 1)){
    q[t] <- 1 - p[t]
    pr[t,t] <- phi[t]*p[t]
  # Above main diagonal
    for (j in (t + 1):(noccasions - 1)){
      pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
  } # j
    # Below main diagonal
  for (j in 1:(t - 1)){
    pr[t,j] <- 0
  } # j
} # t
  # Last column: probability of non-recapture
  for (t in 1:(noccasions - 1)){
    pr[t,noccasions] <- 1 - sum(pr[t,1:(noccasions - 1)])
  } # t
}
", fill = TRUE)
sink()

# data
jags.data <- list(marr = marr, noccasions = dim(marr)[2], r = releases, min.ice.scale = min.ice.scale, min.ice2.scale = min.ice2.scale)

# inits
inits <- function() {list(mean.phi = runif(1, 0, 1),
                          sigma.phi = runif(1, 0, 4),
                          mean.p = runif(1, 0, 1),
                          sigma.p = runif(1, 0, 4),
                          gamma = runif(1, -5, 5),
                          gamma2 = runif(1, -5, 5))}
# params
parameters <- c("phi.est", "p.est", "sigma2.phi", "sigma2.p", "mu.phi", "gamma", "gamma2", "mu.p")
# settings
ni <- 25000
nt <- 1
nb <- 1000
nc <- 3
cjs.RET.RETcov.mnl.JAGS <- jags(jags.data, inits, parameters, "cjs-RET-RETcov-mnl.mod", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)
print(cjs.RET.RETcov.mnl.JAGS, 3)