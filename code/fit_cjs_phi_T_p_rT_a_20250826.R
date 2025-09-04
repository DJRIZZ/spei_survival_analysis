# Fit Hierarchical Bayesian Cormack-Jolly-Seber model to mark-resight encounter histories from spectacled eiders at Kigigak 1992-2025
# no resight effort in 2015-2018 and 2020, with these years set to NA in the encounter history
# code is based on Kery and Schaub 2012, Chpt 7

# load data (encounter history created with file "create_encounterHistory_1992-2025_20250731.R")

# load packages
library(rjags)
#library(R2jags)
library("jagsUI") #package to bridge R to JAGS and work with JAGS output
library(ggplot2)

set.seed(12) # make random results reproducible

# load encounter history data
CH <- read.csv("data/ec_na_1992-2015_ad_ducklings_2019-2025_adults.csv")
CH <- CH[order(CH$id_metalBand), ] # sort descending by id_metalBand
# remove column names and id_metalBand for JAGS input
ch <- CH # new data frame to modify
colnames(ch) <- NULL # remove col names
ch<- ch[ ,2:35] # remove column id_metalBand

# create vector indicating occasion of first capture/marking for each individual
get.first<- function(x) min(which(x != 0))
f <- apply(ch, 1, get.first)
hist(f,
     breaks = seq(1,ncol(ch),1),
     freq = FALSE)

# create vector indicating individuals marked as ducklings
kt_eh <- read.csv("data/data_christie_paper/Encounter_hist_SPEI.csv") # load capture history from Christie et al. 2018 with age (duckling = 1) data for 1992-2015, when ducklings were marked
kt_keep <- c("WHOLE.BAND.NO", "Duckling") # relevant columns (metal band number, and duckling indicator variable column where ducklings = 1, 0 otherwise)
kt_eh <- kt_eh[ ,kt_keep] # reduce to relevant columns
rm(kt_keep) 
colnames(kt_eh) <- c("id_metalBand", "is_duckling")
ct_ducklings <- sum(kt_eh$is_duckling) # count of ducklings
all_metalBand <- as.data.frame(CH$id_metalBand) # all metal band numbers in the encounter history 1992-2025
colnames(all_metalBand) <- "id_metalBand" 
marked_as_duckling <- merge(all_metalBand, kt_eh, by = "id_metalBand", all.x = TRUE) # merge all metal bands with metal bands 1992-2015 with duckling indicator variable
marked_as_duckling <- marked_as_duckling[order(marked_as_duckling$id_metalBand), ] # sort descending by id_metalBand
is_duckling <- ifelse(is.na(marked_as_duckling$is_duckling), 0, marked_as_duckling$is_duckling) # add zero values indicating adults for all invids without data, which are those from 2019-2025
reftbl_duckling <- cbind.data.frame(marked_as_duckling, is_duckling)
ct2_ducklings <- sum(is_duckling) # check duckling count against original data (kt_eh), should be equal

# create matrix indicating third summer (earliest a female would return to breed) for individuals marked as ducklings and 0 otherwise
# for individs marked as ducklings, p_3rd summer will be estimated and no p will be estimated for 2nd summer (NA in matrix)
nind <- nrow(ch) # number of marked individuals in encounter history 1992-2025
n.occasions <- ncol(ch) # number of occasions 1992-2025 including years with no resight effort (2016-2018, 2020)

age_return <- matrix(0, nind, n.occasions) # create template matrix of all zeros

for (i in 1:nind) { # for each individual
  if (is_duckling[i] == 1) { # if it is indicated to be a duckling
    second_summer <- f[i] + 1 # assign indicator for 2nd summer when individ will not return to breeding area
    third_summer <- f[i] + 2  # assign indicator for third summer (summer of they turn 2 years old) as 2 years after marking occasion (f)
    if(second_summer <= n.occasions) { # if second summer is within the range of encounter occasions
      age_return[i, second_summer] <- 0 # assign NA as individ should not have estimated resight for 2nd summer (actually asign 0 to just not estimate)
    }
    if (third_summer <= n.occasions) { # if third summer is within the range of encounter occasions (i.e., don't project past range of data)
      age_return[i, third_summer] <- 1 # create indicator value of 1 for the third summer of individs marked as ducklings
    }
  }
  # adults stay 0
}

# set capture history data to y
y <- as.matrix(ch)

###########################################################################
# CJS model: time-specific fixed phi and RANDOM time p and fixed age p ####
###########################################################################
# write to JAGS
writeLines("
    model {
  # Priors for Phi: time-varying fixed effect (each phi[t] is independent with its own prior)
  for (t in 1:(n.occasions - 1)) { # indexed to n.occasions-1 because phi and p are defined between t
    logit(phi[t]) <- beta_phi[t]   # assign logit survival with beta_phi (logit phi) for each occasion
    beta_phi[t] ~ dnorm(0, 0.001)  # draw vague prior for time-specific fixed effects for phi using normal distribution with mean 0 and precision 0.001
  }

  # Priors for p: time-varying random effect and fixed age effect
  for (t in 1:(n.occasions - 1)) {
    alpha_p[t] ~ dnorm(0, tau_p)  # draw prior for random effect for time from normal distribution with mean 0 and precision tau_p; each alpha_p[t] is drawn from the same distirbition
  }

  tau_p <- pow(sigma_p, -2)       # assign precision (inverse variance) for time random effect using sigma_p where 1/sigma_p^2 is the precision; sigma_p^2 is the variance 
  mu_p ~ dnorm(0, 0.001)          # draw prior for global intercept for p from normal distribution with mean 0 and precision 0.001
  beta_age ~ dnorm(0, 0.001)      # draw priors for fixed effect of age from normal distribution with mean 0 and precision 0.001
  sigma_p ~ dunif(0, 5)           # draw SD hyperparameter (prior for a prior) of time random effect from uniform distribution between 0 and 5; influences how much p varies over t

# priors: beta_phi[t], alpha_p[t], tau_p, mu_p, beta_age, sigma_p (where t is 34 occasions)

  # Likelihood
  for (i in 1:nind) {
    z[i, f[i]] <- 1               # Latent state conditional on first capture 
     # ---- State process ----
    for (t in (f[i] + 1):n.occasions) { # indexed to n.occasions because individs are observed at each t
      mu1[i,t] <- phi[t - 1] * z[i, t - 1] # assign phi for preceding occasion conditional on z of preceding occasion (dead individs stay dead)
      z[i,t] ~ dbern(mu1[i,t])             # draw state z from Bernoulli distrib with probability mu1 for each individ and occasion
      # ---- Observation process ----
      logit(p[i, t]) <- mu_p + alpha_p[t - 1] + beta_age * age_return[i, t] # assign model for detection prob (p)
      mu2[i, t] <- p[i, t] * z[i, t] # assign effective detection prob (i.e., prob that individ i is observed at time t) conditional on individ is alive (z = 1; can't detect dead individs)
      y[i, t] ~ dbern(mu2[i, t])     # draw y (obsevation of capture 1 or non-capture 0) with probability mu2, same as p[i, t] if individ is alive; links observation process to the state process
    } # close time
  } # close individuals
    } # close model statement
    ", con ="model.txt") # save model script as model.txt file (referenced in jags command later) and close writeLines function
# parameters estimated: mu_p, alpha_p[t-1], beta_age...

phi_T_p_randT_age <- 'model.txt' # save model as text file

# create matrix of known z states to include as data for the model (z = 1 after first capture to last encounter, NA otherwise, including NA for f)
z_known <- function(ch){ # encounter history ch as input
  ch_temp <- ch # change object name for tweaking here
  state <- as.matrix(ch_temp) # create state from ch as starting matrix to modify with known z information
  state[is.na(state)] <- 0 # turn columns of NAs in years of no resight effort to 0 for this function because all zeros will eventually be turned to NA unless they occur between two 1's in a row
  for (i in 1:dim(ch_temp)[1]){ # for each individual (row in state)
    n1 <- min(which(state[i, ] == 1)) # assign n1 as the minimum occasion with a 1 
    n2 <- max(which(state[i, ] == 1)) # assign n2 as the maximum occasion with a 1
    state[i, n1:n2] <- 1 # assign 1s between the min and max occasions with 1 (individual known alive even if not seen)
    state[i, n1] <- NA # assign NAs from the first occasion to n1
  }
  state[state == 0] <- NA # all occasions with 0 in ch assign NA in state
  return(state)
}

# Initial values ####
# z initial values based on ch, all known z changed to NA, keep ch values otherwise (i.e., 0, or consider changing to 1's?)
z_init <- function(ch) { # capture history as input
  ch[is.na(ch)] <- 0  # Replace NA with 0 (occasions with no resight effort)
  ch <- as.matrix(ch) # convert to matrix to avoid errors caused by data frame subsetting
  for (i in 1:nrow(ch)) { # for each individual, row
    if (sum(ch[i, ]) == 1) { # for individuals without any resights, skip
      next
    }
    n1 <- min(which(ch[i, ] == 1)) # assign n1 for min occasion with 1 in the row
    n2 <- max(which(ch[i, ] == 1)) # assign n2 for max occasion with 1 in the row
   ch[i, n1:n2] <- NA # assign NA between n1 and n2
    ch[i, n1] <- NA # assign NA for first capture
  }
  for(i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA # assign NA for occasions before first capture
  }
  return(ch) # output the modified ch_temp matrix as z_init
}

  # all required initial values for estimated parameters as starting points for the estimation
  inits <- function() {
    list(
    z = z_init(ch),
    #z = z.init,
    beta_phi = rnorm(n.occasions-1, 0, 1),   # logit-scale survival effects, random init values centered on 0
    mu_p = rnorm(1, 0, 1), # logit-scale intercept for detection model, random init values centered on 0
    beta_age = rnorm(1, 0, 1), # age effect for detection model
    alpha_p = rnorm(n.occasions - 1, 0, 1), # random time effect on detection
    sigma_p = runif(1, 0.1, 1) # SD of random time effect on detection, init must be positive, set as small positive value
  )
}

# parameters to save
params <- c("phi", "beta_phi", "mu_p", "beta_age", "sigma_p", "alpha_p") # "deviance", "z"

# run the model ####
# set MCMC parameters for test run and run model in JAGS with jags function
#results <- jags(data = list(y = ch,
#                            age_return = age_return,
#                            n.occasions = n.occasions, 
#                            nind = nind,
#                            f = f,
#                            z = z_known(ch)),
#                inits = inits,
#                parameters.to.save = params,# what parameters to monitor for each iteration 
#                model.file = "model.txt", # where to find the model code
#                n.chains = 2, # 2 independent chains to permit assessing R hat
#                n.iter = 1000, # run each chain for 10k iterations
#                n.burnin = 200, # discard the first 2k iterations
#                n.thin = 1, # keep every 5th value
#                parallel = TRUE) # run chains in parallel with multiple processor cores

# set MCMC parameters for inferential run and run model in JAGS with jags function
results <- jags(data = list(y = ch,
                            age_return = age_return,
                            n.occasions = n.occasions, 
                            nind = nind,
                            f = f,
                            z = z_known(ch)),
                inits = inits,
                parameters.to.save = params,# what parameters to monitor for each iteration 
                model.file = "model.txt", # where to find the model code
                n.chains = 4, # 3 independent chains to permit assessing R hat
                n.iter = 100000, # run each chain for 10k iterations, 50K
                n.burnin = 20000, # discard the first 2k iterations
                n.thin = 5, # keep every 5th value
                parallel = TRUE) # run chains in parallel with multiple processor cores
#                progress.bar = "text") # show progress bar for simulation in the console

# get summaries for phi
phi_summary <- results$summary[grep("^phi\\[", rownames(results$summary)), ]

# format phi summary results as data frame
phi_summary_df <- as.data.frame(phi_summary)
phi_summary_df$year <- 1:(n.occasions - 1)

# get mean and 95% credible interval
phi_summary_df[, c("year", "mean", "2.5%", "97.5%")]

# plot annual phi estimates and credible intervals using ggplot2
p.phi <- ggplot(phi_summary_df, aes(x = year, y = mean)) +
  #geom_bar() +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2, color = "black") + # Add 95% CI error bars
  #geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.2) +
  labs(x = "Year", y = "Survival Probability (phi)",
       title = "Annual Apparent Survival Estimates") +
  theme_minimal()
p.phi


summary(results)
results
traceplot(results)

today <- Sys.Date()
save(results, file = paste0("output/cjs_model_output", today,".RData"))

summary_stats <- results$summary  

write.csv(summary_stats, paste0("output/cjs_model_summary", today, ".csv"), row.names = TRUE)

update1_results <- update(object = results, n.iter = 1000)
