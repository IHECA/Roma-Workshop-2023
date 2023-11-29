## IHE Workshop 2023 Toy ## VOI ### 
## Institute of Health Economics ##

#### HOUSEKEEPING ####
#getwd()
load("IHEworkshop2023ToyMarkov.RData")

#### STRUCTURAL PARAMETERS ####
pop <- 1   ## Define the  size of the population of interest, so the total value of perfect information is scaled to the population-level magnitude of the intervention. 

#### EXPECTED VALUE OF PERFECT INFORMATION ####
opt.strat.nb <- max.nb <- evpi <- matrix(0,              #Optimal monetary benefit, by willingness to pay
                                         nrow = nsim, ncol = n.wtp,
                                         dimnames = list(1:nsim,wtp))
# Based on the PSA, determine the optimal (most likely to CE) strategy for every willingness to pay threshold
# This optimal strategy is acting as the "chosen" strategy, assuming the decision maker chose the best strategy for a given WTP threshold
for (i in 1:n.wtp) {
  opt.strat.nb[,i] <- cea[1,max.col(ceac)[i],,i]
}

# Unless the probability of the optimal strategy being C-E is 100%, there is some possibility the chosen strategy was incorrect, due to parameter uncertainty.
## In some runs of our PSA, the non-optimal strategy's Net Benefit was higher than the chosen option,
## What's more, the proportion of these instances are likely changing based on WTP (see CEAC)

# So we also determine what the highest Net Benefit is for each run, for each WTP threshold
## Each 'Maximum NB' value represents the best strategy selected for each run - where the parameter values are known, so this is 'perfect' information
for (n in 1:nsim) {
  for (i in 1:n.wtp) {
    max.nb[n,i] <- max(cea[1,,n,i])
  }
}

# We can now calculate the loss due to uncertainty (choosing the optimal instead of the max)
### (max.nb - opt.strategy.nb) = the loss from imperfect information, per Run and per WTP.
# But if we summed up all those values as they are, we're 'double counting' the total loss for however  many runs we ran, so we divide by nsim
## This also has the benefit that certain loss results are more likely than others because of the distribution around each input parameter.
## So by dividing the full matrix by nsim, we're also 'weighting' the result by the joint probability density function of all input parameters.
## Finally, we multiply the results by the expected size of the population of interest
evpi <- colSums((max.nb - opt.strat.nb)/nsim * pop)

#That's it! Very little additional work was needed once we ensured our PSA results were properly formatted for VOI analysis
#Just the visuals remain.
matplot(evpi, type = 'l',
        xlab = "Willingness to Pay Threshold per QALY gained (,000 $)",
        ylab = "EVPI ($)",
        main = "Expected Value of Perfect Information",
        col= c("black"),
        lwd = 3)



#### EXPECTED VALUE OF PARTIAL PERFECT INFORMATION ####
n1sim <- 500   ##number of runs, first integration
n2sim <- 1000  ##number of runs, second integration
clock.evppi <- Sys.time()

## EVPPI-User Duplicate of Markov Trace & Transition Matrix
#### Note Sim parameter changed to match second integration Monte Carlo
evp.M.soc <- evp.M.sub <- array(0,
                                dim = c(n.t+1, n.state, n2sim),     
                                dimnames = list(0:n.t, state, 1:n2sim))
evp.M.soc[1,,] <- evp.M.sub[1,,] <- t0_dist    # set initial state distribution

evp.P.soc <- evp.P.sub <- array(0,
                                dim = c(n.state, n.state, n.t, n2sim),
                                dimnames = list(state,state,1:n.t, 1:n2sim))

## EVPPI-User Duplicate of EVPI calculators
evp.opt.strat.nb <- evp.max.nb <- matrix(0,              #Optimal monetary benefit, by willingness to pay
                                         nrow = n2sim, ncol = n.wtp,
                                         dimnames = list(1:n2sim,wtp))

## EVPPI Results matrix; now a 3D array to allow both MC integration
#### Build an array for each parameter / group for which you intend to derive EVPPI
#### Here we derive for each parameter group (Transition Probabilities, QALYs) and the Submitter's Treatment Effect parameter.
groups.evppi <- c("Transition Probs","QALYs","Submitter Tx")  ##Name the parameter groups for which you intend to generate EVPPI results

evppi <- matrix(0,     ##This will be our final results table, but there's some work to get to here.
                nrow = length(groups.evppi), ncol = n.wtp,
                dimnames = list(groups.evppi,wtp))

evppi.f.trans <- evppi.f.qalys <- evppi.f.subtx <- array(0,     ##Results array; for each second integration there is an EVPI matrix of first integration
                                                         dim = c(n1sim,n2sim,n.wtp),
                                                         dimnames = list(1:n1sim, 1:n2sim, wtp))

evppi.trans <- evppi.qalys <- evppi.subtx <- matrix(0,          ##For transforming the first integration results into the second integration EVPI; the EVPPI    
                                                    nrow = n1sim, ncol = n.wtp,
                                                    dimnames = list(1:n1sim,wtp))

# Each EVPPI analysis will require some minor adjustments from the base case build to allow two integration MC to run 
# The approach below is known as 'Double MCS' and is hardly efficient but does allow you to see precisely what is happening to the base case model and the effects on the EVP(P)I.
# Double MCS' biggest pro is it is very easy to validate and understand each step. The more elegant approaches include:
#### (1) Quadrature; Samples your PSA results according to each parameter's (binned) probability density function. Requires confidence in pdf integration
#### (2) Shefield's SAVI online tool that generates results based on a .csv file of your PSA results: https://savi.shef.ac.uk/SAVI/
#### (3) The DARTH Package employs a non-parametric regression method based on your PSA results, but requires confidence with R code: https://rdrr.io/github/DARTH-git/dampack/src/R/evppi.R

## EVPPI: TRANSITION PROBABILITIES ####
for (f in 1:n1sim) {     ## 'f' is First integration loop, generating second integration results for each first integration value.Note only adjustments are in parameter setting and the EVPPI array.
  #For Standard of Care  
  ##from Stable
  evp.P.soc["Stable","Stable",,]     <- 1-(p.stable_moderate[f]+p.stable_severe[f]+p.stable_dead[f]) 
  evp.P.soc["Stable","Moderate",,]   <- p.stable_moderate[f]
  evp.P.soc["Stable","Severe",,]     <- p.stable_severe[f]
  evp.P.soc["Stable","Dead",,]       <- p.stable_dead[f]
  ##from Moderate
  evp.P.soc["Moderate","Moderate",,] <- 1-(p.moderate_stable[f]+p.moderate_severe[f]+p.moderate_dead[f]) 
  evp.P.soc["Moderate","Stable",,]   <- p.moderate_stable[f]
  evp.P.soc["Moderate","Severe",,]   <- p.moderate_severe[f]
  evp.P.soc["Moderate","Dead",,]     <- p.moderate_dead[f]
  ##from severe
  evp.P.soc["Severe","Severe",,]     <- 1-(p.severe_moderate[f]+p.severe_stable[f]+p.severe_dead[f]) 
  evp.P.soc["Severe","Moderate",,]   <- p.severe_moderate[f]
  evp.P.soc["Severe","Stable",,]     <- p.severe_stable[f]
  evp.P.soc["Severe","Dead",,]       <- p.severe_dead[f]
  ##from Dead
  evp.P.soc["Dead","Dead",,]         <- 1
  
  #For Submitter Tx
  ##from Stable
  evp.P.sub["Stable","Stable",,]     <- 1-(p.stable_moderate[f]+p.stable_severe[f]+p.stable_dead[f]) 
  evp.P.sub["Stable","Moderate",,]   <- p.stable_moderate[f]
  evp.P.sub["Stable","Severe",,]     <- p.stable_severe[f]
  evp.P.sub["Stable","Dead",,]       <- p.stable_dead[f]
  ##from Moderate
  evp.P.sub["Moderate","Moderate",,] <- 1-(p.moderate_stable[f]+(p.moderate_severe[f]*tx.moderate_severe[1:n2sim])+p.moderate_dead[f]) 
  evp.P.sub["Moderate","Stable",,]   <- p.moderate_stable[f]
  evp.P.sub["Moderate","Severe",,]   <- p.moderate_severe[f]*tx.moderate_severe[1:n2sim]
  evp.P.sub["Moderate","Dead",,]     <- p.moderate_dead[f]
  ##from severe
  evp.P.sub["Severe","Severe",,]     <- 1-(p.severe_moderate[f]+p.severe_stable[f]+p.severe_dead[f]) 
  evp.P.sub["Severe","Moderate",,]   <- p.severe_moderate[f]
  evp.P.sub["Severe","Stable",,]     <- p.severe_stable[f]
  evp.P.sub["Severe","Dead",,]       <- p.severe_dead[f]
  ##from Dead
  evp.P.sub["Dead","Dead",,]         <- 1
  
  for (n in 1:n2sim) {                                      # loop through the number of simulation runs
    for (t in 1:n.t)  {                                     # loop through the number of cycles
      evp.M.soc[t+1,,n] <- evp.M.soc[t,,n] %*% evp.P.soc[,,t,n]   # estimate the state vector for the next cycle (t+1)
    }
    for (t in 1:n.t){                                      
      evp.M.sub[t+1,,n] <- evp.M.sub[t,,n] %*% evp.P.sub[,,t,n]    
    }
  }
  
  #PSA Costs and QALYs
  evp.tc.soc <- evp.tc.sub <- evp.tu.soc <- evp.tu.sub <- matrix(0, nrow = n.t+1, ncol = n2sim, dimnames = list(0:n.t,1:n2sim))
  for (n in 1:n2sim) {
    evp.tc.soc[,n] <- evp.M.soc[,,n] %*% c(c.stable,c.moderate,c.severe,c.dead)
    evp.tc.sub[,n] <- evp.M.sub[,,n] %*% c(c.stable,c.moderate+c.tx,c.severe,c.dead)
    evp.tu.soc[,n] <- evp.M.soc[,,n] %*% c(u.stable[n],u.moderate[n],u.severe[n],u.dead[n])
    evp.tu.sub[,n] <- evp.M.sub[,,n] %*% c(u.stable[n],u.moderate[n],u.severe[n],u.dead[n])
  }
  ##Discounting Costs and QALYs
  for (t in 1:n.t) {
    evp.tc.soc[t+1,] <- evp.tc.soc[t+1,] * disc.c[t]
    evp.tc.sub[t+1,] <- evp.tc.sub[t+1,] * disc.c[t]
    evp.tu.soc[t+1,] <- evp.tu.soc[t+1,] * disc.c[t]
    evp.tu.sub[t+1,] <- evp.tu.sub[t+1,] * disc.c[t]
  }
  
  evp.psa <- array(0,
                   dim = c(n.strat,n2sim,n.params),
                   dimnames = list(strategy, 1:n2sim, params))
  
  evp.psa["Standard of Care",,"Cost"] <- colSums(evp.tc.soc)
  evp.psa["Standard of Care",,"QALY"] <- colSums(evp.tu.soc)/c.t
  evp.psa["Submitted Tx",,"Cost"] <- colSums(evp.tc.sub)
  evp.psa["Submitted Tx",,"QALY"] <- colSums(evp.tu.sub)/c.t
  
  evp.cea <- array(0,
                   dim = c(length(v.cea),n.strat,n2sim,n.wtp),
                   dimnames = list(v.cea,strategy,1:n2sim,wtp))
  for (i in 1:n.wtp) {
    evp.cea["NB","Standard of Care",,i] <- evp.psa["Standard of Care",,"QALY"] * wtp[i] - evp.psa["Standard of Care",,"Cost"]
    evp.cea["NB","Submitted Tx",,i] <- evp.psa["Submitted Tx",,"QALY"] * wtp[i] - evp.psa["Submitted Tx",,"Cost"]
  }
  for (n in 1:n2sim) {
    for (i in 1:n.wtp) {
      if (evp.cea["NB","Submitted Tx",n,i] > evp.cea["NB","Standard of Care",n,i]) {
        evp.cea["Optimal","Submitted Tx",n,i] <- 1
      } else {
        evp.cea["Optimal","Standard of Care",n,i] <- 1
      }
    }
  }
  
  evp.ceac <- matrix(0, 
                     nrow = n.wtp, ncol = n.strat,
                     dimnames = list(wtp,strategy))
  for (i in 1:n.wtp) {
    evp.ceac[i,"Standard of Care"] <- mean(evp.cea["Optimal","Standard of Care",,i])
    evp.ceac[i,"Submitted Tx"]     <- mean(evp.cea["Optimal","Submitted Tx",,i])
  }
  
  ##Generating EVPI results, with first-integration fixed values on our parameter group
  for (i in 1:n.wtp) {
    evp.opt.strat.nb[,i] <- evp.cea[1,max.col(evp.ceac)[i],,i]
  }
  
  for (n in 1:n2sim) {
    for (i in 1:n.wtp) {
      evp.max.nb[n,i] <- max(evp.cea[1,,n,i])
    }
  }
  
  evppi.f.trans[f,,] <- (evp.max.nb - evp.opt.strat.nb)/n2sim * pop
}

#Combine second integration results to make value differences based on parameter group values only
for(n in 1:n1sim) {
  for(i in 1:n.wtp) {
    evppi.trans[n,i] <- sum(evppi.f.trans[n,,i])  
  }  
}

#Mean first integration results to get probabilistic distribution results
evppi["Transition Probs",] <- colMeans(evppi.trans)

matplot(evppi["Transition Probs",], type = 'l',
        xlab = "Willingness to Pay Threshold per QALY gained (,000 $)",
        ylab = "EVPI ($)",
        main = "Expected Value of Partial Perfect Information: Transition Probabilities",
        col= c("black"),
        lwd = 3)

run.time_evppi.trans <- Sys.time() - clock.evppi

print(run.time_evppi.trans)

#Saving for posterity
save(evppi.f.trans, file="evppi__trans.Rdata")



#### EVPPI: UTILITIES ####
for (f in 1:n1sim) {  
  for (n in 1:n2sim) {
    evp.tc.soc[,n] <- m.M.soc[,,n] %*% c(c.stable,c.moderate,c.severe,c.dead)
    evp.tc.sub[,n] <- m.M.sub[,,n] %*% c(c.stable,c.moderate+c.tx,c.severe,c.dead)
    evp.tu.soc[,n] <- m.M.soc[,,n] %*% c(u.stable[f],u.moderate[f],u.severe[f],u.dead[f])
    evp.tu.sub[,n] <- m.M.sub[,,n] %*% c(u.stable[f],u.moderate[f],u.severe[f],u.dead[f])
  }
  
  ##Discounting Costs and QALYs
  for (t in 1:n.t) {
    evp.tc.soc[t+1,] <- evp.tc.soc[t+1,] * disc.c[t]
    evp.tc.sub[t+1,] <- evp.tc.sub[t+1,] * disc.c[t]
    evp.tu.soc[t+1,] <- evp.tu.soc[t+1,] * disc.c[t]
    evp.tu.sub[t+1,] <- evp.tu.sub[t+1,] * disc.c[t]
  }
  
  evp.psa <- array(0,
                   dim = c(n.strat,n2sim,n.params),
                   dimnames = list(strategy, 1:n2sim, params))
  
  evp.psa["Standard of Care",,"Cost"] <- colSums(evp.tc.soc[,1:n2sim])
  evp.psa["Standard of Care",,"QALY"] <- colSums(evp.tu.soc[,1:n2sim])/c.t
  evp.psa["Submitted Tx",,"Cost"] <- colSums(evp.tc.sub[,1:n2sim])
  evp.psa["Submitted Tx",,"QALY"] <- colSums(evp.tu.sub[,1:n2sim])/c.t
  
  evp.cea <- array(0,
                   dim = c(length(v.cea),n.strat,n2sim,n.wtp),
                   dimnames = list(v.cea,strategy,1:n2sim,wtp))
  for (i in 1:n.wtp) {
    evp.cea["NB","Standard of Care",,i] <- evp.psa["Standard of Care",,"QALY"] * wtp[i] - evp.psa["Standard of Care",,"Cost"]
    evp.cea["NB","Submitted Tx",,i] <- evp.psa["Submitted Tx",,"QALY"] * wtp[i] - evp.psa["Submitted Tx",,"Cost"]
  }
  for (n in 1:n2sim) {
    for (i in 1:n.wtp) {
      if (evp.cea["NB","Submitted Tx",n,i] > evp.cea["NB","Standard of Care",n,i]) {
        evp.cea["Optimal","Submitted Tx",n,i] <- 1
      } else {
        evp.cea["Optimal","Standard of Care",n,i] <- 1
      }
    }
  }
  
  evp.ceac <- matrix(0, 
                     nrow = n.wtp, ncol = n.strat,
                     dimnames = list(wtp,strategy))
  for (i in 1:n.wtp) {
    evp.ceac[i,"Standard of Care"] <- mean(evp.cea["Optimal","Standard of Care",,i])
    evp.ceac[i,"Submitted Tx"]     <- mean(evp.cea["Optimal","Submitted Tx",,i])
  }
  
  ##Generating EVPI results, with first-integration fixed values on our parameter group
  for (i in 1:n.wtp) {
    evp.opt.strat.nb[,i] <- evp.cea[1,max.col(evp.ceac)[i],,i]
  }
  
  for (n in 1:n2sim) {
    for (i in 1:n.wtp) {
      evp.max.nb[n,i] <- max(evp.cea[1,,n,i])
    }
  }
  
  evppi.f.qalys[f,,] <- (evp.max.nb - evp.opt.strat.nb)/n2sim * pop
}

#Combine second integration results to make value differences based on parameter group values only
for(n in 1:n1sim) {
  for(i in 1:n.wtp) {
    evppi.qalys[n,i] <- sum(evppi.f.qalys[n,,i])  
  }  
}

#Mean first integration results to get probabilistic distribution results
evppi["QALYs",] <- colMeans(evppi.qalys[])


matplot(evppi["QALYs",], type = 'l',
        xlab = "Willingness to Pay Threshold per QALY gained (,000 $)",
        ylab = "EVPI ($)",
        main = "Expected Value of Partial Perfect Information: QALYs",
        col= c("black"),
        lwd = 3)

run.time_evppi.qalys <- Sys.time() - clock.evppi
print(run.time_evppi.qalys)

#Saving for posterity
save(evppi.f.qalys, file="evppi__qalys.Rdata")



#### EVPPI: TREATMENT EFFECT ####
for (f in 1:n1sim) {
  #For Standard of Care  
  ##from Stable
  evp.P.soc["Stable","Stable",,]     <- 1-(p.stable_moderate[1:n2sim]+p.stable_severe[1:n2sim]+p.stable_dead[1:n2sim]) 
  evp.P.soc["Stable","Moderate",,]   <- p.stable_moderate[1:n2sim]
  evp.P.soc["Stable","Severe",,]     <- p.stable_severe[1:n2sim]
  evp.P.soc["Stable","Dead",,]       <- p.stable_dead[1:n2sim]
  ##from Moderate
  evp.P.soc["Moderate","Moderate",,] <- 1-(p.moderate_stable[1:n2sim]+p.moderate_severe[1:n2sim]+p.moderate_dead[1:n2sim]) 
  evp.P.soc["Moderate","Stable",,]   <- p.moderate_stable[1:n2sim]
  evp.P.soc["Moderate","Severe",,]   <- p.moderate_severe[1:n2sim]
  evp.P.soc["Moderate","Dead",,]     <- p.moderate_dead[1:n2sim]
  ##from severe
  evp.P.soc["Severe","Severe",,]     <- 1-(p.severe_moderate[1:n2sim]+p.severe_stable[1:n2sim]+p.severe_dead[1:n2sim]) 
  evp.P.soc["Severe","Moderate",,]   <- p.severe_moderate[1:n2sim]
  evp.P.soc["Severe","Stable",,]     <- p.severe_stable[1:n2sim]
  evp.P.soc["Severe","Dead",,]       <- p.severe_dead[1:n2sim]
  ##from Dead
  evp.P.soc["Dead","Dead",,]         <- 1
  
  #For Submitter Tx
  ##from Stable
  evp.P.sub["Stable","Stable",,]     <- 1-(p.stable_moderate[1:n2sim]+p.stable_severe[1:n2sim]+p.stable_dead[1:n2sim]) 
  evp.P.sub["Stable","Moderate",,]   <- p.stable_moderate[1:n2sim]
  evp.P.sub["Stable","Severe",,]     <- p.stable_severe[1:n2sim]
  evp.P.sub["Stable","Dead",,]       <- p.stable_dead[1:n2sim]
  ##from Moderate
  evp.P.sub["Moderate","Moderate",,] <- 1-(p.moderate_stable[1:n2sim]+(p.moderate_severe[1:n2sim]*tx.moderate_severe[f])+p.moderate_dead[1:n2sim]) 
  evp.P.sub["Moderate","Stable",,]   <- p.moderate_stable[1:n2sim]
  evp.P.sub["Moderate","Severe",,]   <- p.moderate_severe[1:n2sim]*tx.moderate_severe[f]
  evp.P.sub["Moderate","Dead",,]     <- p.moderate_dead[1:n2sim]
  ##from severe
  evp.P.sub["Severe","Severe",,]     <- 1-(p.severe_moderate[1:n2sim]+p.severe_stable[1:n2sim]+p.severe_dead[1:n2sim]) 
  evp.P.sub["Severe","Moderate",,]   <- p.severe_moderate[1:n2sim]
  evp.P.sub["Severe","Stable",,]     <- p.severe_stable[1:n2sim]
  evp.P.sub["Severe","Dead",,]       <- p.severe_dead[1:n2sim]
  ##from Dead
  evp.P.sub["Dead","Dead",,]         <- 1
  
  for (n in 1:n2sim) {                                      # loop through the number of simulation runs
    for (t in 1:n.t)  {                                     # loop through the number of cycles
      evp.M.soc[t+1,,n] <- evp.M.soc[t,,n] %*% evp.P.soc[,,t,n]   # estimate the state vector for the next cycle (t+1)
    }
    for (t in 1:n.t){                                      
      evp.M.sub[t+1,,n] <- evp.M.sub[t,,n] %*% evp.P.sub[,,t,n]    
    }
  }
  
  #PSA Costs and QALYs
  evp.tc.soc <- evp.tc.sub <- evp.tu.soc <- evp.tu.sub <- matrix(0, nrow = n.t+1, ncol = n2sim, dimnames = list(0:n.t,1:n2sim))
  for (n in 1:n2sim) {
    evp.tc.soc[,n] <- evp.M.soc[,,n] %*% c(c.stable,c.moderate,c.severe,c.dead)
    evp.tc.sub[,n] <- evp.M.sub[,,n] %*% c(c.stable,c.moderate+c.tx,c.severe,c.dead)
    evp.tu.soc[,n] <- evp.M.soc[,,n] %*% c(u.stable[n],u.moderate[n],u.severe[n],u.dead[n])
    evp.tu.sub[,n] <- evp.M.sub[,,n] %*% c(u.stable[n],u.moderate[n],u.severe[n],u.dead[n])
  }
  ##Discounting Costs and QALYs
  for (t in 1:n.t) {
    evp.tc.soc[t+1,] <- evp.tc.soc[t+1,] * disc.c[t]
    evp.tc.sub[t+1,] <- evp.tc.sub[t+1,] * disc.c[t]
    evp.tu.soc[t+1,] <- evp.tu.soc[t+1,] * disc.c[t]
    evp.tu.sub[t+1,] <- evp.tu.sub[t+1,] * disc.c[t]
  }
  
  evp.psa <- array(0,
                   dim = c(n.strat,n2sim,n.params),
                   dimnames = list(strategy, 1:n2sim, params))
  
  evp.psa["Standard of Care",,"Cost"] <- colSums(evp.tc.soc)
  evp.psa["Standard of Care",,"QALY"] <- colSums(evp.tu.soc)/c.t
  evp.psa["Submitted Tx",,"Cost"] <- colSums(evp.tc.sub)
  evp.psa["Submitted Tx",,"QALY"] <- colSums(evp.tu.sub)/c.t
  
  evp.cea <- array(0,
                   dim = c(length(v.cea),n.strat,n2sim,n.wtp),
                   dimnames = list(v.cea,strategy,1:n2sim,wtp))
  for (i in 1:n.wtp) {
    evp.cea["NB","Standard of Care",,i] <- evp.psa["Standard of Care",,"QALY"] * wtp[i] - evp.psa["Standard of Care",,"Cost"]
    evp.cea["NB","Submitted Tx",,i] <- evp.psa["Submitted Tx",,"QALY"] * wtp[i] - evp.psa["Submitted Tx",,"Cost"]
  }
  for (n in 1:n2sim) {
    for (i in 1:n.wtp) {
      if (evp.cea["NB","Submitted Tx",n,i] > evp.cea["NB","Standard of Care",n,i]) {
        evp.cea["Optimal","Submitted Tx",n,i] <- 1
      } else {
        evp.cea["Optimal","Standard of Care",n,i] <- 1
      }
    }
  }
  
  evp.ceac <- matrix(0, 
                     nrow = n.wtp, ncol = n.strat,
                     dimnames = list(wtp,strategy))
  for (i in 1:n.wtp) {
    evp.ceac[i,"Standard of Care"] <- mean(evp.cea["Optimal","Standard of Care",,i])
    evp.ceac[i,"Submitted Tx"]     <- mean(evp.cea["Optimal","Submitted Tx",,i])
  }
  
  ##Generating EVPI results, with first-integration fixed values on our parameter group
  for (i in 1:n.wtp) {
    evp.opt.strat.nb[,i] <- evp.cea[1,max.col(evp.ceac)[i],,i]
  }
  
  for (n in 1:n2sim) {
    for (i in 1:n.wtp) {
      evp.max.nb[n,i] <- max(evp.cea[1,,n,i])
    }
  }
  
  evppi.f.subtx[f,,] <- (evp.max.nb - evp.opt.strat.nb)/n2sim * pop
}

#Combine second integration results to make value differences based on parameter group values only
for(n in 1:n1sim) {
  for(i in 1:n.wtp) {
    evppi.subtx[n,i] <- sum(evppi.f.subtx[n,,i])  
  }  
}

#Mean first integration results to get probabilistic distribution results
evppi["Submitter Tx",] <- colMeans(evppi.subtx)

matplot(evppi["Submitter Tx",], type = 'l',
        xlab = "Willingness to Pay Threshold per QALY gained (,000 $)",
        ylab = "EVPI ($)",
        main = "Expected Value of Partial Perfect Information: Submitter Tx Effect",
        col= c("black"),
        lwd = 3)

run.time_evppi.subtx <- Sys.time() - clock.evppi
print(run.time_evppi.subtx)

#Saving for posterity
save(evppi.f.subtx, file="evppi__subtx.Rdata")


#### GENERATE PLOTS ####
#df_evppi <- data.frame(evppi)
matplot(evppi[1,], type = 'l',
        ylim = c(0,2200),
        xlab = "Willingness to Pay Threshold per QALY gained (,000 $)",
        ylab = "EVPI ($)",
        main = "Expected Value of Partial Perfect Information",
        col  = "firebrick3",
        lwd  = 3) + 
  lines(evppi[2,], col = "royalblue1", lwd = 3) +
  lines(evppi[3,], col = "springgreen4", lwd = 3)
legend("top", groups.evppi, col= c("firebrick3", "royalblue1","springgreen4"),lty = 1:4, bty = "n")

# lines(evpi, col = "black", lwd = 3)

# Save full environment to avert re-runs
#save.image(file = "IHEworkshop2023 VOI.RData")
#load("IHEworkshop2023 VOI.RData")

