## IHE Workshop 2023 # Toy # Markov Model ### 
## Institute of Health Economics ##

#### HOUSEKEEPING ####
getwd()  ##Tells you what Working Directory you are currently saving to. Change with setwd, eg:
setwd("C:/Users/skatwyk/Desktop/VOIworkshop2023")

##Intall / Load packages with some useful (and common) custom commands that will make our lives easier later
#install.packages("dplyr")
library(dplyr)
2
#install.packages("devtools")
#library(devtools)
#install.packages("ggplot2")
library(ggplot2)

#### STRUCTURAL PARAMETERS ####
options(scipen = 999)
nsim   <- 10000      #Number of iterations run for probabilistic sensitivity analysis
n.t    <- 60        #Cycles - Time = Months
c.t    <- 12        #Number of Cycles in a Year
disc   <- 0.015/c.t  #Discount Rate
disc.c <- 1/(1 + disc)^(1:n.t) 
disc.u <- 1/(1 + disc)^(1:n.t)
wtp    <- seq(0,150000,by=1000)
wtp.str<- sprintf("%s",seq(0,150000,by=1000))
n.wtp  <- length(wtp)
clock  <- Sys.time()
set.seed(8709816)


#### INPUT PARAMETERS ####
strategy <- c("Standard of Care","Submitted Tx")   #Strategy Names
state    <- c("Stable","Moderate","Severe","Dead") #State Names
#Initial State Distribution
t0_dist  <- c("Stable"  = 0,
              "Moderate"= 0.70,
              "Severe"  = 0.30,
              "Dead"    = 0)
n.state  <- length(state)
n.strat  <- length(strategy)

#Transition Probability Input Parameters
##User-defined input parameters on the progression of observed persons, according to transition probability per cycle
##For simplification, all  other parameters are time-independent
mu.stable_moderate <- 0.075
mu.stable_severe   <- 0.01
mu.moderate_stable <- 0.04
mu.moderate_severe <- 0.12
mu.severe_stable   <- 0.08
mu.severe_moderate <- 0.001

var.p <- 0.1 #Toy simplification; all transition probabilities have same variance

##Beta Distribution Functions
estAlpha <- function(mu,var) {
  alpha <- ((1-mu)/var-1/mu)*mu^2
  return(alpha)
}
estBeta <- function(mu,var) {
  alpha <- ((1-mu)/var-1/mu)*mu^2
  beta  <- alpha*(1/mu-1)
  return(beta)
}

##Generate stochastic parameters for probabilistic sensitivity analysis
p.stable_moderate <- rbeta(nsim,estAlpha(mu.stable_moderate,mu.stable_moderate*var.p),estBeta(mu.stable_moderate,mu.stable_moderate*var.p))
p.stable_severe   <- rbeta(nsim,estAlpha(mu.stable_severe,mu.stable_severe*var.p),estBeta(mu.stable_severe,mu.stable_severe*var.p))
p.moderate_stable <- rbeta(nsim,estAlpha(mu.moderate_stable,mu.moderate_stable*var.p),estBeta(mu.moderate_stable,mu.moderate_stable*var.p))
p.moderate_severe <- rbeta(nsim,estAlpha(mu.moderate_severe,mu.moderate_severe*var.p),estBeta(mu.moderate_severe,mu.moderate_severe*var.p))
p.severe_stable   <- rbeta(nsim,estAlpha(mu.severe_stable,mu.severe_stable*var.p),estBeta(mu.severe_stable,mu.severe_stable*var.p))
p.severe_moderate <- rbeta(nsim,estAlpha(mu.severe_moderate,mu.severe_moderate*var.p),estBeta(mu.severe_moderate,mu.severe_moderate*var.p))

##For Toy demonstration, mortality is a time-dependent parameter (probability changes over time). 
mu.stable_dead    <- 0.004
mu.moderate_dead  <- 0.005
mu.severe_dead    <- 0.04
p.stable_dead <- p.moderate_dead <- p.severe_dead <- matrix(0,nrow=n.t,ncol=nsim)
p.stable_dead[1,]    <- rbeta(nsim,estAlpha(mu.stable_dead,mu.stable_dead*var.p),estBeta(mu.stable_dead,mu.stable_dead*var.p))
p.moderate_dead[1,]  <- rbeta(nsim,estAlpha(mu.moderate_dead,mu.moderate_dead*var.p),estBeta(mu.moderate_dead,mu.moderate_dead*var.p))
p.severe_dead[1,]    <- rbeta(nsim,estAlpha(mu.severe_dead,mu.severe_dead*var.p),estBeta(mu.severe_dead,mu.severe_dead*var.p))
for (t in 2:n.t) {
  p.stable_dead[t,]   <- p.stable_dead[t-1,]*(1.016^(0.024*t))
  p.moderate_dead[t,] <- p.moderate_dead[t-1,]*(1.016^(0.024*t))
  p.severe_dead[t,]   <- p.severe_dead[t-1,]*(1.016^(0.024*t))
}

#Treatment Effect 
#Submitter treatment reports an average 20% reduction in progression from moderate state to severe state
mu.tx.moderate_severe  <- 0.20
var.tx.moderate_severe <- 0.02
tx.moderate_severe     <- rnorm(nsim,mu.tx.moderate_severe,var.tx.moderate_severe)
#tx.moderate_severe     <- rbeta(nsim,estAlpha(mu.tx.moderate_severe,var.tx.moderate_severe),estBeta(mu.tx.moderate_severe,var.tx.moderate_severe))
hist(tx.moderate_severe)

#Cost Input Parameters
##Costs per cycle per health state
##Toy simplification; costs by health state are fixed, which is sometimes the case
c.stable   <- 50
c.moderate <- 150
c.severe   <- 200
c.dead     <- 0
c.tx       <- 100  #Submitter cost for all Moderate patients per cycle to avert Severe outcome

#Utility (QALY) Input Parameters
mu.u.stable   <- 0.8
mu.u.moderate <- 0.5
mu.u.severe   <- 0.25
u.dead        <- rep(0,nsim) #QALY value for dead always 0
var.u         <- 0.02 #Toy simplification; all utility parameters have same variance

u.stable   <- rbeta(nsim,estAlpha(mu.u.stable,mu.u.stable*var.u),estBeta(mu.u.stable,mu.u.stable*var.u))  
u.moderate <- rbeta(nsim,estAlpha(mu.u.moderate,mu.u.moderate*var.u),estBeta(mu.u.moderate,mu.u.moderate*var.u))  
u.severe   <- rbeta(nsim,estAlpha(mu.u.severe,mu.u.severe*var.u),estBeta(mu.u.severe,mu.u.severe*var.u))  

####MODEL BUILD####
#Create the cohort trace; Array = {Time, States, Run}
##Create trace for each Strategy
m.M.soc <- m.M.sub <- array(0,
                            dim = c(n.t+1, n.state, nsim),
                            dimnames = list(0:n.t, state, 1:nsim))
m.M.soc[1,,] <- m.M.sub[1,,] <- t0_dist    # set initial state distribution

#Transition Matrix 
## Basic matrix is 2D: {State x State}
## To allow time-dependent probabilities add a dimension: {State x State x Time}
## To allow PSA stochastic probabilities add a dimension: {State x State x Runs} 
## To allow both we're employing a 4D array = {State x State x Time x Runs}
m.P.soc <- m.P.sub <- array(0,
                            dim = c(n.state, n.state, n.t, nsim),
                            dimnames = list(state,state,1:n.t,1:nsim))

#Build transition matrix for each strategy
##In this Toy, most of the transition matrix is identical but we will lay out all separately for clarity
#For Standard of Care
##from Stable
m.P.soc["Stable","Stable",,]     <- 1-(p.stable_moderate+p.stable_severe+p.stable_dead) 
m.P.soc["Stable","Moderate",,]   <- p.stable_moderate
m.P.soc["Stable","Severe",,]     <- p.stable_severe
m.P.soc["Stable","Dead",,]       <- p.stable_dead
##from Moderate
m.P.soc["Moderate","Moderate",,] <- 1-(p.moderate_stable+p.moderate_severe+p.moderate_dead) 
m.P.soc["Moderate","Stable",,]   <- p.moderate_stable
m.P.soc["Moderate","Severe",,]   <- p.moderate_severe
m.P.soc["Moderate","Dead",,]     <- p.moderate_dead
##from severe
m.P.soc["Severe","Severe",,]     <- 1-(p.severe_moderate+p.severe_stable+p.severe_dead) 
m.P.soc["Severe","Moderate",,]   <- p.severe_moderate
m.P.soc["Severe","Stable",,]     <- p.severe_stable
m.P.soc["Severe","Dead",,]       <- p.severe_dead
##from Dead
m.P.soc["Dead","Dead",,]         <- 1

#For Submitter Tx
##from Stable
m.P.sub["Stable","Stable",,]     <- 1-(p.stable_moderate+p.stable_severe+p.stable_dead) 
m.P.sub["Stable","Moderate",,]   <- p.stable_moderate
m.P.sub["Stable","Severe",,]     <- p.stable_severe
m.P.sub["Stable","Dead",,]       <- p.stable_dead
##from Moderate
m.P.sub["Moderate","Moderate",,] <- 1-(p.moderate_stable+(p.moderate_severe*tx.moderate_severe)+p.moderate_dead) 
m.P.sub["Moderate","Stable",,]   <- p.moderate_stable
m.P.sub["Moderate","Severe",,]   <- p.moderate_severe*tx.moderate_severe
m.P.sub["Moderate","Dead",,]     <- p.moderate_dead
##from severe
m.P.sub["Severe","Severe",,]     <- 1-(p.severe_moderate+p.severe_stable+p.severe_dead) 
m.P.sub["Severe","Moderate",,]   <- p.severe_moderate
m.P.sub["Severe","Stable",,]     <- p.severe_stable
m.P.sub["Severe","Dead",,]       <- p.severe_dead
##from Dead
m.P.sub["Dead","Dead",,]         <- 1

####MODEL RUN####
for (n in 1:nsim) {                                      # loop through the number of simulation runs
  for (t in 1:n.t)  {                                    # loop through the number of cycles
    m.M.soc[t+1,,n] <- m.M.soc[t,,n] %*% m.P.soc[,,t,n]  # estimate the state vector for the next cycle (t+1)
  }
  for (t in 1:n.t){                                      
    m.M.sub[t+1,,n] <- m.M.sub[t,,n] %*% m.P.sub[,,t,n]    
  }
}  

#Generate summary results from PSA outcomes
d.M.soc <- d.M.sub <- matrix(0, nrow = n.t, ncol = n.state, dimnames = list(1:n.t,state))
for (s in 1:n.state) {                                      
  for (t in 0:n.t)  { 
    d.M.soc[t,s] <- mean(m.M.soc[t,s,])
    d.M.sub[t,s] <- mean(m.M.sub[t,s,])
  }
}

#sum(rowSums(d.M.soc)) == n.t ##Check, TRUE means Markov did not drop persons
#sum(rowSums(d.M.sub)) == n.t ##Check, TRUE means Markov did not drop persons

####PLOT RESULTS####
matplot(d.M.soc[1:n.t,], type = 'l',
        xlab = "Cycle",
        ylab = "Distribution across health states",
        ylim = c(0,1),
        main = "Cohort Trace: Standard of Care",
        col= c("springgreen4","royalblue1","firebrick4","black"),
        lwd = 3)
legend("topright", state, col= c("springgreen4","royalblue1","firebrick4","black"),lty = 1:4, bty = "n")
#Plot vertical line where Severe health state is highest
abline(v = which.max(d.M.soc[,"Severe"]), col = "azure4")

matplot(d.M.sub[1:n.t,], type = 'l',
        xlab = "Cycle",
        ylab = "Distribution across health states",
        ylim = c(0,1),
        main = "Cohort Trace: Submitter's Intervention",
        col= c("springgreen4","royalblue1","firebrick4","black"),
        lwd = 3)
legend("topright", state, col= c("springgreen4","royalblue1","firebrick4","black"),lty = 1:4, bty = "n")
#Plot vertical line where Severe health state is highest
abline(v = which.max(d.M.sub[,"Severe"]), col = "azure4")

#Generating 'Base Case' Aggregated PSA CEA Results 
#Mean Costs & QALYs
v.tc.soc <- d.M.soc %*% c(c.stable,c.moderate,c.severe,c.dead)
v.tc.sub <- d.M.sub %*% c(c.stable,c.moderate+c.tx,c.severe,c.dead)
v.tu.soc <- d.M.soc %*% c(mean(u.stable),mean(u.moderate),mean(u.severe),mean(u.dead))
v.tu.sub <- d.M.sub %*% c(mean(u.stable),mean(u.moderate),mean(u.severe),mean(u.dead))

#Discounted Mean Costs & QALYs
tc_soc <- v.tc.soc * disc.c
tc_sub <- v.tc.sub * disc.c
tu_soc <- v.tu.soc * disc.u
tu_sub <- v.tu.sub * disc.u

#Store Each Strategy into cost and QALY vectors
v.tc <- c(sum(tc_soc),sum(tc_sub))
v.tu <- c(sum(tu_soc),sum(tu_sub))

df_cu <- data.frame(Strategy = strategy,
                    Cost = v.tc,
                    Utility = v.tu/c.t,
                    Inc_Cost = 0,
                    Inc_Utility = 0,
                    INB = 0)

df_cu$Inc_Cost[df_cu$Strategy=="Submitted Tx"]    <- df_cu$Cost[df_cu$Strategy=="Submitted Tx"] - df_cu$Cost[df_cu$Strategy=="Standard of Care"]
df_cu$Inc_Utility[df_cu$Strategy=="Submitted Tx"] <- df_cu$Utility[df_cu$Strategy=="Submitted Tx"] - df_cu$Utility[df_cu$Strategy=="Standard of Care"]
df_cu$INB[df_cu$Strategy=="Submitted Tx"]         <- df_cu$Inc_Utility[df_cu$Strategy=="Submitted Tx"]*100000 - df_cu$Inc_Cost[df_cu$Strategy=="Submitted Tx"]

df_cu

#PSA Costs and QALYs; keeping results in matrix for plottting and probabilistic analysis
## Generate results tables
m.tc.soc <- m.tc.sub <- m.tu.soc <- m.tu.sub <- matrix(0, nrow = n.t+1, ncol = nsim, dimnames = list(0:n.t,1:nsim))
for (n in 1:nsim) {
  m.tc.soc[,n] <- m.M.soc[,,n] %*% c(c.stable,c.moderate,c.severe,c.dead)
  m.tc.sub[,n] <- m.M.sub[,,n] %*% c(c.stable,c.moderate+c.tx,c.severe,c.dead)
  m.tu.soc[,n] <- m.M.soc[,,n] %*% c(u.stable[n],u.moderate[n],u.severe[n],u.dead[n])
  m.tu.sub[,n] <- m.M.sub[,,n] %*% c(u.stable[n],u.moderate[n],u.severe[n],u.dead[n])
}
##Discounting Costs and QALYs
for (t in 1:n.t) {
  m.tc.soc[t+1,] <- m.tc.soc[t+1,] * disc.c[t]
  m.tc.sub[t+1,] <- m.tc.sub[t+1,] * disc.c[t]
  m.tu.soc[t+1,] <- m.tu.soc[t+1,] * disc.c[t]
  m.tu.sub[t+1,] <- m.tu.sub[t+1,] * disc.c[t]
}

## Generate PSA table. This will be all important to everything to follow!
params   <- c("Cost","QALY")
n.params <- length(params)

psa <- array(0,
             dim = c(n.strat,nsim,n.params),
             dimnames = list(strategy, 1:nsim, params))

psa["Standard of Care",,"Cost"] <- colSums(m.tc.soc)
psa["Standard of Care",,"QALY"] <- colSums(m.tu.soc)/c.t
psa["Submitted Tx",,"Cost"] <- colSums(m.tc.sub)
psa["Submitted Tx",,"QALY"] <- colSums(m.tu.sub)/c.t

##If you don't want a 3D array, you can create a results data set for each Strategy instead
#psa_soc <- psa_sub <- array(0,
#                            dim = c(nsim,n.params),
#                            dimnames = list(1:nsim, params))
#psa_soc[,"Cost"] <- colSums(m.tc.soc)
#psa_soc[,"QALY"] <- colSums(m.tu.soc)/c.t
#psa_sub[,"Cost"] <- colSums(m.tc.sub)
#psa_sub[,"QALY"] <- colSums(m.tu.sub)/c.t

## Determine what the optimal option is for all WTP thresholds
## 'Optimal' = The highest probability cost-effective strategy
v.cea <- c("NB","Optimal")
cea <- array(0,
             dim = c(length(v.cea),n.strat,nsim,n.wtp),
             dimnames = list(v.cea,strategy,1:nsim,wtp))
for (i in 1:n.wtp) {       ##Determine the net benefit of every strategy, for each WTP threshold
  cea["NB","Standard of Care",,i] <- psa["Standard of Care",,"QALY"] * wtp[i] - psa["Standard of Care",,"Cost"]
  cea["NB","Submitted Tx",,i] <- psa["Submitted Tx",,"QALY"] * wtp[i] - psa["Submitted Tx",,"Cost"]
}
for (n in 1:nsim) {        ##Set the Optimal Strategy's NB in its own valiable so we can refer to it going forward.
  for (i in 1:n.wtp) {
    if (cea["NB","Submitted Tx",n,i] > cea["NB","Standard of Care",n,i]) {
      cea["Optimal","Submitted Tx",n,i] <- 1
    } else {
      cea["Optimal","Standard of Care",n,i] <- 1
    }
  }
}
1 == mean(cea["Optimal","Submitted Tx",,]) + mean(cea["Optimal","Standard of Care",,]) #Check

## Generate Cost Effectiveness Acceptablility Curve
ceac <- matrix(0, 
               nrow = n.wtp, ncol = n.strat,
               dimnames = list(wtp,strategy))
for (i in 1:n.wtp) {
  ceac[i,"Standard of Care"] <- mean(cea["Optimal","Standard of Care",,i])
  ceac[i,"Submitted Tx"]     <- mean(cea["Optimal","Submitted Tx",,i])
}

matplot(ceac[4:n.wtp,], type = 'l',
        xlab = "Willingness to Pay Threshold per QALY gained (,000 $)",
        ylab = "Probability Strategy Cost-Effective",
        main = "Cost Effectivenes Acceptability Curve",
        col= c("black","blue3"),
        lwd = 3)
legend("right", strategy, col= c("black","blue3"),lty = 1:4, bty = "n")
#Plot vertical line at closest WTP value at which strategies are 50/50
#abline(v = which.min(abs(ceac[,2]-ceac[,1])), col = "azure4")

#Submitter *Net* Monetary Benefit Plotting
sub_inc_cost <- psa["Submitted Tx",,"Cost"]-psa["Standard of Care",,"Cost"]
sub_inc_qaly <- psa["Submitted Tx",,"QALY"]-psa["Standard of Care",,"QALY"]

sub_nmb <- matrix(0,
                  nrow = nsim, ncol = n.wtp,
                  dimnames = list(1:nsim,wtp))

for (i in 1:n.wtp) {
  sub_nmb[,i] <- sub_inc_qaly * wtp[i] - sub_inc_cost
}

##Plot some NMB and NHB results - Not especially interesting visually but good reference
matplot(colMeans(sub_nmb), type = 'l',
        xlab = "Willingness to Pay Threshold",
        ylab = "Net Monetary Benefit",
        main = "NMB, by willingness to pay threshold",
        col= c("blue3"),
        lwd = 3)

##We'll get momentarily distracted and look at the PSA results for NMB at given WTP thresholds
df_sub_nmb <- data.frame(sub_nmb)
#colnames(df_nmb) <- wtp.str
nmb.wtp.50k <- ggplot(df_sub_nmb,aes(x = X50000))
nmb.wtp.50k +  geom_density(aes(color = X50000), size = 1, fill="gray") + 
  geom_vline(aes(xintercept = 0),linetype = "dashed", size = 0.8) +
  ggtitle("Monetary benefit of 'Submitter Tech'","WTP = 50,000") +
  xlab("Net Monetary Benefit")

nmb.wtp.100k <- ggplot(df_sub_nmb,aes(x = X100000))
nmb.wtp.100k +  geom_density(aes(color = X100000), size = 1, fill="gray") + 
  geom_vline(aes(xintercept = 0),linetype = "dashed", size = 0.8) +
  ggtitle("Monetary benefit of 'Submitter Tech'","WTP = 100,000") +
  xlab("Net Monetary Benefit")

nmb.wtp.150k <- ggplot(df_sub_nmb,aes(x = X150000))
nmb.wtp.150k +  geom_density(aes(color = X150000), size = 1, fill="gray") + 
  geom_vline(aes(xintercept = 0),linetype = "dashed", size = 0.8) +
  ggtitle("Monetary benefit of 'Submitter Tech'","WTP = 150,000") +
  xlab("Net Monetary Benefit")

##Clean Up & Save Workspace
run.time  <- Sys.time()-clock
print(run.time)

save.image(file = "IHEworkshop2023ToyMarkov.RData")
