## CADTH 2023 Toy # Markov Model ###
## Institute of Health Economics ##

#install.packages("dplyr")
library(dplyr)
2
#install.packages("devtools")
#library(devtools)
#install.packages("ggplot2")
library(ggplot2)

####STRUCTURAL PARAMETERS####
options(scipen = 999)
n.t    <- 60        #Cycles - Time = Months
c.t    <- 12        #Number of Cycles in a Year
disc   <- 0.015/c.t  #Discount Rate
disc.c <- 1/(1 + disc)^(1:(n.t+1)) 
disc.u <- 1/(1 + disc)^(1:(n.t+1))
wtp.def<- 100000
wtp    <- seq(0,150000,by=1000)
wtp.str<- sprintf("%s",seq(0,150000,by=1000))
n.wtp  <- length(wtp)

####INPUT PARAMETERS####
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
p.stable_moderate <- 0.075
p.stable_severe   <- 0.01
p.moderate_stable <- 0.04
p.moderate_severe <- 0.12
p.severe_stable   <- 0.08
p.severe_moderate <- 0.001

##For Toy demonstration, mortality is a time-dependent parameter (probability changes over time). 
p.stable_dead    <- rep(0.004,n.t)
p.moderate_dead  <- rep(0.005,n.t)
p.severe_dead    <- rep(0.04,n.t)
for (t in 2:n.t) {
  p.stable_dead[t]   <- p.stable_dead[t-1]*(1.016^(0.024*t))
  p.moderate_dead[t] <- p.moderate_dead[t-1]*(1.016^(0.024*t))
  p.severe_dead[t]   <- p.severe_dead[t-1]*(1.016^(0.024*t))
}

#Treatment Effect 
#Submitter treatment reports an average 20% reduction in progression from moderate state to severe state
tx.moderate_severe  <- 0.20

#Cost Input Parameters
##Costs per cycle per health state
##Toy simplification; costs by health state are fixed, which is sometimes the case
c.stable   <- 50
c.moderate <- 150
c.severe   <- 200
c.dead     <- 0
c.tx       <- 100  #Submitter cost for all Moderate patients per cycle to avert Severe outcome

#Utility (QALY) Input Parameters
u.stable   <- 0.8
u.moderate <- 0.5
u.severe   <- 0.25
u.dead     <- 0 #QALY value for dead always 0

####MODEL BUILD####
#Create the cohort trace; Array = {Time, Run}
##Create trace for each Strategy
d.M.soc <- d.M.sub <- array(0,
                            dim = c(n.t+1, n.state),
                            dimnames = list(0:n.t, state))
d.M.soc[1,] <- d.M.sub[1,] <- t0_dist    # set initial state distribution

#Transition Matrix 
## Basic matrix is 2D: {State x State}
## To allow time-dependent probabilities add a dimension: {State x State x Time}
m.P.soc <- m.P.sub <- array(0,
                            dim = c(n.state, n.state, n.t),
                            dimnames = list(state,state,1:n.t))

#Build transition matrix for each strategy
##In this Toy, most of the transition matrix is identical but we will lay out all separately for clarity
#For Standard of Care
##from Stable
m.P.soc["Stable","Stable",]     <- 1-(p.stable_moderate+p.stable_severe+p.stable_dead) 
m.P.soc["Stable","Moderate",]   <- p.stable_moderate
m.P.soc["Stable","Severe",]     <- p.stable_severe
m.P.soc["Stable","Dead",]       <- p.stable_dead
##from Moderate
m.P.soc["Moderate","Moderate",] <- 1-(p.moderate_stable+p.moderate_severe+p.moderate_dead) 
m.P.soc["Moderate","Stable",]   <- p.moderate_stable
m.P.soc["Moderate","Severe",]   <- p.moderate_severe
m.P.soc["Moderate","Dead",]     <- p.moderate_dead
##from severe
m.P.soc["Severe","Severe",]     <- 1-(p.severe_moderate+p.severe_stable+p.severe_dead) 
m.P.soc["Severe","Moderate",]   <- p.severe_moderate
m.P.soc["Severe","Stable",]     <- p.severe_stable
m.P.soc["Severe","Dead",]       <- p.severe_dead
##from Dead
m.P.soc["Dead","Dead",]         <- 1

#For Standard of Care
##from Stable
m.P.sub["Stable","Stable",]     <- 1-(p.stable_moderate+p.stable_severe+p.stable_dead) 
m.P.sub["Stable","Moderate",]   <- p.stable_moderate
m.P.sub["Stable","Severe",]     <- p.stable_severe
m.P.sub["Stable","Dead",]       <- p.stable_dead
##from Moderate
m.P.sub["Moderate","Moderate",] <- 1-(p.moderate_stable+(p.moderate_severe*tx.moderate_severe)+p.moderate_dead) 
m.P.sub["Moderate","Stable",]   <- p.moderate_stable
m.P.sub["Moderate","Severe",]   <- p.moderate_severe*tx.moderate_severe
m.P.sub["Moderate","Dead",]     <- p.moderate_dead
##from severe
m.P.sub["Severe","Severe",]     <- 1-(p.severe_moderate+p.severe_stable+p.severe_dead) 
m.P.sub["Severe","Moderate",]   <- p.severe_moderate
m.P.sub["Severe","Stable",]     <- p.severe_stable
m.P.sub["Severe","Dead",]       <- p.severe_dead
##from Dead
m.P.sub["Dead","Dead",]         <- 1

####MODEL RUN####
for (t in 1:n.t)  {                                    # loop through the number of cycles
  d.M.soc[t+1,] <- d.M.soc[t,] %*% m.P.soc[,,t]  # estimate the state vector for the next cycle (t+1)
  d.M.sub[t+1,] <- d.M.sub[t,] %*% m.P.sub[,,t]
}

####PLOT RESULTS####
matplot(d.M.soc[], type = 'l',
        xlab = "Cycle",
        ylab = "Distribution across health states",
        ylim = c(0,1),
        main = "Cohort Trace: Standard of Care",
        col= c("springgreen4","royalblue1","firebrick4","black"),
        lwd = 3)
legend("topright", state, col= c("springgreen4","royalblue1","firebrick4","black"),lty = 1:4, bty = "n")
#Plot vertical line where Severe health state is highest
abline(v = which.max(d.M.soc[,"Severe"]), col = "azure4")

matplot(d.M.sub[], type = 'l',
        xlab = "Cycle",
        ylab = "Distribution across health states",
        ylim = c(0,1),
        main = "Cohort Trace: Submitter's Intervention",
        col= c("springgreen4","royalblue1","firebrick4","black"),
        lwd = 3)
legend("topright", state, col= c("springgreen4","royalblue1","firebrick4","black"),lty = 1:4, bty = "n")
#Plot vertical line where Severe health state is highest
abline(v = which.max(d.M.sub[,"Severe"]), col = "azure4")

#Mean Costs & QALYs
v.tc.soc <- d.M.soc %*% c(c.stable,c.moderate,c.severe,c.dead)
v.tc.sub <- d.M.sub %*% c(c.stable,c.moderate+c.tx,c.severe,c.dead)
v.tu.soc <- d.M.soc %*% c(u.stable,u.moderate,u.severe,u.dead)
v.tu.sub <- d.M.sub %*% c(u.stable,u.moderate,u.severe,u.dead)

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
df_cu$INB[df_cu$Strategy=="Submitted Tx"]         <- df_cu$Inc_Utility[df_cu$Strategy=="Submitted Tx"]*wtp.def - df_cu$Inc_Cost[df_cu$Strategy=="Submitted Tx"]

df_cu


