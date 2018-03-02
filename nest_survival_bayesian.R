# install/load necessary packages
my.packs <- c('circular','rgdal','sp','RColorBrewer','rgeos','spatstat','jagsUI','ggplot2')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

rm(list=ls()) #Clear R's working memory

#Read inp file
inp = read.csv("../../../Data/coei_inp.csv", header=TRUE)

#Column 1: $nestIDyear; nest ID
#Column 2: $FirstFound; Calendar date first found
#Column 3: $LastPresent; Calendar date last seen active (or expected hatch date, whichever is smaller)
#Column 4: $LastChecked; Calendar date last checked (if hatched, this should be equal to last active)
#Column 5: $Fate; Nest Fate (0 = hatched, 1 = failed)
#Column 6: $Year; Year of study

#Sort from earliest to latest year
inp = inp[order(inp$Year),]

inp$FirstFound = as.integer(inp$FirstFound)
inp$LastPresent = as.integer(inp$LastPresent)
inp$LastChecked = as.integer(inp$LastChecked)

#Length of encounter history for each nest
encounter_length = inp$LastChecked - inp$FirstFound + 1
last_active = inp$LastPresent - inp$FirstFound + 1

#Create matrix for nest encounters
#For Bayesian analysis, each nest must be observed on day 1
#And then subsequently fail or hatch
N_mat = matrix(NA, nrow=nrow(inp), ncol = max(encounter_length))

for (i in 1:nrow(inp)){

    N_mat[i,1:last_active[i]] = 1
    N_mat[i,encounter_length[i]] = abs(inp$Fate[i] - 1) #0 if failed, 1 if hatched

}

# JAGS SCRIPT
sink("nest_surv.jags")
cat("model{

    ########
    # Priors
    ########

    for (a in 1:Years){
        phi.intercept[a] ~ dunif(0,1)
    }

    ########
    # Likelihood
    ########
    for (i in 1:N){

        z[i,1]<-1

        for (j in 2:encounter_length[i]){

            #Observed encounter history
            y[i,j] ~ dbern(p[i,j])

            #Status on this day is status previous day * survival rate
            p[i,j] <- (y[i,j-1])*q[i,j]

            #Covariates
            q[i,j]<-phi.intercept[nestYear[i]]


        } #d

    } #i

    }
    ",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = N_mat,
                 N = nrow(N_mat),
                 encounter_length = encounter_length,
                 Years = max(inp$Year) - min(inp$Year) + 1,
                 nestYear = inp$Year - min(inp$Year) + 1)

# Initial values
inits <- function() list()

# Parameters monitored
params <- c("phi.intercept")

#-----------------------------------------------------------------------
#MCMC settings
ni <- 5000 # Number of iterations
nb <- 2000  # Burn-in length
nt <- 5
nc <- 5
#-----------------------------------------------------------------------

# Call JAGS from R (BRT 47 min)
out <- jags(win.data, inits, params, "./nest_surv.jags", n.chains = nc,
          n.thin = nt, n.iter = ni, n.burnin = nb)

#Create dataframe to store model output
annual_NS = data.frame(Year = 1:win.data$Years)

#Nest survival across 28 day incubation is daily survival (phi)^28
annual_NS$phi.50 = apply(out$sims.list$phi.intercept ^ 28,2,function(x) quantile(x,0.5))
annual_NS$phi.025 = apply(out$sims.list$phi.intercept ^ 28,2,function(x) quantile(x,0.025))
annual_NS$phi.975 = apply(out$sims.list$phi.intercept ^ 28,2,function(x) quantile(x,0.975))

#Only select years with actual observations
annual_NS$Year = min(inp$Year):max(inp$Year)
annual_NS = annual_NS[which(annual_NS$Year %in% inp$Year),]

#Print plot
summary_plot = ggplot(data = annual_NS) +
    geom_point(aes(x = Year, y = phi.50))+
    geom_errorbar(aes(x = Year, ymin = phi.025, ymax = phi.975), width = 0.05)+

    ylab("Nest Survival (DSR^28)")+
    theme_bw()

print(summary_plot)

dev.copy(pdf,"nest_surv.pdf", width=7, height=4)
dev.off()
