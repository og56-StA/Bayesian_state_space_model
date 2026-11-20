cat("\f")

library(statsecol)

library(jagsUI)
library(coda)

library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)

data("wildebeest")

#Path to BUGS model
modelPath <- "BUGS_models/N_logNorm_y_norm.txt" #Other model is N_logNorm_y_logNorm.txt

sigma2 <- wildebeest$sehat^2
tauObs <- 1/sigma2

obsIdx <- which(!is.na(wildebeest$Nhat))

#Add data to 1994 - RC
nproj <- 5
yrs   <- c(wildebeest$year, max(wildebeest$year) + 1:nproj)      #1960-1994

#Center rainfall around mean, so that B0 is informative - RC
rainMean <- mean(wildebeest$rain)
rain     <- c(wildebeest$rain, rep(rainMean, nproj)) - rainMean

#Harvest assumption that removals continue at current rate - RC
removals <- c(wildebeest$Catch, rep(tail(wildebeest$Catch, 1), nproj))

#No observations in projection years, fill-in precision is arbitrary and ignored - RC
y      <- c(wildebeest$Nhat, rep(NA, nproj))
tauObs <- c(tauObs,          rep(1,  nproj))

#Model data (changed to match updated vars - RC)
data <- list(n        = length(yrs),
             y        = y,
             rain     = rain,
             removals = removals,
             tauObs   = tauObs,
             obsIdx   = obsIdx)

#Initial values for MCMC chains (changed so that chains start from different positions - RC)
chainInits <- function(){
  list(beta0    = runif(1, -0.2, 0.3),
       beta1    = runif(1, -0.2, 0.2),
       tauState = runif(1, 4, 400),
       N        = c(runif(1, 0.1, 0.5), rep(NA, length(yrs) - 1)))
}

#Parameters to monitor (added beto0,beta1,lambda,lambdaAvf,pBetaPos - RC)
monitors <- c("N", "lambda", "beta0", "beta1", "tauState", "lambdaAvg", "pBeta1Pos")

#MCMC parameters
niter   <- 300000 # upped to bump ESS - RC
nburnin <- 10000 # upped - RC
thin    <- 1
nchains <- 3

#JAGS model
model <- jags(data  = data,
              inits = chainInits,
              parameters.to.save = monitors,
              model.file = modelPath,
              n.chains = nchains,
              n.iter = niter,
              n.burnin = nburnin,
              n.thin = thin)

#Diagnostic checks
RHat <- model$Rhat  #All sufficiently low, all Rhat <1.01
ESS  <- model$n.eff #All above 10k 

#Get model samples
samples <- model$samples

#Get samples from each chain
chain1Samples <- samples[[1]]
chain2Samples <- samples[[2]]

#Function to get only stochastic nodes. E.g. not the known observations
getStoch <- function(samples){
  stochNodes <- samples %>%
    as.data.frame() %>%
    select_if(function(col){length(unique(col)) > 1})
}

#Get samples from stochastic nodes for each chain
chain1Stoch <- getStoch(chain1Samples)
chain2Stoch <- getStoch(chain2Samples)

#Combine them
combinedStoch <- mcmc.list(as.mcmc(chain1Stoch),
                           as.mcmc(chain2Stoch))

#Get the summary
stochSummary <- summary(combinedStoch)

#Get time series std. error and posterior std. deviation
timeSeriesSE <- stochSummary$statistics[,4]
stdDev       <- stochSummary$statistics[,2]

#Time series SE / SD. 
MCerrThumb <- timeSeriesSE/stdDev #All less than 0.005. So the changes in the 
#MC estimates of the posterior means due to 
#the MC simulation are negligible

#Convert chains to data frames
chain1Df <- as.data.frame(chain1Samples)
chain2Df <- as.data.frame(chain2Samples)

#Combine them
combinedChainsDf <- rbind(chain1Df,
                          chain2Df) %>%
  as.data.frame()

#Report process s.d - easier to interpret than tau 
sigmaDraws <- 1 / sqrt(combinedChainsDf$tauState)
quantile(sigmaDraws, c(0.025, 0.5, 0.975))

#Parameter summary (added to check predictions for brief - RC)
paramSummary <- combinedChainsDf %>%
  select(beta0, beta1, lambdaAvg, pBeta1Pos, "N[35]") %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(mean  = mean(value),
            sd    = sd(value),
            lower = quantile(value, 0.025),
            upper = quantile(value, 0.975))
paramSummary
#pBeta1Pos mean = posterior probability that rainfall increases growth

#Function to summarise parameters and return data in long format
getSummary <- function(var){
  varDf <- combinedChainsDf %>%
    select(starts_with(var))
  
  names(varDf) <- 1:ncol(varDf)
  
  varDf <- varDf %>%
    pivot_longer(everything())
  
  summ <- varDf %>%
    group_by(name) %>%
    summarise(mean = mean(value),
              sd   = sd(value),
              lower = quantile(value, 0.025),
              upper = quantile(value, 0.975)) %>%
    mutate(name = as.numeric(name)) %>%
    arrange(name)
  
  names(summ) <- paste0(names(summ), "_", var)
  
  summ <- summ %>%
    mutate(year = yrs) %>% #changed to match added predicted years - RC
    select(-names(summ)[1])
  
  return(summ)
}

#Get state summaries
summState <- getSummary("N")

#Apped summaries to wildebeest data (changed to match new params - RC)
wildebeestSumm <- data.frame(year = yrs) %>%
  left_join(wildebeest, by = "year") %>%
  left_join(summState,  by = "year") %>%
  mutate(observed = !is.na(Nhat))

#Plot populations
ggplot(wildebeestSumm, aes(x = year)) +
  #NEW: shade projection period
  annotate("rect", xmin = 1989.5, xmax = 1994.5, ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "orange") +
  #True state
  geom_line(aes(y = mean_N, colour = "State")) +
  geom_point(aes(y = mean_N, colour = "State")) +
  geom_ribbon(aes(ymin = lower_N, ymax = upper_N, fill = "State"), alpha = 0.25) +
  #Observed data
  geom_point(aes(y = Nhat, colour = "Observation")) +
  geom_errorbar(aes(ymin = lci, ymax = uci, colour = "Observation")) +
  labs(x = "Year",
       y = "Population size (millions)",
       title = "Estimated wildebeest population 1960-1989, projected to 1994") +
  scale_color_manual(name = "Process",
                     values = c("State" = "#000000", "Observation" = "#0261ff"),
                     limits = c("State", "Observation")) +
  scale_fill_manual(name = "Process",
                    values = c("State" = "#000000", "Observation" = NA),
                    limits = c("State", "Observation")) +
  theme_bw()
