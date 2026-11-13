cat("\f")

library(statsecol)

library(jagsUI)
library(coda)

library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)

data("wildebeest")

windowsFonts(Font = windowsFont("CMU Serif"))

#Path to BUGS model
modelPath <- "BUGS_models/N_logNorm_y_norm.txt" #Other model is N_logNorm_y_logNorm.txt

if (modelPath == "BUGS_models/N_logNorm_y_logNorm.txt"){
  cv     <- wildebeest$cvhat
  sigma2 <- log(1 + cv^2)
  tauObs <- 1/sigma2
  
  tauObs[is.na(tauObs)] <- 1
} else if (modelPath == "BUGS_models/N_logNorm_y_norm.txt"){
  sigma2 <- wildebeest$sehat^2
  tauObs <- 1/sigma2
  
  tauObs[is.na(tauObs)] <- 1
}

#Model data
data <- list(n        = nrow(wildebeest),
             y        = wildebeest$Nhat,
             rain     = wildebeest$rain,
             removals = wildebeest$Catch,
             tauObs   = tauObs)

#Initial values for MCMC chains
chainInits <- list(list(N = rep(2, nrow(wildebeest)), beta0 = 0, beta1 = 0, tauState = 1),
                   list(N = rep(2, nrow(wildebeest)), beta0 = 0, beta1 = 0, tauState = 1))

#Parameters to monitor
monitors <- c("N", "y", "lambda", "tauState")

#MCMC parameters
niter   <- 100000
nburnin <- 5000
thin    <- 1
nchains <- 2

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
RHat <- model$Rhat  #All sufficiently low. Highest is 1.014071
ESS  <- model$n.eff #Smallest is 835. I believe this is sufficient, but we can increase niter

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
MCerrThumb <- timeSeriesSE/stdDev #All less than 0.05. So the changes in the 
                                  #MC estimates of the posterior means due to 
                                  #the MC simulation are negligible

#Convert chains to data frames
chain1Df <- as.data.frame(chain1Samples)
chain2Df <- as.data.frame(chain2Samples)

#Combine them
combinedChainsDf <- rbind(chain1Df,
                          chain2Df) %>%
  as.data.frame()

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
    mutate(year = wildebeest$year) %>%
    select(-names(summ)[1])
  
  return(summ)
}

#Get observation and state summaries
summObs   <- getSummary("y")
summState <- getSummary("N")

#Apped summaries to wildebeest data
wildebeestSumm <- wildebeest %>%
  left_join(summObs,   by = "year") %>%
  left_join(summState, by = "year") %>%
  mutate(observed = !is.na(Nhat))

#Use observed values to fill in std. errors and credible intervals
wildebeestSumm$sd_y   [!is.na(wildebeest$Nhat)] <- wildebeest$sehat[!is.na(wildebeest$Nhat)]
wildebeestSumm$lower_y[!is.na(wildebeest$Nhat)] <- wildebeest$lci  [!is.na(wildebeest$Nhat)]
wildebeestSumm$upper_y[!is.na(wildebeest$Nhat)] <- wildebeest$uci  [!is.na(wildebeest$Nhat)]

wildebeestSumm$mean_y [is.na(wildebeest$Nhat)] <- NA
wildebeestSumm$sd_y   [is.na(wildebeest$Nhat)] <- NA
wildebeestSumm$lower_y[is.na(wildebeest$Nhat)] <- NA
wildebeestSumm$upper_y[is.na(wildebeest$Nhat)] <- NA

#Plot populations
ggplot(wildebeestSumm, aes(x = year)) +
  #True state data
  geom_line(aes(y = mean_N, colour = "State")) +
  geom_point(aes(y = mean_N, colour = "State")) +
  geom_ribbon(aes(ymin = lower_N, ymax = upper_N, fill = "State"), alpha = 0.25) +
  #Observed data
  geom_point(aes(y = mean_y, colour = "Observation")) +
  geom_errorbar(aes(ymin = lower_y, ymax = upper_y, colour = "Observation")) +
  #Formatting
  labs(x = "Year",
       y = "Estimate",
       title = "Estimate of population and observations over study period") +
  scale_color_manual(name = "Process",
                     values = c("State"       = "#000000",
                                "Observation" = "#0261ff"),
                     limits = c("State",
                                "Observation")) +
  scale_fill_manual(name = "Process",
                    values = c("State"       = "#000000",
                               "Observation" = NA),
                    limits = c("State",
                               "Observation")) +
  theme_bw() +
  theme(text = element_text(family = "Font"))
