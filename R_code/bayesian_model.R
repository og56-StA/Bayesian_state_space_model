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
  list(beta0      = runif(1, -0.2, 0.3),
       beta1      = runif(1, -0.2, 0.2),
       sigmaState = runif(1, 0.05, 0.5),
       N          = c(runif(1, 0.1, 0.5), rep(NA, length(yrs) - 1)))
}

#Parameters to monitor (added beto0,beta1,lambda,lambdaAvf,pBetaPos - RC)
monitors <- c("N", "lambda", "beta0", "beta1", "sigmaState", "tauState",
              "lambdaAvg", "pBeta1Pos")
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
combinedStoch    <- mcmc.list(lapply(samples, function(ch) as.mcmc(getStoch(ch))))

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

#Combine them, 3 chains RC
combinedChainsDf <- do.call(rbind, lapply(samples, as.data.frame))

#Report process s.d - easier to interpret than tau, done using sigmaState so tauState not needed RC 
sigmaDraws <- combinedChainsDf$sigmaState
quantile(sigmaDraws, c(0.025, 0.5, 0.975))

#Parameter summary (added to check predictions for brief - RC)
paramSummary <- combinedChainsDf %>%
  select(beta0, beta1, lambdaAvg, pBeta1Pos, "N[35]") %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(mean   = mean(value),
            median = median(value),     
            sd     = sd(value),
            lower  = quantile(value, 0.025),
            upper  = quantile(value, 0.975))
paramSummary
#pBeta1Pos mean = posterior probability that rainfall increases growth

#Function to summarise parameters and return data in long format
getSummary <- function(var){
  varDf <- combinedChainsDf %>%
    select(starts_with(var) & ends_with("]")) #ends_with is to only select parameters with an index, i.e. NOT lambdaAvg
  
  names(varDf) <- 1:ncol(varDf)
  
  varDf <- varDf %>%
    pivot_longer(everything())
  
  summ <- varDf %>%
    group_by(name) %>%
    summarise(mean = mean(value),
              median = median(value),
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
summState  <- getSummary("N")
summLambda <- getSummary("lambda")

#Apped summaries to wildebeest data (changed to match new params - RC)
wildebeestSumm <- data.frame(year = yrs) %>%
  left_join(wildebeest, by = "year") %>%
  left_join(summState,  by = "year") %>%
  left_join(summLambda, by = "year") %>%
  mutate(observed = !is.na(Nhat))

theme_report <- theme_classic(base_size = 12) +
  theme(legend.position   = "top",
        legend.title      = element_blank(),
        legend.background = element_blank())

#Subsample draws for heavy calculations (870,000 draws is more than needed)
set.seed(5767)
drawSub <- slice_sample(combinedChainsDf, n = 10000)

# Fig 1. 
fit  <- filter(wildebeestSumm, year <= 1989)
proj <- filter(wildebeestSumm, year >= 1989)

fig1 <- ggplot(wildebeestSumm, aes(x = year)) +
  annotate("rect", xmin = 1989.5, xmax = 1994.5, ymin = 0, ymax = Inf,
           fill = "grey92") +
  annotate("text", x = 1992, y = Inf, label = "Projection",
           vjust = 1.5, size = 3.5, colour = "grey30") +
  geom_ribbon(aes(ymin = lower_N, ymax = upper_N,
                  fill = "Estimated abundance (95% CrI)"), alpha = 0.4) +
  geom_line(data = fit,  aes(y = median_N, colour = "Estimated abundance (95% CrI)"),
            linewidth = 0.8) +
  geom_line(data = proj, aes(y = median_N, colour = "Estimated abundance (95% CrI)"),
            linewidth = 0.8, linetype = "dashed") +
  geom_errorbar(aes(ymin = lci, ymax = uci, colour = "Survey estimate (95% CI)"),
                width = 0.4, na.rm = TRUE) +
  geom_point(aes(y = Nhat, colour = "Survey estimate (95% CI)"),
             size = 2, na.rm = TRUE) +
  scale_colour_manual(values = c("Estimated abundance (95% CrI)" = "black",
                                 "Survey estimate (95% CI)"      = "#0072B2")) +
  scale_fill_manual(values = c("Estimated abundance (95% CrI)" = "grey55")) +
  scale_x_continuous(breaks = seq(1960, 1995, 5), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Year", y = "Wildebeest abundance (millions)") +
  theme_report

ggsave("fig1_trajectory.png", fig1, width = 8, height = 4.5, dpi = 300)

#Fig 2. Growth rate against rainfall
rainGrid <- seq(min(wildebeest$rain), max(wildebeest$rain), length.out = 100)

lamDraws <- exp(outer(drawSub$beta0, rep(1, length(rainGrid))) +
                  outer(drawSub$beta1, rainGrid - rainMean))

lamDf <- data.frame(rain   = rainGrid,
                    median = apply(lamDraws, 2, median),
                    lower  = apply(lamDraws, 2, quantile, 0.025),
                    upper  = apply(lamDraws, 2, quantile, 0.975))

fig2 <- ggplot(lamDf, aes(rain)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
  geom_line(aes(y = median), linewidth = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_rug(data = wildebeest, aes(x = rain), inherit.aes = FALSE,
           sides = "b", alpha = 0.7) +
  annotate("text", x = max(rainGrid), y = Inf, hjust = 1, vjust = 1.5, size = 3.5,
           label = sprintf("Pr(rainfall increases growth) = %.2f",
                           mean(combinedChainsDf$pBeta1Pos))) +
  labs(x = "Dry-season rainfall (cm)",
       y = expression("Annual growth multiplier (" * lambda * ")")) +
  theme_report

ggsave("fig2_growth_rainfall.png", fig2, width = 6.5, height = 4.5, dpi = 300)

#Figure 3, posterior dist.
projDf <- data.frame(
  value = c(combinedChainsDf$`N[30]`, combinedChainsDf$`N[35]`),
  year  = rep(c("1989 (last survey year)", "1994 (projected)"),
              each = nrow(combinedChainsDf)))

projMed <- projDf %>% group_by(year) %>% summarise(med = median(value))

fig3 <- ggplot(projDf, aes(value, fill = year)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = projMed, aes(xintercept = med, colour = year),
             linetype = "dashed", show.legend = FALSE) +
  scale_fill_manual(values = c("#0072B2", "#E69F00")) +
  scale_colour_manual(values = c("#0072B2", "#E69F00")) +
  coord_cartesian(xlim = c(0.6, 4.5)) +
  labs(x = "Wildebeest abundance (millions)", y = "Posterior density") +
  theme_report
ggsave("fig3_projection.png", fig3, width = 6.5, height = 4.5, dpi = 300)

traceDf <- do.call(rbind, lapply(seq_along(samples), function(k) {
  ch <- as.data.frame(samples[[k]])
  keep <- seq(1, nrow(ch), by = 100)          #every 100th draw for plotting
  data.frame(iter  = keep, chain = factor(k),
             beta0 = ch$beta0[keep], beta1 = ch$beta1[keep],
             sigma = ch$sigmaState[keep])
})) %>%
  pivot_longer(c(beta0, beta1, sigma))

fig4 <- ggplot(traceDf, aes(iter, value, colour = chain)) +
  geom_line(linewidth = 0.2, alpha = 0.7) +
  facet_wrap(~ name, ncol = 1, scales = "free_y",
             labeller = as_labeller(c(beta0 = "beta[0]", beta1 = "beta[1]",
                                      sigma = "sigma[proc]"), label_parsed)) +
  scale_colour_manual(values = c("#0072B2", "#E69F00", "#009E73")) +
  labs(x = "Iteration (post burn-in)", y = NULL, colour = "Chain") +
  theme_report +
  theme(legend.title = element_text())

ggsave("fig4_traces.png", fig4, width = 7, height = 6, dpi = 300)

#Matrix of N draws (rows = draws, columns = years 1960-1994)
Nmat <- as.matrix(drawSub[, paste0("N[", seq_along(yrs), "]")])

#Realised growth N[t]/N[t-1] and harvest as % of previous year's population
growth  <- Nmat[, -1] / Nmat[, -ncol(Nmat)]
harvest <- sweep(1 / Nmat[, -ncol(Nmat)], 2, removals[-1], "*") * 100

#Summarise each year (median and 95% CrI)
summQ <- function(m) data.frame(year   = yrs[-1],
                                median = apply(m, 2, median),
                                lower  = apply(m, 2, quantile, 0.025),
                                upper  = apply(m, 2, quantile, 0.975))

growthDf <- summQ(growth)  %>% mutate(panel = "Annual growth multiplier")
harvDf   <- summQ(harvest) %>% mutate(panel = "Harvest (% of population)")
plotDf   <- bind_rows(growthDf, harvDf)

#Rainfall-expected growth (lambda) and the no-change reference line, top panel only
lamLine <- summLambda %>%
  filter(year > 1960) %>%
  transmute(year, median = median_lambda, panel = "Annual growth multiplier")
refLine <- data.frame(panel = "Annual growth multiplier", y = 1)

fig5 <- ggplot(plotDf, aes(x = year)) +
  annotate("rect", xmin = 1989.5, xmax = 1994.5, ymin = -Inf, ymax = Inf,
           fill = "grey92") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey55", alpha = 0.4) +
  geom_line(data = growthDf, aes(y = median, colour = "Realised growth (95% CrI)"),
            linewidth = 0.8) +
  geom_line(data = lamLine, aes(y = median, colour = "Expected from rainfall"),
            linewidth = 0.7, linetype = "dashed") +
  geom_line(data = harvDf, aes(y = median), colour = "black", linewidth = 0.8) +
  geom_hline(data = refLine, aes(yintercept = y), linetype = "dotted") +
  facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_colour_manual(values = c("Realised growth (95% CrI)" = "black",
                                 "Expected from rainfall"    = "#D55E00")) +
  scale_x_continuous(breaks = seq(1960, 1995, 5), expand = c(0.01, 0)) +
  labs(x = "Year", y = NULL) +
  theme_report +
  theme(strip.placement  = "outside",
        strip.background = element_blank(),
        strip.text       = element_text(size = 11))

ggsave("fig5_growth_harvest.png", fig5, width = 7, height = 6, dpi = 300)

#Numbers for the text
mean(combinedChainsDf$`N[35]` > combinedChainsDf$`N[30]`)   #Pr(larger in 1994 than 1989)
