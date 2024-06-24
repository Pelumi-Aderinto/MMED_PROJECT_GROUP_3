# Load necessary libraries
library(boot)
library(deSolve)
# library(ellipse)
library(ggplot2)

# Function to set disease parameters with default values
disease_params <- function(lambda = 0.9, a = 5.5, cMax = 0.85, cRate = 1, cHalf = 1990, delta = 0.1, b = 0.02, mu = 0.01) {
  return(as.list(environment()))
}

original_paramters <- c(lambda = 0.9, a = 5.5, cMax = 0.85, cRate = 1, cHalf = 1990, delta = 0.1, b = 0.02, mu = 0.01)

# Define the SI ODE model
SImod <- function(tt, yy, parms) with(as.list(c(parms, yy)), {
  I <- I1 + I2 + I3 + I4  # Total infected
  N <- S + I  # Total population
  lambdaHat <- lambda * exp(-a * I / N)  # Effective transmission rate
  C_t <- 1 - cMax / (1 + exp(-cRate * (tt - cHalf)))  # Control function C(t)
  FOI <- lambdaHat * I / N * C_t  # Force of infection
  g <- 4 * delta  # Transition rate between stages
  
  deriv <- rep(NA, 7)
  deriv[1] <- b * N - FOI * S - mu * S  # Susceptibles
  deriv[2] <- FOI * S - g * I1 - mu * I1  # Infection class I1
  deriv[3] <- g * I1 - g * I2 - mu * I2  # Infection class I2
  deriv[4] <- g * I2 - g * I3 - mu * I3  # Infection class I3
  deriv[5] <- g * I3 - g * I4 - mu * I4  # Infection class I4
  deriv[6] <- FOI * S  # Cumulative incidence
  deriv[7] <- g * I4  # Cumulative mortality
  
  return(list(deriv))
})

# Function to run the deterministic model simulation
simEpidemic <- function(init, tseq, modFunction, parms) {
  simDat <- as.data.frame(ode(y = init, times = tseq, func = modFunction, parms = parms))
  simDat$I <- rowSums(simDat[, paste0("I", 1:4)])
  simDat$N <- rowSums(simDat[, c("S", paste0("I", 1:4))])
  simDat$P <- simDat$I / simDat$N
  return(simDat)
}

# Initial conditions and time sequence
# initPrev <- exp(-7)  # Infected at start
init <- c(S = 999, I1 = 5, I2 = 3, I3 = 1, I4 = 1, CI = 0, CD = 0)
tseqMonth <- seq(1980, 2020, by = 1)

# Simulate the epidemic with default parameters
trueParms <- disease_params()
simDat <- simEpidemic(init, tseq = tseqMonth, modFunction = SImod, parms = trueParms)

# Plot simulated prevalence through time
ggplot(simDat, aes(x = time, y = P)) +
  geom_line(color = "red") +
  labs(title = "Simulated HIV Prevalence Over Time",
       x = "Year",
       y = "Prevalence") +
  theme_minimal()



# Function to sample the epidemic data
sampleEpidemic <- function(simDat, sampleDates = seq(1980, 2020, by = 3), numSamp = rep(200, length(sampleDates))) {
  prev_at_sample_times <- simDat[simDat$time %in% sampleDates, "P"]
  numPos <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
  lci <- mapply(function(x, n) binom.test(x, n)$conf.int[1], x = numPos, n = numSamp)
  uci <- mapply(function(x, n) binom.test(x, n)$conf.int[2], x = numPos, n = numSamp)
  return(data.frame(time = sampleDates, numPos, numSamp, sampPrev = numPos / numSamp, lci = lci, uci = uci))
}

# Sample the epidemic data
set.seed(1)
myDat <- sampleEpidemic(simDat)

# Plot sampled prevalence with confidence intervals
ggplot(myDat, aes(x = time, y = sampPrev)) +
  geom_point(color = "red", size = 2) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2, color = "red") +
  labs(title = "Sampled HIV Prevalence with Confidence Intervals",
       x = "Year",
       y = "Sampled Prevalence") +
  theme_minimal() +
  geom_line(data = simDat, aes(x = time, y = P), color = "blue", linetype = "dashed")


# Negative log-likelihood function
nllikelihood <- function(parms = disease_params(), obsDat = myDat) {
  simDat <- simEpidemic(init, tseq = tseqMonth, modFunction = SImod, parms = parms)
  matchedTimes <- simDat$time %in% obsDat$time
  nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P[matchedTimes], log = TRUE)
  return(sum(nlls))
}

# Calculate log-likelihood for true parameters
nllikelihood(trueParms)

# Calculate log-likelihood for some random guessed parameters
nllikelihood(disease_params(lambda = 0.1, a = 0.05))


# Function to substitute parameters
subsParms <- function(fit.params, fixed.params = disease_params()) {
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl("log_", names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl("log_", names(fit.params))]
    for (nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for (nm in loggedParms) assign(gsub("log_", "", nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })
}

# Make likelihood a function of fixed and fitted parameters
objFXN <- function(fit.params, fixed.params = disease_params(), obsDat = myDat) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat)
}

# Initial parameter guesses
init.pars <- c(log_lambda = log(0.5), log_a = log(3), log_cMax = log(0.65), log_cRate = log(2), log_cHalf = log(2000), log_delta = log(0.02), log_b = log(0.03), log_mu = log(0.02))

# Optimize using SANN
optim.vals <- optim(par = init.pars, objFXN, fixed.params = disease_params(), obsDat = myDat, control = list(trace = 3, maxit = 150), method = "SANN")
exp(unname(optim.vals$par))

# Optimize using Nelder-Mead
optim.vals <- optim(par = optim.vals$par, objFXN, fixed.params = disease_params(), obsDat = myDat, control = list(trace = 3, maxit = 800, reltol = 10^-7), method = "Nelder-Mead", hessian = TRUE)
optim.vals
MLEfits <- optim.vals$par
final_parameters <- c(exp(unname(MLEfits)))


# Plot MLE fit time series
fitDat <- simEpidemic(init, tseq = tseqMonth, modFunction = SImod, parms = subsParms(optim.vals$par, trueParms))

ggplot() +
  geom_line(data = simDat, aes(x = time, y = P), color = "red") +
  geom_line(data = fitDat, aes(x = time, y = P), color = "blue") +
  geom_point(data = myDat, aes(x = time, y = sampPrev), color = "red", size = 2) +
  geom_errorbar(data = myDat, aes(x = time, ymin = lci, ymax = uci), width = 0.2, color = "red") +
  labs(title = "True and Fitted HIV Prevalence",
       x = "Year",
       y = "Prevalence") +
  theme_minimal()


# # Calculate Fisher Information Matrix
# fisherInfMatrix <- solve(optim.vals$hessian)
# 
# # Initialize plot of parameters
# plot(1, 1, type = 'n', log = 'xy', xlim = c(0.05, 1), ylim = c(0.01, 1),
#      xlab = expression(lambda), ylab = expression(a), main = "-log(likelihood) contours", bty = "n")
# 
# # Add true parameter values to the plot
# points(exp(MLEfits['log_lambda']), exp(MLEfits['log_a']), pch = 16, cex = 2, col = 'red')
# # Add MLE to the plot
# points(exp(MLEfits['log_lambda']), exp(MLEfits['log_a']), pch = 16, cex = 2, col = 'black')
# # Add 95% contour ellipse from Hessian
# lines(exp(ellipse(fisherInfMatrix, centre = MLEfits, level = .95)))
# 
# legend("topleft", c('truth', 'MLE', '95% Confidence Region'), lty = c(NA, NA, 1), pch = c(16, 16, NA),
#        col = c('red', 'black', 'black'), bty = 'n')
