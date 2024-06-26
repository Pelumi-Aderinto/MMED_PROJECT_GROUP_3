# Load necessary libraries
library(boot)
library(deSolve)
library(ggplot2)

# Function to set disease parameters with default values
disease_params <- function(lambda = 0.9, a = 5.5, cMax = 0.85, cRate = 1, cHalf = 1990, delta = 0.1, b = 0.02, mu = 0.01) {
  return(as.list(environment()))
}

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



# Define the S-I-I-I-I-I-I-I ODE model
SI7mod <- function(tt, yy, parms) with(as.list(c(parms, yy)), {
  I <- I1 + I2 + I3 + I4 + I5 + I6 + I7  # Total infected
  N <- S + I  # Total population
  lambdaHat <- lambda * exp(-a * I / N)  # Effective transmission rate
  C_t <- 1 - cMax / (1 + exp(-cRate * (tt - cHalf)))  # Control function C(t)
  FOI <- lambdaHat * I / N * C_t  # Force of infection
  g <- 7 * delta  # Transition rate between stages
  
  deriv <- rep(NA, 10)
  deriv[1] <- b * N - FOI * S - mu * S  # Susceptibles
  deriv[2] <- FOI * S - g * I1 - mu * I1  # Infection class I1
  deriv[3] <- g * I1 - g * I2 - mu * I2  # Infection class I2
  deriv[4] <- g * I2 - g * I3 - mu * I3  # Infection class I3
  deriv[5] <- g * I3 - g * I4 - mu * I4  # Infection class I4
  deriv[6] <- g * I4 - g * I5 - mu * I5  # Infection class I5
  deriv[7] <- g * I5 - g * I6 - mu * I6  # Infection class I6
  deriv[8] <- g * I6 - g * I7 - mu * I7  # Infection class I7
  deriv[9] <- FOI * S  # Cumulative incidence
  deriv[10] <- g * I7  # Cumulative mortality
  
  return(list(deriv))
})

# Define the S-I ODE model
SImod_simple <- function(tt, yy, parms) with(as.list(c(parms, yy)), {
  I <- I1  # Total infected
  N <- S + I  # Total population
  lambdaHat <- lambda * exp(-a * I / N)  # Effective transmission rate
  C_t <- 1 - cMax / (1 + exp(-cRate * (tt - cHalf)))  # Control function C(t)
  FOI <- lambdaHat * I / N * C_t  # Force of infection
  g <- delta  # Transition rate between stages
  
  deriv <- rep(NA, 4)
  deriv[1] <- b * N - FOI * S - mu * S  # Susceptibles
  deriv[2] <- FOI * S - g * I1 - mu * I1  # Infection class I1
  deriv[3] <- FOI * S  # Cumulative incidence
  deriv[4] <- g * I1  # Cumulative mortality
  
  return(list(deriv))
})


# Define the SI ODE model without intervention
SImod_no_intervention <- function(tt, yy, parms) with(as.list(c(parms, yy)), {
  I <- I1 + I2 + I3 + I4  # Total infected
  N <- S + I  # Total population
  lambdaHat <- lambda * exp(-a * I / N)  # Effective transmission rate
  FOI <- lambdaHat * I / N  # Force of infection without C_t
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

# Define the S-I-I-I-I-I-I-I ODE model without intervention
SI7mod_no_intervention <- function(tt, yy, parms) with(as.list(c(parms, yy)), {
  I <- I1 + I2 + I3 + I4 + I5 + I6 + I7  # Total infected
  N <- S + I  # Total population
  lambdaHat <- lambda * exp(-a * I / N)  # Effective transmission rate
  FOI <- lambdaHat * I / N 
  # FOI <- lambdaHat * I / N  # Force of infection without C_t
  g <- 7 * delta  # Transition rate between stages
  
  deriv <- rep(NA, 10)
  deriv[1] <- b * N - FOI * S - mu * S  # Susceptibles
  deriv[2] <- FOI * S - g * I1 - mu * I1  # Infection class I1
  deriv[3] <- g * I1 - g * I2 - mu * I2  # Infection class I2
  deriv[4] <- g * I2 - g * I3 - mu * I3  # Infection class I3
  deriv[5] <- g * I3 - g * I4 - mu * I4  # Infection class I4
  deriv[6] <- g * I4 - g * I5 - mu * I5  # Infection class I5
  deriv[7] <- g * I5 - g * I6 - mu * I6  # Infection class I6
  deriv[8] <- g * I6 - g * I7 - mu * I7  # Infection class I7
  deriv[9] <- FOI * S  # Cumulative incidence
  deriv[10] <- g * I7  # Cumulative mortality
  
  return(list(deriv))
})

# Define the S-I ODE model without intervention
SImod_simple_no_intervention <- function(tt, yy, parms) with(as.list(c(parms, yy)), {
  I <- I1  # Total infected
  N <- S + I  # Total population
  lambdaHat <- lambda * exp(-a * I / N)  # Effective transmission rate
  FOI <- lambdaHat * I / N
  # FOI <- lambdaHat * I / N  # Force of infection without C_t
  g <- delta  # Transition rate between stages
  
  deriv <- rep(NA, 4)
  deriv[1] <- b * N - FOI * S - mu * S  # Susceptibles
  deriv[2] <- FOI * S - g * I1 - mu * I1  # Infection class I1
  deriv[3] <- FOI * S  # Cumulative incidence
  deriv[4] <- g * I1  # Cumulative mortality
  
  return(list(deriv))
})



# Function to run the deterministic model simulation
simEpidemic <- function(init, tseq, modFunction, parms) {
  # simDat <- as.data.frame(ode(y = init, times = tseq, func = modFunction, parms = parms))
  simDat <- as.data.frame(lsoda(y = init, times = tseq, func = modFunction, parms = parms, rtol = 1e-9, atol = 1e-9))
  
  # Check if there is only one infection compartment
  infection_cols <- grep("^I", names(simDat))
  
  if (length(infection_cols) == 1) {
    simDat$I <- simDat[, infection_cols]
  } else {
    simDat$I <- rowSums(simDat[, infection_cols])
  }
  
  simDat$N <- rowSums(simDat[, c("S", names(simDat)[infection_cols])])
  simDat$P <- simDat$I / simDat$N
  return(simDat)
}


# Initial conditions and time sequence
# initPrev <- exp(-7)  # Infected at start
init_SI7 <- c(S = 990, I1 = 5, I2 = 3, I3 = 1, I4 = 1, I5 = 0, I6 = 0, I7 = 0, CI = 0, CD = 0)
init_SI <- c(S = 990, I1 = 10, CI = 0, CD = 0)
init <- c(S = 990, I1 = 5, I2 = 3, I3 = 1, I4 = 1, CI = 0, CD = 0)
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
sampleEpidemic <- function(simDat, sampleDates = seq(1980, 2020, by = 2), numSamp = rep(200, length(sampleDates))) {
  prev_at_sample_times <- simDat[simDat$time %in% sampleDates, "P"]
  numPos <- rbinom(length(numSamp), numSamp, prev_at_sample_times)
  lci <- mapply(function(x, n) binom.test(x, n)$conf.int[1], x = numPos, n = numSamp)
  uci <- mapply(function(x, n) binom.test(x, n)$conf.int[2], x = numPos, n = numSamp)
  return(data.frame(time = sampleDates, numPos, numSamp, sampPrev = numPos / numSamp, lci = lci, uci = uci))
}


num_simulations <- 100
vector <- c(1, 2, 6, 7, 8)
# Initial parameter guesses
init.pars <- c(log_lambda = log(0.5), log_a = log(3), log_cMax = log(0.65), log_cRate = log(2), log_cHalf = log(2000), log_delta = log(0.02), log_b = log(0.03), log_mu = log(0.02))

# Data frame to store final parameters for each simulation
results_df_SI7 <- data.frame(matrix(ncol = length(init.pars), nrow = num_simulations))
colnames(results_df_SI7) <- names(init.pars)

results_df_SI <- data.frame(matrix(ncol = length(init.pars), nrow = num_simulations))
colnames(results_df_SI) <- names(init.pars)

results_df_SImod <- data.frame(matrix(ncol = length(init.pars), nrow = num_simulations))
colnames(results_df_SImod) <- names(init.pars)

# cum_results  data.frame(matrix(ncol = length(init.pars), nrow = num_simulations))

cum_results <- data.frame((matrix(ncol= 6, nrow = num_simulations)))
colnames(cum_results) <- c("CI_SI4", "CI_SI4_no_inter", "CI_SI7", "CI_SI7_no_inter", "CI_SI", "CI_SI_no_inter")

# Run the optimization in a loop
for (i in 1:num_simulations) {
  # Sample the epidemic data
  # set.seed(i)
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
  nllikelihood <- function(parms = disease_params(), obsDat = myDat, modFunction, init) {
    simDat <- simEpidemic(init, tseq = tseqMonth, modFunction = modFunction, parms = parms)
    matchedTimes <- simDat$time %in% obsDat$time
    nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P[matchedTimes], log = TRUE)
    return(sum(nlls))
  }
  
  
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
  objFXN <- function(fit.params, fixed.params = disease_params(), obsDat = myDat, modFunction, init) {
    parms <- subsParms(fit.params, fixed.params)
    nllikelihood(parms, obsDat = obsDat, modFunction = modFunction, init = init)
  }
  
  
  
  # Optimizing for the S-I-I-I-I-I-I-I model
  optim.vals_SI7 <- optim(par = init.pars, objFXN, fixed.params = disease_params(), obsDat = myDat, modFunction = SI7mod, init = init_SI7, control = list(trace = 3, maxit = 800, reltol = 10^-7), method = "Nelder-Mead", hessian = TRUE)
  final_parameters_SI7 <- c(exp(optim.vals_SI7$par))
  names(final_parameters_SI7) <- c("lambda", "a", "cMax", "cRate", "cHalf", "delta", "b", "mu")
  final_parameters_SI7_no_intervention <- final_parameters_SI7[vector]
  
  # # Optimizing for the S-I model
  optim.vals_SI <- optim(par = init.pars, objFXN, fixed.params = disease_params(), obsDat = myDat, modFunction = SImod_simple, init = init_SI, control = list(trace = 3, maxit = 800, reltol = 10^-7), method = "Nelder-Mead", hessian = TRUE)
  final_parameters_SI <- c(exp(optim.vals_SI$par))
  names(final_parameters_SI) <- c("lambda", "a", "cMax", "cRate", "cHalf", "delta", "b", "mu")
  final_parameters_SI_no_intervention <- final_parameters_SI[vector]
  
  # Optimize using Nelder-Mead
  optim.vals <- optim(par = init.pars, objFXN, fixed.params = disease_params(), obsDat = myDat, modFunction = SImod, init = init, control = list(trace = 3, maxit = 800, reltol = 10^-7), method = "Nelder-Mead", hessian = TRUE)
  final_parameters <- c(exp(optim.vals$par))
  names(final_parameters) <- c("lambda", "a", "cMax", "cRate", "cHalf", "delta", "b", "mu")
  final_parameters_no_intervention <- final_parameters[vector]
  
  # Plot MLE fit time series for SI7 model
  fitDat_SI7 <- simEpidemic(init = init_SI7, tseq = tseqMonth, modFunction = SI7mod, parms = subsParms(optim.vals_SI7$par, trueParms))
  
  # Plot MLE fit time series for SI7 model without intervention
  fitDat_SI7_no_intervention <- simEpidemic(init = init_SI7, tseq = tseqMonth, modFunction = SI7mod_no_intervention, parms = subsParms(final_parameters_SI7_no_intervention, trueParms))
  
  # # Plot MLE fit time series for SI model
  fitDat_SI <- simEpidemic(init = init_SI, tseq = tseqMonth, modFunction = SImod_simple, parms = subsParms(optim.vals_SI$par, trueParms))
  
  # Plot MLE fit time series for SI model without intervention
  fitDat_SI_no_intervention <- simEpidemic(init = init_SI, tseq = tseqMonth, modFunction = SImod_simple_no_intervention, parms = subsParms(final_parameters_SI_no_intervention, trueParms))
  
  
  fitDat <- simEpidemic(init, tseq = tseqMonth, modFunction = SImod, parms = subsParms(optim.vals$par, trueParms))
  
  # Plot MLE fit time series for S-I-I-I-I model without intervention
  fitDat_no_intervention <- simEpidemic(init = init, tseq = tseqMonth, modFunction = SImod_no_intervention, parms = subsParms(final_parameters_no_intervention, trueParms))
  
  
  # Store the final parameters in the results data frame
  results_df_SI7[i, ] <- final_parameters_SI7
  results_df_SI[i, ] <- final_parameters_SI
  results_df_SImod[i, ] <- final_parameters
  cum_results[i, ] <- c(fitDat$CI[nrow(fitDat)], fitDat_no_intervention$CI[nrow(fitDat)], fitDat_SI7$CI[nrow(fitDat_SI7)], fitDat_SI7_no_intervention$CI[nrow(fitDat_SI)], fitDat_SI$CI[nrow(fitDat_SI)], fitDat_SI_no_intervention$CI[nrow(fitDat_SI)])
}

# # Optimizing for the S-I-I-I-I-I-I-I model without intervention
# optim.vals_SI7_no_intervention <- optim(par = init.pars, objFXN, fixed.params = disease_params(), obsDat = myDat, modFunction = SI7mod_no_intervention, init = init_SI7, control = list(trace = 3, maxit = 800, reltol = 10^-7), method = "Nelder-Mead", hessian = TRUE)
# MLEfits_SI7_no_intervention <- optim.vals_SI7_no_intervention$par
# final_parameters_SI7_no_intervention <- c(exp(unname(MLEfits_SI7_no_intervention)))
# 
# # Optimizing for the S-I model without intervention
# optim.vals_SI_no_intervention <- optim(par = init.pars, objFXN, fixed.params = disease_params(), obsDat = myDat, modFunction = SImod_simple_no_intervention, init = init_SI, control = list(trace = 3, maxit = 800, reltol = 10^-7), method = "Nelder-Mead", hessian = TRUE)
# MLEfits_SI_no_intervention <- optim.vals_SI_no_intervention$par
# final_parameters_SI_no_intervention <- c(exp(unname(MLEfits_SI_no_intervention)))
# 
# # Optimizing for the S-I-I-I-I model without intervention
# optim.vals_no_intervention <- optim(par = init.pars, objFXN, fixed.params = disease_params(), obsDat = myDat, modFunction = SImod_no_intervention, init = init, control = list(trace = 3, maxit = 800, reltol = 10^-7), method = "Nelder-Mead", hessian = TRUE)
# MLEfits_no_intervention <- optim.vals_no_intervention$par
# final_parameters_no_intervention <- c(exp(unname(MLEfits_no_intervention)))



# Plot MLE fit time series for SI7 model
fitDat_SI7 <- simEpidemic(init = init_SI7, tseq = tseqMonth, modFunction = SI7mod, parms = subsParms(optim.vals_SI7$par, trueParms))

# Plot MLE fit time series for SI7 model without intervention
fitDat_SI7_no_intervention <- simEpidemic(init = init_SI7, tseq = tseqMonth, modFunction = SI7mod_no_intervention, parms = subsParms(final_parameters_SI7_no_intervention, trueParms))

ggplot() +
  geom_line(data = simDat, aes(x = time, y = P), color = "red") +
  geom_line(data = fitDat_SI7, aes(x = time, y = P), color = "blue") +
  geom_line(data = fitDat_SI7_no_intervention, aes(x = time, y = P), color = "green") +
  geom_point(data = myDat, aes(x = time, y = sampPrev), color = "red", size = 2) +
  geom_errorbar(data = myDat, aes(x = time, ymin = lci, ymax = uci), width = 0.2, color = "red") +
  labs(title = "True and Fitted HIV Prevalence (SI7 Model)",
       x = "Year",
       y = "Prevalence") +
  theme_minimal() +
  scale_color_manual(values = c("red" = "red", "with_intervention" = "blue", "without_intervention" = "green"), labels = c("True", "With Intervention", "Without Intervention"))



# # Plot MLE fit time series for SI model
fitDat_SI <- simEpidemic(init = init_SI, tseq = tseqMonth, modFunction = SImod_simple, parms = subsParms(optim.vals_SI$par, trueParms))

# Plot MLE fit time series for SI model without intervention
fitDat_SI_no_intervention <- simEpidemic(init = init_SI, tseq = tseqMonth, modFunction = SImod_simple_no_intervention, parms = subsParms(final_parameters_SI_no_intervention, trueParms))


ggplot() +
  geom_line(data = simDat, aes(x = time, y = P), color = "red") +
  geom_line(data = fitDat_SI, aes(x = time, y = P), color = "blue") +
  geom_line(data = fitDat_SI_no_intervention, aes(x = time, y = P), color = "green") +
  geom_point(data = myDat, aes(x = time, y = sampPrev), color = "red", size = 2) +
  #geom_errorbar(data = myDat, aes(x = time, ymin = lci, ymax = uci), width = 0.2, color = "red") +
  labs(title = "True and Fitted HIV Prevalence (SI Model)",
       x = "Year",
       y = "Prevalence") +
  theme_minimal() +
  scale_color_manual(values = c("red" = "red", "with_intervention" = "blue", "without_intervention" = "green"), labels = c("True", "With Intervention", "Without Intervention"))


fitDat <- simEpidemic(init, tseq = tseqMonth, modFunction = SImod, parms = subsParms(optim.vals$par, trueParms))

# Plot MLE fit time series for S-I-I-I-I model without intervention
fitDat_no_intervention <- simEpidemic(init = init, tseq = tseqMonth, modFunction = SImod_no_intervention, parms = subsParms(final_parameters_no_intervention, trueParms))

ggplot() +
  geom_line(data = simDat, aes(x = time, y = P), color = "red") +
  geom_line(data = fitDat, aes(x = time, y = P), color = "blue") +
  geom_line(data = fitDat_no_intervention, aes(x = time, y = P), color = "green") +
  geom_point(data = myDat, aes(x = time, y = sampPrev), color = "red", size = 2) +
  geom_errorbar(data = myDat, aes(x = time, ymin = lci, ymax = uci), width = 0.2, color = "red") +
  labs(title = "True and Fitted HIV Prevalence (SI4 Model)",
       x = "Year",
       y = "Prevalence") +
  theme_minimal() +
  scale_color_manual(values = c("red" = "red", "with_intervention" = "blue", "without_intervention" = "green"), labels = c("True", "With Intervention", "Without Intervention"))



# Calculate log-likelihood for true parameters using the original S-I-I-I-I model
nll_true <- nllikelihood(trueParms, obsDat = myDat, modFunction = SImod, init = init)

# Calculate log-likelihood for estimated parameters for the SI7 model
nll_SI7 <- nllikelihood(subsParms(optim.vals_SI7$par, trueParms), obsDat = myDat, modFunction = SI7mod, init = init_SI7)

# Calculate log-likelihood for estimated parameters for the SI model
nll_SI <- nllikelihood(subsParms(optim.vals_SI$par, trueParms), obsDat = myDat, modFunction = SImod_simple, init = init_SI)

# Calculate log-likelihood for estimated parameters for the SI model
nll_SImod <- nllikelihood(subsParms(optim.vals$par, trueParms), obsDat = myDat, modFunction = SImod, init = init)

# Calculate log-likelihood for estimated parameters for the SI7 model without intervention
nll_SI7_no_intervention <- nllikelihood(subsParms(final_parameters_SI7_no_intervention, trueParms), obsDat = myDat, modFunction = SI7mod_no_intervention, init = init_SI7)

# Calculate log-likelihood for estimated parameters for the SI model without intervention
nll_SI_no_intervention <- nllikelihood(subsParms(final_parameters_SI_no_intervention, trueParms), obsDat = myDat, modFunction = SImod_simple_no_intervention, init = init_SI)

# Calculate log-likelihood for estimated parameters for the S-I-I-I-I model without intervention
nll_no_intervention <- nllikelihood(subsParms(final_parameters_no_intervention, trueParms), obsDat = myDat, modFunction = SImod_no_intervention, init = init)

# Print log-likelihood values
cat("Log-likelihood for true parameters (original model):", nll_true, "\n")
cat("Log-likelihood for estimated parameters (SI7 model with intervention):", nll_SI7, "\n")
cat("Log-likelihood for estimated parameters (SI7 model without intervention):", nll_SI7_no_intervention, "\n")
cat("Log-likelihood for estimated parameters (SI model with intervention):", nll_SI, "\n")
cat("Log-likelihood for estimated parameters (SI model without intervention):", nll_SI_no_intervention, "\n")
cat("Log-likelihood for estimated parameters (SI4 model with intervention):", nll_SImod, "\n")
cat("Log-likelihood for estimated parameters (SI4 model without intervention):", nll_no_intervention, "\n")


##################################################################################################
# Load the ggplot2 library
library(ggplot2)

# Assuming myDat is your dataframe and myColumn is the column you want to plot
# You can adjust the name of the dataframe and column accordingly

# Cumulative incidence ratio ( CIR )   over 100 simulations
avt_results <- data.frame((matrix(ncol=0, nrow = num_simulations)))
avt_results$av_SI <-cum_results$CI_SI/cum_results$CI_SI_no_inter
avt_results$av_SI4 <- cum_results$CI_SI4/cum_results$CI_SI4_no_inter
avt_results$av_SI7 <- cum_results$CI_SI7/cum_results$CI_SI7_no_inter 

####### distribution Cumulative incidence Ratio ("1-proportion of averted cases")

## pair plot
pairs(avt_results,main="correlations checking")

par(mfrow = c(1, 3))
# Plot the density
#plot(density(avt_results$av_SI), main = "distribution of CIR over simulations",ylim=c(0,3))
# Overlay the histogram
hist(avt_results$av_SI, breaks = 50, col = "lightblue", border = "black",alpha=0.5,xlab = "CIR value",main="S-I model stucture")

#plot(density(avt_results$av_SI4), main = "distribution of CIR over simulations",ylim=c(0,3))
hist(avt_results$av_SI4, breaks = 50, col = "lightblue", border = "black",alpha=0.5,xlab = "AR value",main="S-I-I-I-I model stucture")

#plot(density(avt_results$av_SI7), main = "distribution of CIR over simulations",ylim=c(0,3))
hist(avt_results$av_SI7, breaks = 50, col = "lightblue", border = "black",alpha=0.5,xlab = "AR value",main="S-I-I-I-I-I-I-I model stucture")
# Add a legend



########################  parameters evaluations

index <- c("log_lambda","log_a","log_cMax","log_cRate")
# true values
true_params <- c(trueParms$lambda,trueParms$a,trueParms$cMax,trueParms$cRate)

## residuals ( true-estimate)

# SI model 
residuals1 <- true_params-results_df_SI[,index]
names(residuals1) <- c("lambda","a","cMax","cRate") 
#summary statistics of residuals.
summary(residuals1)

# SI4 model 
residuals2 <- true_params-results_df_SImod[,index]
names(residuals2) <- c("lambda","a","cMax","cRate") 
#summary statistics of residuals.
summary(residuals2)

# SI7 model 
residuals3 <- true_params-results_df_SI7[,index]
names(residuals3) <- c("lambda","a","cMax","cRate") 
#summary statistics of residuals.
summary(residuals3)

saveRDS(residuals1,"residuals1.rds")
saveRDS(residuals2,"residuals2.rds")
saveRDS(residuals3,"residuals3.rds")
save(avt_results,"averted_res_2.rds")
# distributions of the residuals
plotter <- function(residuals=residuals2){
  
  #pairs plot (residuals)
  pairs(residuals,main="correlation checking")
  
  # distributions plot
  par(mfrow=c(1,3))
  
  plot(density(residuals$lambda), main = "lambda",ylim=c(0,1),col="red")
  hist(residuals$lambda, breaks = 50, col = "lightblue", border = "black",probability = TRUE, add = TRUE)
  
  plot(density(residuals$a), main = "a",ylim=c(0,1),col="red")
  hist(residuals$a, breaks = 70, col = "lightblue", border = "black",probability = TRUE, add = TRUE)
  
  plot(density(residuals$cMax), main = "cMax",ylim=c(0,1),col="red")
  hist(residuals$cMax, breaks = 50, col = "lightblue", border = "black",probability = TRUE, add = TRUE)
  
  legend("topright", legend = c("estimate density", "Histogram"), fill = c("red", "lightblue"))
  
  #plot(density(residuals1$cRate), main = "cRate",ylim=c(0,1))
  #hist(residuals1$cRate, breaks = 50, col = "lightblue", border = "black",probability = TRUE, add = TRUE)
  
}

# SI
plotter(residuals1)
plotter(residuals2)
plotter(residuals3)







