
## Modeling Learning Rate 
## Computational Modeling of Cognition and Behavior

# Libraries
library(ggplot2)

# Seed
set.seed(123)

# Power Law Function model
power_law = function(N, beta, scale = 6000) 
{ # Where N = number of trial, beta = learning rate
    scale * N^(-beta)
}

# Exponential Function model
exp_model = function(N, alpha, scale = 6000)
{ # Where alpha = learning rate
  scale * exp(-alpha*N)
}

# Trials
trials = 1:150

# Parameters
beta = 0.4 # Power Law Function, learning rate
alpha = 0.5 # Exponential Function, learning rate

# Model predictions
power_predictions = power_law(trials, beta)
exp_predictions = exp_model(trials, alpha)

# Generate data
noise_factor = 200
observed_data = power_predictions + rnorm(length(trials), mean = 0, sd = noise_factor)


# Ensure there are no negative RTs
observed_data[observed_data < 100] = 100 + abs(rnorm(sum(observed_data < 100), 0, 50))


# Create dataframe
df = data.frame(
  Trial = trials,
  Observed = observed_data,
  Power_Law = power_predictions,
  Exponential = exp_predictions
)

# Create plot
ggplot(df, aes(x = Trial)) +
  geom_point(aes(y = Observed), size = 1.5) +
  geom_line(aes(y = Power_Law), linetype = "dashed", size = 1) +
  geom_line(aes(y = Exponential), linetype = "solid", size = 1) +
  labs(x = "Trial Number",
       y = "Response Time (ms)",
       title = "Learning Cruves: Power Law vs Exponential") + 
  theme_minimal() +
  ylim(0, 7000)

# Calculate Residual Mean-Squared Deviation (RMSD)
rmsd_power = sqrt(mean((df$Observed - df$Power_Law)^2))
rmsd_exp = sqrt(mean((df$Observed - df$Exponential)^2))


cat("Power Law RMSD:", rmsd_power, "\n")
cat("Exponential RMSD:", rmsd_exp, "\n")





