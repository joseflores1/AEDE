# Load required libraries
library(GeoModels)
library(fields)


set.seed(123)

# Load data
data <- read.table("data_temp.txt", header = TRUE)

# Extract spatial coordinates and variables
coords <- as.matrix(data[, c("x", "y")])
altitudes <- data$alt
temperature <- data$temp

# Spatial regression model
X <- cbind(1, altitudes)  # Regression matrix (intercept and altitude)
corr_model <- "Wend0"     # Wendland correlation model

# Configure initial parameters
start_params <- list(mean = mean(temperature), 
                     mean1 = 0,         # Initial effect of altitude
                     scale = 500,       # Initial spatial range
                     sill = var(temperature))  # Initial variance

# Fix specific parameters for the Wendland model
fixed_params <- list(power2 = 4)

# Set bounds for parameter values
lower_bounds <- list(mean = -Inf, mean1 = -Inf, scale = 0.01, sill = 0.01)
upper_bounds <- list(mean = Inf, mean1 = Inf, scale = Inf, sill = Inf)

# Fit the model using GeoFit
fit <- GeoFit(
  data = temperature,        # Temperature
  coordx = coords,           # Spatial coordinates
  X = X,                     # Regression matrix
  corrmodel = corr_model,    # Correlation model
  likelihood = "Conditional",# Conditional likelihood
  type = "Pairwise",         # Pairwise likelihood
  start = start_params,      # Initial parameters
  fixed = fixed_params,      # Fixed parameters
  lower = lower_bounds,      # Lower bounds
  upper = upper_bounds,      # Upper bounds
  neighb = 3,                # Neighborhood order 3
  optimizer = "BFGS",        # Optimization method
  sensitivity = TRUE         # Sensitivity matrix for bootstrap
)

# Display parameter estimates
print("Estimated parameters:")
print(fit$param)

# Perform bootstrap for standard error estimation
bootstrap_fit <- GeoVarestbootstrap(fit, K = 100)  # 100 bootstrap iterations
print("Bootstrap standard errors:")
print(bootstrap_fit$stderr)

emp_variogram <- GeoVariogram(
  data = temperature,
  coordx = coords,
  maxdist = 400,
  numbins = 200    # Mayor resolución en los intervalos
)

# Plot the empirical variogram
plot(emp_variogram, pch = 20, main = "Empirical Semivariogram of Residuals")

# Compare with the fitted (theoretical) variogram
GeoCovariogram(
  fitted = fit,              # Fitted model
  show.vario = TRUE,         # Show variogram
  vario = emp_variogram,     # Compare with empirical variogram
  pch = 20
)


