# Libraries
library(GeoModels)
library(fields)

# Load necessary data
coords <- as.matrix(read.table("coords1.txt", header = TRUE))
data <- read.table("data1.txt", header = TRUE)$data
loc_to_pred <- as.matrix(read.table("loc_to_pred.txt", header = TRUE))

# Plot semivariogram (assuming isotropy)
semivariogram <- GeoVariogram(coordx = coords, data = data, maxdist = 0.5)
plot(semivariogram, xlab = "h", ylab = expression(gamma(h)), main = "Isotropic Semivariogram")

# Plot data distribution
hist(data, main = "Data Distribution", xlab = "Values", ylab = "Frequency")

# Fit the first model (Exponential correlation)
model_1 <- GeoFit(
  data = data,
  coordx = coords,
  corrmodel = "Exponential",
  likelihood = "Marginal",
  type = "Pairwise",
  start = list(mean = 0, sill = 1, scale = 0.1),
  lower = list(mean = -Inf, sill = 0, scale = 0),
  upper = list(mean = Inf, sill = Inf, scale = Inf)
)

# Fit the second model (Matern correlation with ν = 3/2)
model_2 <- GeoFit(
  data = data,
  coordx = coords,
  corrmodel = "Matern",
  likelihood = "Marginal",
  type = "Pairwise",
  start = list(mean = 0, sill = 1, scale = 0.1),
  lower = list(mean = -Inf, sill = 0, scale = 0),
  upper = list(mean = Inf, sill = Inf, scale = Inf),
  fixed = list(smooth = 1.5)
)

# Print parameter estimates for both models
cat("Model 1 (Exponential) Parameters:\n")
print(model_1$param)
cat("Model 2 (Matern) Parameters:\n")
print(model_2$param)

# Compare models using AIC values
cat("AIC Values:\n")
cat("Model 1 (Exponential):", model_1$logCompLik, "\n")
cat("Model 2 (Matern):", model_2$logCompLik, "\n")

# Cross-validation for RMSE
set.seed(123)  
n_iterations <- 100
test_fraction <- 0.1
rmse_model_1 <- numeric(n_iterations)
rmse_model_2 <- numeric(n_iterations)

for (i in 1:n_iterations) {
  print(paste("Iteration", i, "of", n_iterations))
  n <- length(data)
  test_indices <- sample(1:n, size = floor(test_fraction * n))
  train_indices <- setdiff(1:n, test_indices)
  
  train_data <- data[train_indices]
  train_coords <- coords[train_indices, ]
  test_data <- data[test_indices]
  test_coords <- coords[test_indices, ]
  
  # Fit models on training data
  model_1 <- GeoFit(
    data = train_data,
    coordx = train_coords,
    corrmodel = "Exponential",
    likelihood = "Marginal",
    type = "Pairwise",
    start = list(mean = 0, sill = 1, scale = 0.1),
    lower = list(mean = -Inf, sill = 0, scale = 0),
    upper = list(mean = Inf, sill = Inf, scale = Inf)
  )
  
  model_2 <- GeoFit(
    data = train_data,
    coordx = train_coords,
    corrmodel = "Matern",
    likelihood = "Marginal",
    type = "Pairwise",
    start = list(mean = 0, sill = 1, scale = 0.1),
    lower = list(mean = -Inf, sill = 0, scale = 0),
    upper = list(mean = Inf, sill = Inf, scale = Inf),
    fixed = list(smooth = 1.5)
  )
  
  # Predict test data
  predictions_1 <- GeoKrig(estobj = model_1, loc = test_coords, mse = FALSE)$pred
  predictions_2 <- GeoKrig(estobj = model_2, loc = test_coords, mse = FALSE)$pred
  
  # Compute RMSE
  rmse_model_1[i] <- sqrt(mean((test_data - predictions_1)^2))
  rmse_model_2[i] <- sqrt(mean((test_data - predictions_2)^2))
}

# Calculate average RMSE
mean_rmse_model_1 <- mean(rmse_model_1)
mean_rmse_model_2 <- mean(rmse_model_2)

cat("Average RMSE Values:\n")
cat("Model 1 (Exponential):", mean_rmse_model_1, "\n")
cat("Model 2 (Matern):", mean_rmse_model_2, "\n")

# Choose best model based on RMSE
if (mean_rmse_model_1 < mean_rmse_model_2) {
  chosen_model <- model_1
  best_model <- "Exponential"
} else {
  chosen_model <- model_2
  best_model <- "Matern (ν = 3/2)"
}
cat("Best Model Based on RMSE:", best_model, "\n")

# Perform kriging predictions on loc_to_pred
kriging_result <- GeoKrig(estobj = chosen_model, loc = loc_to_pred, mse = TRUE)

# Extract predictions and MSE
predictions <- kriging_result$pred
mse <- kriging_result$mse

# Reshape predictions and MSE for plotting
x_coords <- unique(loc_to_pred[, 1])
y_coords <- unique(loc_to_pred[, 2])
grid_size_x <- length(x_coords)
grid_size_y <- length(y_coords)

prediction_matrix <- matrix(predictions, nrow = grid_size_x, ncol = grid_size_y)
mse_matrix <- matrix(mse, nrow = grid_size_x, ncol = grid_size_y)

# Plot predictions and MSE
par(mfrow = c(1, 2))
image.plot(x_coords, y_coords, prediction_matrix,
           main = "Kriging Prediction",
           xlab = "X", ylab = "Y")
image.plot(x_coords, y_coords, mse_matrix,
           main = "Kriging MSE",
           xlab = "X", ylab = "Y")



