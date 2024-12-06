# Cargar librerías necesarias
library(GeoModels)
library(geoR)

# Cargar coordenadas espaciales desde el archivo coords1.txt
coords <- read.table("coords1.txt", header = TRUE)

# Definir tiempos
times <- c(1, 2, 3, 4)

# Número de ubicaciones espaciales y tiempos
n_spatial <- nrow(coords)   # 500 ubicaciones espaciales
n_temporal <- length(times) # 4 tiempos únicos

# Parámetro inicial del modelo
alpha_s <- 2/3  # Escala espacial según lo que nos indica el problema
alpha_t <- 0.2 / 3  # Escala temporal
nugget <- 0.1    # Nugget inicial
sill <- 1        # Sill inicial

# Crear matriz de correlación espacio-temporal usando el producto Kronecker
dist_spatial <- as.matrix(dist(coords)) # Distancias entre ubicaciones espaciales
dist_temporal <- as.matrix(dist(times)) # Distancias entre tiempos

# Correlación espacio-temporal, usando el parámetro alpha_s para la escala espacial y alpha_t para la temporal
rho_matrix <- kronecker(
  exp(-dist_temporal / alpha_t),  # Correlación temporal
  exp(-dist_spatial / alpha_s)   # Correlación espacial (ahora usa alpha_s)
)

# Verificar que rho_matrix sea positiva definida
if (any(eigen(rho_matrix)$values <= 0)) {
  stop("La matriz rho_matrix no es positiva definida.")
}

# Simular datos espacio-temporales usando descomposición de Cholesky
set.seed(123) # Asegurar reproducibilidad
Z <- t(chol(rho_matrix)) %*% rnorm(n_spatial * n_temporal, mean = 0, sd = 1)

# Convertir los datos simulados en una matriz t x d (temporal x espacial)
Y <- matrix(Z, nrow = n_temporal, ncol = n_spatial)

# Validar dimensiones de los datos simulados
print(paste("Dimensiones de Y (datos simulados):", dim(Y)[1], "x", dim(Y)[2])) # Debe ser 4 x 500

# Parámetros iniciales para GeoFit
start_params <- list(scale = alpha_s, nugget = nugget, sill = sill)

# Ajustar el modelo espacio-temporal usando GeoFit
fit <- GeoFit(
  coordx = coords,       # Coordenadas espaciales
  coordt = times,        # Coordenadas temporales
  data = Y,              # Matriz t x d de observaciones espacio-temporales
  corrmodel = "exponential",  # Modelo de correlación
  likelihood = "Marginal",    # Verosimilitud marginal
  type = "Pairwise",          # Tipo de aproximación
  start = start_params        # Parámetros iniciales
)

# Mostrar resultados del ajuste
print(fit)

# Definir coordenadas para la predicción
coords_pred <- matrix(c(0.3, 0.2,  # Primera ubicación
                        0.5, 0.5), # Segunda ubicación
                      ncol = 2, byrow = TRUE)

# Definir tiempos para predicción
time_pred <- c(2, 2) # Tiempo asociado a cada coordenada en coords_pred

# Validar dimensiones y datos
if (dim(coords)[1] != ncol(Y)) stop("El número de filas en 'coords' debe coincidir con las columnas de 'Y'")
if (length(times) != nrow(Y)) stop("La longitud de 'times' debe coincidir con las filas de 'Y'")
if (nrow(coords_pred) != length(time_pred)) stop("El número de filas en 'coords_pred' debe coincidir con la longitud de 'time_pred'")

# GeoKrig con el objeto 'fit' de GeoFit
krig_result <- GeoKrig(
  estobj = fit,                # Usar el objeto de ajuste directamente
  loc = coords_pred,           # Coordenadas para predicción
  time = time_pred,            # Tiempos para predicción
  mse = TRUE                   # Calcular errores cuadrados medios
)

# Mostrar resultados de la predicción
print("Predicción en las ubicaciones:")
print(krig_result$pred)
print("Varianza de predicción:")
print(krig_result$mse)

# Graficar las ubicaciones y predicciones
plot(coords, col = "blue", pch = 16, main = "Ubicaciones de predicción",
     xlab = "Coordenada X", ylab = "Coordenada Y")
points(coords_pred, col = "red", pch = 17)
legend("topright", legend = c("Datos observados", "Puntos de predicción"),
       col = c("blue", "red"), pch = c(16, 17))
