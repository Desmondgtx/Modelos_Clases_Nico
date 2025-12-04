
## Generalized Context Model (GCM)
# Diego Garrido - Nicolás Marchant
# Viña, 2025

########################
# Seed para reproducibilidad
set.seed(42)

########################
# DEFINICIÓN DEL MODELO GCM

# Predicción del GCM

gcm_predict = function(stim, exemplars, c, w, gamma, cat_A = c(3, 4, 7, 8)) {
  
  # stim:      Matriz de estímulos (n x d)
  # exemplars: Matriz de ejemplares almacenados (m x d)
  # c:         Parámetro de sensibilidad/especificidad
  # w:         Vector de pesos atencionales
  # gamma:     Parámetro de sesgo de respuesta
  # cat_A:     Índices de ejemplares de categoría A
  
  
  
  n_stim = nrow(stim)
  n_ex = nrow(exemplars)
  w = w / sum(w)  # Normalizar pesos
  
  # Calcular matriz de similitudes
  similarities = matrix(0, n_stim, n_ex)
  for (i in 1:n_stim) {
    for (j in 1:n_ex) {
      dist_ij = sum(w * abs(stim[i, ] - exemplars[j, ]))  # City-block
      similarities[i, j] = exp(-c * dist_ij)  # Similitud exponencial
    }
  }
  
  # Sumar similitudes por categoría
  cat_B = setdiff(1:n_ex, cat_A)
  sum_A = rowSums(similarities[, cat_A, drop = FALSE])
  sum_B = rowSums(similarities[, cat_B, drop = FALSE])
  
  # Regla de Luce con gamma
  prob_A = (sum_A^gamma) / (sum_A^gamma + sum_B^gamma)
  
  return(prob_A) # Vector de probabilidades P(A) para cada estímulo
}

########################
# FUNCIONES DE AJUSTE

# Negative Log-Likelihood (binomial)
gcm_nll = function(theta, stim, exemplars, data, N, cat_A = c(3, 4, 7, 8)) {
  c = theta[1]
  w = theta[2:4]
  gamma = theta[5]
  
  prob_pred = gcm_predict(stim, exemplars, c, w, gamma, cat_A)
  prob_pred = pmax(pmin(prob_pred, 1 - 1e-10), 1e-10)  # Evitar log(0)
  
  k = data * N
  nll = -sum(k * log(prob_pred) + (N - k) * log(1 - prob_pred))
  return(nll)
}

# RMSE
gcm_rmse = function(theta, stim, exemplars, data, cat_A = c(3, 4, 7, 8)) {
  c = theta[1]
  w = theta[2:4]
  gamma = theta[5]
  
  prob_pred = gcm_predict(stim, exemplars, c, w, gamma, cat_A)
  return(sqrt(mean((data - prob_pred)^2)))
}

########################
# FUNCIÓN DE AJUSTE PRINCIPAL

# Ajustar GCM a datos
gcm_fit = function(data, stim, exemplars, N = 10, cat_A = c(3, 4, 7, 8),
                   method = "nll") {
  # data:      Vector de proporciones observadas
  # stim:      Matriz de estímulos
  # exemplars: Matriz de ejemplares
  # N:         Número de ensayos por estímulo
  # method:   "nll" (default) o "rmse"
  
  theta_init = c(1, 0.33, 0.33, 0.33, 1)
  lower = c(0.1, 0.01, 0.01, 0.01, 0.1)
  upper = c(10, 1, 1, 1, 5)
  
  if (method == "nll") {
    fit = optim(par = theta_init, fn = gcm_nll, method = "L-BFGS-B",
                lower = lower, upper = upper,
                stim = stim, exemplars = exemplars, data = data, 
                N = N, cat_A = cat_A)
  } else {
    fit = optim(par = theta_init, fn = gcm_rmse, method = "L-BFGS-B",
                lower = lower, upper = upper,
                stim = stim, exemplars = exemplars, data = data, 
                cat_A = cat_A)
  }
  
  # Extraer parámetros
  c_fit = fit$par[1]
  w_fit = fit$par[2:4]
  w_norm = w_fit / sum(w_fit)
  gamma_fit = fit$par[5]
  
  # Predicciones
  pred = gcm_predict(stim, exemplars, c_fit, w_norm, gamma_fit, cat_A)
  
  # Métricas
  rmse = sqrt(mean((data - pred)^2))
  ss_res = sum((data - pred)^2)
  ss_tot = sum((data - mean(data))^2)
  r2 = ifelse(ss_tot > 0, 1 - ss_res/ss_tot, NA)
  
  n_params = 5
  n_obs = length(data) * N
  nll = gcm_nll(fit$par, stim, exemplars, data, N, cat_A)
  AIC = 2 * nll + 2 * n_params
  BIC = 2 * nll + n_params * log(n_obs)
  
  return(list(
    c = c_fit,
    w = w_norm,
    gamma = gamma_fit,
    predictions = pred,
    nll = nll,
    rmse = rmse,
    r2 = r2,
    AIC = AIC,
    BIC = BIC,
    convergence = fit$convergence
  )) # Lista con parámetros ajustados y métricas
}

########################
# DATOS SINTÉTICOS Y DEMOSTRACIÓN

# Definir ejemplares (cubo 3D)
exemplars = matrix(c(
  -1, -1, -1,  # 1 - Cat B
  -1, -1,  1,  # 2 - Cat B
  -1,  1, -1,  # 3 - Cat A
  -1,  1,  1,  # 4 - Cat A
  1, -1, -1,  # 5 - Cat B
  1, -1,  1,  # 6 - Cat B
  1,  1, -1,  # 7 - Cat A
  1,  1,  1   # 8 - Cat A
), nrow = 8, byrow = TRUE)

cat_A = c(3, 4, 7, 8)  # Categoría A: dim 2 = 1
stim = exemplars

# Parámetros verdaderos
true_c = 2.5
true_w = c(0.2, 0.6, 0.2)
true_gamma = 1.0

# Generar datos
true_probs = gcm_predict(stim, exemplars, true_c, true_w, true_gamma, cat_A)
N_trials = 10
observed_counts = rbinom(8, N_trials, true_probs)
observed_probs = observed_counts / N_trials


cat("Probabilidades verdaderas:\n")
print(round(true_probs, 3))
cat("\nProporciones observadas:\n")
print(round(observed_probs, 3))


########################
# AJUSTE DEL MODELO

result = gcm_fit(observed_probs, stim, exemplars, N_trials, cat_A, method = "nll")

cat("\n========== RESULTADOS DEL AJUSTE ==========\n")
cat("\nParámetros recuperados:\n")
cat("  c (sensibilidad):", round(result$c, 3), "\n")
cat("  w (pesos):", round(result$w, 3), "\n")
cat("  gamma (sesgo):", round(result$gamma, 3), "\n")

cat("\nMétricas de ajuste:\n")
cat("  RMSE:", round(result$rmse, 4), "\n")
cat("  R²:", round(result$r2, 4), "\n")
cat("  AIC:", round(result$AIC, 2), "\n")
cat("  BIC:", round(result$BIC, 2), "\n")

cat("\n========== COMPARACIÓN ==========\n")
cat("\nParámetro | Verdadero | Recuperado\n")
cat("----------|-----------|------------\n")
cat(sprintf("c         | %9.3f | %10.3f\n", true_c, result$c))
cat(sprintf("w1        | %9.3f | %10.3f\n", true_w[1], result$w[1]))
cat(sprintf("w2        | %9.3f | %10.3f\n", true_w[2], result$w[2]))
cat(sprintf("w3        | %9.3f | %10.3f\n", true_w[3], result$w[3]))
cat(sprintf("gamma     | %9.3f | %10.3f\n", true_gamma, result$gamma))


########################
# AJUSTE PARA MÚLTIPLES SUJETOS (EJEMPLO)

n_subjects = 5
subjects_data = matrix(0, n_subjects, 8)
subjects_results = vector("list", n_subjects)

# Simular datos para cada sujeto
for (s in 1:n_subjects) {
  subjects_data[s, ] = rbinom(8, N_trials, true_probs) / N_trials
  subjects_results[[s]] = gcm_fit(subjects_data[s, ], stim, exemplars, 
                                  N_trials, cat_A)
}

# Resumen de parámetros
params_summary = data.frame(
  Subject = 1:n_subjects,
  c = sapply(subjects_results, function(x) x$c),
  w1 = sapply(subjects_results, function(x) x$w[1]),
  w2 = sapply(subjects_results, function(x) x$w[2]),
  w3 = sapply(subjects_results, function(x) x$w[3]),
  gamma = sapply(subjects_results, function(x) x$gamma),
  RMSE = sapply(subjects_results, function(x) x$rmse),
  R2 = sapply(subjects_results, function(x) x$r2)
)

print(round(params_summary, 3))

cat("\nPromedios:\n")
cat("  c:", round(mean(params_summary$c), 3), "\n")
cat("  w:", round(colMeans(params_summary[, c("w1", "w2", "w3")]), 3), "\n")
cat("  gamma:", round(mean(params_summary$gamma), 3), "\n")
cat("  RMSE:", round(mean(params_summary$RMSE), 4), "\n")

########################
# VISUALIZACIÓN

# Guardar gráfico como PNG
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# Gráfico de Barras comparativas
barplot_matrix = rbind(observed_probs, result$predictions)
bp = barplot(barplot_matrix, 
             beside = TRUE, 
             names.arg = 1:8,
             col = c("steelblue", "coral"),
             ylim = c(0, 1.1),
             xlab = "Estímulo",
             ylab = "P(Categoría A)",
             main = "GCM: Observado vs Predicho")
legend("topright", 
       legend = c("Observado", "Predicho"), 
       fill = c("steelblue", "coral"),
       cex = 0.9)
x_pos = colMeans(bp)
points(x_pos, true_probs, pch = 4, cex = 1.5, lwd = 2)
legend("topleft", legend = "Verdadero", pch = 4, pt.lwd = 2, cex = 0.9)

# Gráfico de dispersión
plot(observed_probs, result$predictions,
     pch = 19, cex = 2, col = "steelblue",
     xlim = c(0, 1), ylim = c(0, 1),
     xlab = "P(A) Observado",
     ylab = "P(A) Predicho",
     main = "Ajuste del Modelo",
     asp = 1)
abline(0, 1, lty = 2, col = "gray40", lwd = 2)
text(observed_probs, result$predictions, labels = 1:8, pos = 3, cex = 0.8)

mtext(paste0("RMSE = ", round(result$rmse, 4), 
             " | R² = ", round(result$r2, 4)),
      side = 1, line = 4, cex = 0.9)
