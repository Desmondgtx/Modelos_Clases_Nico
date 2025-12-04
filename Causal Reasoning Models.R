
# Modelos de razonamiento causal (normativo y Mutation Sampler)
# Diego Garrido - Nicolás Marchant
# Viña, 2025

########################
# Seed para reproducibilidad
set.seed(123)
options(scipen = 999)


########################
# SECCIÓN 1: FUNCIONES HELPER

# Crear el DAG base con todos los estados posibles
create_base_dag <- function() {
  dag <- data.frame(
    X1 = c(1,1,1,1,0,0,0,0),
    Y = c(1,1,0,0,1,1,0,0),
    X2 = c(1,0,1,0,1,0,1,0)
  )
  rownames(dag) <- c("111", "110", "101", "100", "011", "010", "001", "000")
  return(dag) # dataframe con 8 estados (combinaciones de X1, Y, X2)
}

# Obtener estructura de padres según el tipo de DAG
get_parent_structure <- function(dag_type, dag_df) {
  # dag_type: Tipo de estructura ("commoncause", "chain", o "commoneffect")
  # dag_df: Dataframe del DAG base

  n <- nrow(dag_df)
  
  if (dag_type == "commoncause") {
    # X1 <- Y -> X2  (Y exógeno)
    dag_df$parentsX1 <- lapply(dag_df$Y, function(y) c(y))
    dag_df$parentsY  <- replicate(n, numeric(0), simplify = FALSE)
    dag_df$parentsX2 <- lapply(dag_df$Y, function(y) c(y))
    
  } else if (dag_type == "chain") {
    # X1 -> Y -> X2  (X1 exógeno)
    dag_df$parentsX1 <- replicate(n, numeric(0), simplify = FALSE)
    dag_df$parentsY  <- lapply(dag_df$X1, function(x) c(x))
    dag_df$parentsX2 <- lapply(dag_df$Y, function(y) c(y))
    
  } else if (dag_type == "commoneffect") {
    # X1 -> Y <- X2  (X1, X2 exógenos)
    dag_df$parentsX1 <- replicate(n, numeric(0), simplify = FALSE)
    dag_df$parentsY  <- mapply(function(x1, x2) c(x1, x2), 
                               dag_df$X1, dag_df$X2, SIMPLIFY = FALSE)
    dag_df$parentsX2 <- replicate(n, numeric(0), simplify = FALSE)
    
  } else {
    stop("Tipo de DAG inválido. Opciones: 'commoncause', 'chain', o 'commoneffect'")
  }
  
  dag_df # DAG con columnas de padres agregadas
}

# Función de probabilidad Noisy-OR
prob <- function(var, parents = c(), c, m, b) {
  # var:     Valor de la variable (0 o 1)
  # parents: Vector con valores de los padres
  # c:       Probabilidad base (causa exógena)
  # m:       Fuerza causal (power of the cause)
  # b:       Probabilidad de fondo (background rate)
  
  exogenous <- length(parents) == 0
  if (exogenous) {
    p <- c
  } else {
    cn <- sum(parents)
    p <- 1 - (1 - b) * (1 - m)^cn
  }
  return(var * p + (1 - var) * (1 - p)) # Probabilidad
}

# Construir matriz de adyacencia para el Mutation Sampler
build_adjmat <- function() {
  states <- expand.grid(X1 = c(1,0), Y = c(1,0), X2 = c(1,0))
  states <- states[order(-states$X1, -states$Y, -states$X2), ]
  
  n <- nrow(states)
  adj <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      differences <- sum(states[i, ] != states[j, ])
      if (differences == 1) {
        adj[i, j] <- 1
      }
    }
  }
  return(adj) # Matriz 8x8 donde 1 indica estados vecinos
}

# Función auxiliar para calcular probabilidades condicionales
cond_prob <- function(numerator_states, denominator_states, probs) {
  # numerator_states:   Índices de estados en el numerador
  # denominator_states: Índices de estados en el denominador
  # probs:              Vector de probabilidades
  
  num <- sum(probs[numerator_states])
  denom <- sum(probs[denominator_states])
  if (denom == 0 || is.nan(denom) || is.na(denom)) {
    return(0.5)
  }
  result <- num / denom
  if (is.na(result) || is.nan(result) || !is.finite(result)) {
    return(0.5)
  }
  return(result) # Probabilidad condicional
}

########################
# FUNCIONES DE PREDICCIÓN

# Predicciones del Modelo Normativo
PredNormative <- function(c, m, b, dag_structure) {
  # c,m,b:         Parámetros del modelo Noisy-OR
  # dag_structure: Tipo de estructura causal
  
  dag <- create_base_dag()
  dag <- get_parent_structure(dag_structure, dag)
  
  # Calcular probabilidades individuales
  dag$px1 <- mapply(function(x, p) prob(x, p, c = c, m = m, b = b), 
                    dag$X1, dag$parentsX1)
  dag$py  <- mapply(function(x, p) prob(x, p, c = c, m = m, b = b), 
                    dag$Y, dag$parentsY)
  dag$px2 <- mapply(function(x, p) prob(x, p, c = c, m = m, b = b), 
                    dag$X2, dag$parentsX2)
  dag$jointprob <- dag$px1 * dag$py * dag$px2
  
  inf_data <- data.frame(
    cases = c("X2=1|Y=0,X1=0", "X2=1|Y=1,X1=0", "X2=1|X1=0", 
              "X2=1|Y=0,X1=1", "X2=1|Y=1,X1=1", "X2=1|X1=1", 
              "X2=1|Y=0", "X2=1|Y=1"),
    resp = c(cond_prob(7, c(7, 8), dag$jointprob),
             cond_prob(5, c(5, 6), dag$jointprob),
             cond_prob(c(5, 7), c(5, 7, 6, 8), dag$jointprob),
             cond_prob(3, c(3, 4), dag$jointprob),
             cond_prob(1, c(1, 2), dag$jointprob),
             cond_prob(c(1, 3), c(1, 3, 2, 4), dag$jointprob),
             cond_prob(c(3, 7), c(3, 7, 4, 8), dag$jointprob),
             cond_prob(c(1, 5), c(1, 5, 2, 6), dag$jointprob))
  )
  
  return(inf_data) # Dataframe con inferencias y valores predichos
}

# Predicciones del Mutation Sampler

PredMutationSampler <- function(len, c, m, b, dag_structure, 
                                use_poisson = FALSE, bias = 0.5) {
  # len:           Largo de la cadena MCMC
  # c,m,b:         Parámetros del modelo Noisy-OR
  # dag_structure: Tipo de estructura causal
  # use_poisson:   Si TRUE, len es el lambda de una Poisson
  # bias:          Sesgo hacia estado inicial (1 = todo 1s, 8 = todo 0s)
  
  
  if (use_poisson) {
    actual_len <- rpois(1, len - 2) + 2
  } else {
    actual_len <- round(len)
  }
  
  startstate <- sample(c(1, 8), 1, prob = c(bias, 1 - bias))
  
  dag <- create_base_dag()
  dag <- get_parent_structure(dag_structure, dag)
  
  # Calcular probabilidades individuales
  dag$px1 <- mapply(function(x, p) prob(x, p, c = c, m = m, b = b), 
                    dag$X1, dag$parentsX1)
  dag$py  <- mapply(function(x, p) prob(x, p, c = c, m = m, b = b), 
                    dag$Y, dag$parentsY)
  dag$px2 <- mapply(function(x, p) prob(x, p, c = c, m = m, b = b), 
                    dag$X2, dag$parentsX2)
  dag$jointprob <- dag$px1 * dag$py * dag$px2
  
  adjmat <- build_adjmat()
  
  # Algoritmo Metropolis-Hastings
  x <- rep(0, actual_len)
  x[1] <- startstate
  
  for (i in 2:actual_len) {
    currentstate <- x[i-1]
    proposedstate <- which(adjmat[currentstate, ] %in% 1)[sample(1:3, 1)]
    A <- dag$jointprob[proposedstate] / dag$jointprob[currentstate]
    
    if (is.na(A) || is.nan(A) || !is.finite(A)) {
      A <- 0.5
    }
    
    if (runif(1) < A) {
      x[i] <- proposedstate
    } else {
      x[i] <- currentstate
    }
  }
  
  samples <- x
  state_counts <- table(factor(samples, levels = 1:8))
  state_probs <- state_counts / length(samples)
  dag$mcmc_probs <- as.numeric(state_probs)
  
  inf_data <- data.frame(
    cases = c("X2=1|Y=0,X1=0", "X2=1|Y=1,X1=0", "X2=1|X1=0", 
              "X2=1|Y=0,X1=1", "X2=1|Y=1,X1=1", "X2=1|X1=1", 
              "X2=1|Y=0", "X2=1|Y=1"),
    resp = c(cond_prob(7, c(7, 8), dag$mcmc_probs),
             cond_prob(5, c(5, 6), dag$mcmc_probs),
             cond_prob(c(5, 7), c(5, 7, 6, 8), dag$mcmc_probs),
             cond_prob(3, c(3, 4), dag$mcmc_probs),
             cond_prob(1, c(1, 2), dag$mcmc_probs),
             cond_prob(c(1, 3), c(1, 3, 2, 4), dag$mcmc_probs),
             cond_prob(c(3, 7), c(3, 7, 4, 8), dag$mcmc_probs),
             cond_prob(c(1, 5), c(1, 5, 2, 6), dag$mcmc_probs))
  )
  
  return(inf_data) # Dataframe con inferencias y valores predichos
}

########################
# FUNCIONES DE FITTING

# Función objetivo para fittear el Modelo Normativo
FitNormative <- function(par, data, dag_structure) {
  # par:           Vector de parámetros [c, m, b]
  # data:          Dataframe con columnas 'cases' y 'resp'
  # dag_structure: Tipo de estructura causal
  
  c_par = par[1]
  m_par = par[2]
  b_par = par[3]
  
  dag <- create_base_dag()
  dag <- get_parent_structure(dag_structure, dag)
  
  # Calcular probabilidades individuales
  dag$px1 <- mapply(function(x, p) prob(x, p, c = c_par, m = m_par, b = b_par), 
                    dag$X1, dag$parentsX1)
  dag$py  <- mapply(function(x, p) prob(x, p, c = c_par, m = m_par, b = b_par), 
                    dag$Y, dag$parentsY)
  dag$px2 <- mapply(function(x, p) prob(x, p, c = c_par, m = m_par, b = b_par), 
                    dag$X2, dag$parentsX2)
  dag$jointprob <- dag$px1 * dag$py * dag$px2
  
  predictions <- c(
    cond_prob(7, c(7, 8), dag$jointprob),
    cond_prob(5, c(5, 6), dag$jointprob),
    cond_prob(c(5, 7), c(5, 7, 6, 8), dag$jointprob),
    cond_prob(3, c(3, 4), dag$jointprob),
    cond_prob(1, c(1, 2), dag$jointprob),
    cond_prob(c(1, 3), c(1, 3, 2, 4), dag$jointprob),
    cond_prob(c(3, 7), c(3, 7, 4, 8), dag$jointprob),
    cond_prob(c(1, 5), c(1, 5, 2, 6), dag$jointprob)
  )
  
  mse <- mean((data$resp - predictions)^2)
  return(mse) # MSE entre predicciones y datos
} 

#Función objetivo para fittear el Mutation Sampler
FitMutationSampler <- function(par, data, dag_structure, 
                               use_poisson = FALSE, bias = 0.5, n_chains = 10) {
  # par:           Vector de parámetros [c, m, b, len]
  # data:          Dataframe con columnas 'cases' y 'resp'
  # dag_structure: Tipo de estructura causal
  # use_poisson:   Si TRUE, len es el lambda de una Poisson
  # bias:          Sesgo hacia estado inicial
  # n_chains:      Número de cadenas MCMC para promediar (reduce varianza)
  
  tryCatch({
    c_par = par[1]
    m_par = par[2]
    b_par = par[3]
    len = par[4]
    
    if (any(is.na(par)) || any(is.nan(par)) || any(!is.finite(par))) {
      return(1e10)
    }
    if (c_par <= 0 || c_par >= 1 || m_par <= 0 || m_par >= 1 || 
        b_par <= 0 || b_par >= 1 || len < 2) {
      return(1e10)
    }
    
    # Promediar sobre múltiples cadenas para reducir varianza
    all_probs <- matrix(0, nrow = n_chains, ncol = 8)
    
    for (chain in 1:n_chains) {
      if (use_poisson) {
        actual_len <- rpois(1, len - 2) + 2
      } else {
        actual_len <- round(len)
      }
      
      startstate <- sample(c(1, 8), 1, prob = c(bias, 1 - bias))
      
      dag <- create_base_dag()
      dag <- get_parent_structure(dag_structure, dag)
      
      # Calcular probabilidades individuales
      dag$px1 <- mapply(function(x, p) prob(x, p, c = c_par, m = m_par, b = b_par), 
                        dag$X1, dag$parentsX1)
      dag$py  <- mapply(function(x, p) prob(x, p, c = c_par, m = m_par, b = b_par), 
                        dag$Y, dag$parentsY)
      dag$px2 <- mapply(function(x, p) prob(x, p, c = c_par, m = m_par, b = b_par), 
                        dag$X2, dag$parentsX2)
      dag$jointprob <- dag$px1 * dag$py * dag$px2
      
      adjmat <- build_adjmat()
      
      x <- rep(0, actual_len)
      x[1] <- startstate
      
      for (i in 2:actual_len) {
        currentstate <- x[i-1]
        proposedstate <- which(adjmat[currentstate, ] %in% 1)[sample(1:3, 1)]
        A <- dag$jointprob[proposedstate] / dag$jointprob[currentstate]
        
        if (is.na(A) || is.nan(A) || !is.finite(A)) {
          A <- 0.5
        }
        
        if (runif(1) < A) {
          x[i] <- proposedstate
        } else {
          x[i] <- currentstate
        }
      }
      
      state_counts <- table(factor(x, levels = 1:8))
      all_probs[chain, ] <- as.numeric(state_counts / length(x))
    }
    
    # Promediar probabilidades sobre cadenas
    avg_probs <- colMeans(all_probs)
    
    predictions <- c(
      cond_prob(7, c(7, 8), avg_probs),
      cond_prob(5, c(5, 6), avg_probs),
      cond_prob(c(5, 7), c(5, 7, 6, 8), avg_probs),
      cond_prob(3, c(3, 4), avg_probs),
      cond_prob(1, c(1, 2), avg_probs),
      cond_prob(c(1, 3), c(1, 3, 2, 4), avg_probs),
      cond_prob(c(3, 7), c(3, 7, 4, 8), avg_probs),
      cond_prob(c(1, 5), c(1, 5, 2, 6), avg_probs)
    )
    
    mse <- mean((data$resp - predictions)^2, na.rm = TRUE)
    
    if (is.na(mse) || is.nan(mse) || !is.finite(mse)) {
      return(1e10)
    }
    
    return(mse)
    
  }, error = function(e) {
    return(1e10)
  }) # MSE entre predicciones y datos
}

########################
# GENERACIÓN DE DATOS SINTÉTICOS

# Generar datos sintéticos usando el Mutation Sampler
generate_synthetic_data <- function(n_participants, true_c, true_m, true_b, true_len,
                                    dag_structure = "chain", noise_sd = 0.05, 
                                    n_samples = 20) {
  # n_participants:                Número de participantes a simular
  # true_c,true_m,true_b,true_len: Parámetros verdaderos
  # dag_structure:                 Tipo de estructura causal
  # noise_sd:                      Desviación estándar del ruido gaussiano
  # n_samples:                     Número de muestras MCMC para promediar
  
  
  all_data <- data.frame()
  
  for (p in 1:n_participants) {
    # Promediar múltiples corridas del Mutation Sampler para datos más estables
    all_preds <- matrix(0, nrow = n_samples, ncol = 8)
    
    for (s in 1:n_samples) {
      pred <- PredMutationSampler(len = true_len, c = true_c, m = true_m, b = true_b,
                                  dag_structure = dag_structure, use_poisson = TRUE)
      all_preds[s, ] <- pred$resp
    }
    
    # Promediar las predicciones
    avg_resp <- colMeans(all_preds)
    
    # Crear dataframe del participante
    pred_df <- data.frame(
      cases = c("X2=1|Y=0,X1=0", "X2=1|Y=1,X1=0", "X2=1|X1=0", 
                "X2=1|Y=0,X1=1", "X2=1|Y=1,X1=1", "X2=1|X1=1", 
                "X2=1|Y=0", "X2=1|Y=1"),
      resp = avg_resp
    )
    
    # Agregar ruido gaussiano (acotado entre 0 y 1)
    pred_df$resp <- pmax(0, pmin(1, pred_df$resp + rnorm(8, 0, noise_sd)))
    
    # Agregar ID del participante
    pred_df$participant_id <- paste0("synth_", p)
    
    all_data <- rbind(all_data, pred_df)
  }
  
  return(all_data) # Dataframe con datos sintéticos
}

########################
# FUNCIONES DE FITTING CON GRID SEARCH

# Fittear modelo Normativo a un participante
fit_normative_participant <- function(df, dag_structure = "chain") {
  # df:            Dataframe con datos del participante
  # dag_structure: Tipo de estructura causal
  
  grid <- expand.grid(
    c = c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
    m = c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99),
    b = c(0.01, 0.2, 0.4, 0.6, 0.8, 0.99)
  )
  
  grid$mse <- apply(grid, 1, function(par) {
    tryCatch({
      FitNormative(par = par, dag_structure = dag_structure, data = df)
    }, error = function(e) Inf)
  })
  
  best_start <- as.numeric(grid[which.min(grid$mse), 1:3])
  
  # Usar optim con L-BFGS-B para optimización con límites
  result <- optim(
    par = best_start, 
    fn = FitNormative,
    method = "L-BFGS-B",
    lower = c(0.01, 0.01, 0.01), 
    upper = c(0.99, 0.99, 0.99),
    dag_structure = dag_structure,
    data = df
  )
  
  return(list(
    c = result$par[1],
    m = result$par[2],
    b = result$par[3],
    mse = result$value
  )) # Lista con parámetros óptimos y MSE
}

#' Fittear Mutation Sampler a un participante
fit_ms_participant <- function(df, dag_structure = "chain") {
  # df:            Dataframe con datos del participante
  # dag_structure: Tipo de estructura causal
  
  # Grid más pequeño para búsqueda inicial rápida
  grid <- expand.grid(
    c = c(0.2, 0.5, 0.8),
    m = c(0.2, 0.5, 0.8),
    b = c(0.1, 0.3, 0.5),
    len = c(8, 16, 32)
  )
  
  grid$mse <- apply(grid, 1, function(par) {
    tryCatch({
      FitMutationSampler(par = par, dag_structure = dag_structure, data = df,
                         use_poisson = FALSE, n_chains = 3)
    }, error = function(e) Inf)
  })
  
  best_start <- as.numeric(grid[which.min(grid$mse), 1:4])
  
  objective_fn <- function(par) {
    val <- FitMutationSampler(par, dag_structure = dag_structure, data = df,
                              use_poisson = FALSE, n_chains = 5)
    if (length(val) != 1 || is.na(val) || is.nan(val) || !is.finite(val)) {
      return(1e10)
    }
    return(as.numeric(val))
  }
  
  result <- tryCatch({
    optim(
      par = best_start, 
      fn = objective_fn,
      method = "L-BFGS-B",
      lower = c(0.01, 0.01, 0.01, 3), 
      upper = c(0.99, 0.99, 0.99, 64),
      control = list(maxit = 50)
    )
  }, error = function(e) {
    list(par = best_start, value = min(grid$mse))
  })
  
  return(list(
    c = result$par[1],
    m = result$par[2],
    b = result$par[3],
    len = result$par[4],
    mse = result$value
  )) # Lista con parámetros óptimos y MSE
}

########################
# EJECUCIÓN PRINCIPAL

# Parámetros para la generación de datos
TRUE_C <- 0.5         # Probabilidad base
TRUE_M <- 0.8         # Fuerza causal
TRUE_B <- 0.1         # Probabilidad de fondo
TRUE_LEN <- 12        # Largo promedio de la cadena MCMC
N_PARTICIPANTS <- 10  # Reducido para demostración más rápida
DAG_STRUCTURE <- "chain"


# Generar datos sintéticos
synthetic_data <- generate_synthetic_data(
  n_participants = N_PARTICIPANTS,
  true_c = TRUE_C,
  true_m = TRUE_M,
  true_b = TRUE_B,
  true_len = TRUE_LEN,
  dag_structure = DAG_STRUCTURE,
  noise_sd = 0.03
)


########################
# FITTING DE MODELOS

# Separar por participante
participants <- unique(synthetic_data$participant_id)

results_nrm <- data.frame()
results_ms <- data.frame()

for (p in participants) {
  cat(sprintf("Procesando %s...\n", p))
  
  df_p <- synthetic_data[synthetic_data$participant_id == p, ]
  
  # Fit Normativo
  fit_nrm <- fit_normative_participant(df_p, DAG_STRUCTURE)
  results_nrm <- rbind(results_nrm, data.frame(
    participant = p,
    c = fit_nrm$c,
    m = fit_nrm$m,
    b = fit_nrm$b,
    mse = fit_nrm$mse
  ))
  
  # Fit Mutation Sampler
  fit_ms <- fit_ms_participant(df_p, DAG_STRUCTURE)
  results_ms <- rbind(results_ms, data.frame(
    participant = p,
    c = fit_ms$c,
    m = fit_ms$m,
    b = fit_ms$b,
    len = fit_ms$len,
    mse = fit_ms$mse
  ))
}

########################
# COMPARACIÓN DE MODELOS

cat("RESULTADOS: Modelo Normativo\n")
print(results_nrm)
cat(sprintf("\nMSE promedio: %.6f\n", mean(results_nrm$mse)))
cat(sprintf("Parámetros promedio: c=%.3f, m=%.3f, b=%.3f\n", 
            mean(results_nrm$c), mean(results_nrm$m), mean(results_nrm$b)))

cat("RESULTADOS: Mutation Sampler\n")
print(results_ms)
cat(sprintf("\nMSE promedio: %.6f\n", mean(results_ms$mse)))
cat(sprintf("Parámetros promedio: c=%.3f, m=%.3f, b=%.3f, len=%.1f\n", 
            mean(results_ms$c), mean(results_ms$m), mean(results_ms$b), mean(results_ms$len)))

# Calcular AIC
n_trials <- 8
k_nrm <- 3
k_ms <- 4

results_nrm$AIC <- n_trials * log(results_nrm$mse) + 2 * (k_nrm + 1)
results_ms$AIC <- n_trials * log(results_ms$mse) + 2 * (k_ms + 1)

cat("COMPARACIÓN DE MODELOS (AIC)\n")

comparison <- data.frame(
  participant = results_nrm$participant,
  AIC_NRM = results_nrm$AIC,
  AIC_MS = results_ms$AIC
)
comparison$best_model <- ifelse(comparison$AIC_NRM < comparison$AIC_MS, "Normativo", "Mutation Sampler")

print(comparison)

cat(sprintf("Modelo Normativo ganador: %d participantes\n", sum(comparison$best_model == "Normativo")))
cat(sprintf("Mutation Sampler ganador: %d participantes\n", sum(comparison$best_model == "Mutation Sampler")))
cat(sprintf("\nAIC promedio Normativo: %.2f\n", mean(comparison$AIC_NRM)))
cat(sprintf("AIC promedio MS: %.2f\n", mean(comparison$AIC_MS)))



########################
# GUARDAR RESULTADOS

write.csv(synthetic_data, "synthetic_data.csv", row.names = FALSE)
write.csv(results_nrm, "fit_results_normative.csv", row.names = FALSE)
write.csv(results_ms, "fit_results_ms.csv", row.names = FALSE)
write.csv(comparison, "model_comparison.csv", row.names = FALSE)

