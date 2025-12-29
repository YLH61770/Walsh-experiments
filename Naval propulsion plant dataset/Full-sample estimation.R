library(mvtnorm)    # For generating multivariate normal and t-distributed data
library(quantreg)   # For quantile regression (LAD estimation)
library(doParallel) # For parallel computing
library(foreach)    # For parallel loops

start_time <- Sys.time()

# =============== Full-sample matrix construction function ===============
build_tilde_matrix <- function(X, y, pairs_matrix) {
  # Function: Construct X_tilde_full and Y_tilde_full matrices
  # Parameters: X, y, pairs_matrix
  # Returns: A list containing X_tilde_full and Y_tilde_full
  
  n_pairs <- ncol(pairs_matrix)
  n_obs <- length(y)
  p <- ncol(X)
  
  X_tilde_full <- matrix(0, nrow = n_pairs + n_obs, ncol = p + 1)
  Y_tilde_full <- numeric(n_pairs + n_obs)
  chunk1_indices <- 1: n_pairs
  idx1 <- pairs_matrix[1,chunk1_indices]
  idx2 <- pairs_matrix[2,chunk1_indices]
  X_tilde_full[chunk1_indices, 1] <- 2
  for (j in 1:p) {
    col_idx <- j + 1
    X_tilde_full[chunk1_indices, col_idx] <- X[idx1, j] + X[idx2, j]
  }
  Y_tilde_full[chunk1_indices] <- y[idx1] + y[idx2]
  
  
  self_indices <- (n_pairs + 1):(n_pairs + n_obs)
  self_idx <- 1:n_obs 
  
  X_tilde_full[self_indices, 1] <- 2
  
  for (j in 1:p) {
    col_idx <- j + 1
    X_tilde_full[self_indices, col_idx] <- 2 * X[self_idx, j]
  }
  
  
  Y_tilde_full[self_indices] <- 2 * y[self_idx]
  
  list(X_tilde = X_tilde_full, Y_tilde = Y_tilde_full)
}


# =============== Full-sample LAD fitting function ===============
fit_lad_model <- function(X_tilde, Y_tilde) {
  
  methods_to_try <- c("pfn", "fnc", "fnb")
  
  for (method in methods_to_try) {
    fit_full <- tryCatch({
      quantreg::rq.fit(X_tilde, Y_tilde, tau = 0.5, method = method)
    }, error = function(e) NULL)
    
    if (!is.null(fit_full)) {
      theta_hat_full <- as.numeric(coef(fit_full))
      
      if (all(is.finite(theta_hat_full)) && !any(is.na(theta_hat_full))) {
        return(theta_hat_full)
      }
    }
  } 
  warning("All LAD fitting methods failed")
  return(rep(NA, ncol(X_tilde)))
}


# =============== Main Simulation Program ===============

data <- read.csv("data.csv")
y <- data[, ncol(data)] 
X <- data[, 1:(ncol(data)-2)]  
y <- as.numeric(y)
X <- as.matrix(X)

n <- length(y) 
p <- ncol(X) 

pairs_matrix <- combn(n, 2, simplify = TRUE)
storage.mode(pairs_matrix) <- "integer"


matrices <- build_tilde_matrix(X, y, pairs_matrix)
X_tilde_full <- matrices$X_tilde
Y_tilde_full <- matrices$Y_tilde

theta_hat_full <- fit_lad_model(X_tilde_full, Y_tilde_full)
write.csv(theta_hat_full, "theta_hat_full_2.csv", row.names = FALSE)
