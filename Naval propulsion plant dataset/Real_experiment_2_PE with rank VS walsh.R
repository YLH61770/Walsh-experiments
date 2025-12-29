library(mvtnorm)    # For generating multivariate normal and t-distributed data
library(quantreg)   # For quantile regression (LAD estimation)
library(doParallel) # For parallel computing
library(foreach)    # For parallel loops
library(nor1mix)    # For generating mixed normal distributions

start_time <- Sys.time()


# =============== Random disturbance weight generation function ===============
generate_weights <- function(n, r, method) {
  q <- r / n
  u <- rbinom(n, 1, q)
  
  if (method == "exp") {
    v <- rexp(n, rate = q)
  } else if (method == "geom") {
    v <- rgeom(n, q)
  } else if (method == "beta") {
    v <- runif(n, 0, 2 / q)
  } else if (method == "pois") {
    v <- rpois(n, lambda = 1 / q)
  } else if (method == "gamma"){
    v <- rgamma(n, shape = 1 / q) 
  }
  
  w <- u * v
  list(w = w, q = q)
}
# =============== Method 1: Subsample LAD Fitting Function (Rank Regression) ===============
fit_subsample_rank <- function(X, y, w) {
  
  S <- which(w != 0)
  n_S <- length(S)
  
  if (n_S == 0) {
    return(rep(NA, ncol(X)))
  }
  
  X_S <- X[S, , drop = FALSE]
  y_S <- y[S]
  w_S <- w[S]
  pairs_matrix <- combn(n_S, 2, simplify = TRUE)
  storage.mode(pairs_matrix) <- "integer"
  n_pairs <- ncol(pairs_matrix)
  n_obs <- length(y_S)
  X_tilde <- matrix(0, nrow = n_pairs + n_obs, ncol = p)
  Y_tilde <- numeric(n_pairs + n_obs)
  
  chunk1_indices <- 1: n_pairs
  idx1 <- pairs_matrix[1,chunk1_indices]
  idx2 <- pairs_matrix[2,chunk1_indices]
  
  for (j in 1:p) {
    col_idx <- j
    X_tilde[chunk1_indices, col_idx] <- (w_S[idx1] * w_S[idx2]) * (X_S[idx1, j] - X_S[idx2, j])
  }
  Y_tilde[chunk1_indices] <- (w_S[idx1] * w_S[idx2]) * (y_S[idx1] - y_S[idx2])
  
  tryCatch({
    fit <- quantreg::rq.fit(X_tilde, Y_tilde, tau = 0.5, method = "fnb")
    return(as.numeric(coef(fit)))
  }, error = function(e) {
    return(rep(NA, p))
  })
}

# =============== Method 2:Subsample LAD Fitting Function (Pairwise Method) ===============
fit_subsample_lad_optimized <- function(X, y, w) {
  
  S <- which(w != 0)
  n_S <- length(S)
  
  if (n_S == 0) {
    return(rep(NA, ncol(X) + 1))
  }
  
  X_S <- X[S, , drop = FALSE]
  y_S <- y[S]
  w_S <- w[S]
  pairs_matrix <- combn(n_S, 2, simplify = TRUE)
  storage.mode(pairs_matrix) <- "integer"
  n_pairs <- ncol(pairs_matrix)
  n_obs <- length(y_S)
  X_tilde <- matrix(0, nrow = n_pairs + n_obs, ncol = p + 1)
  Y_tilde <- numeric(n_pairs + n_obs)
  
  chunk1_indices <- 1: n_pairs
  idx1 <- pairs_matrix[1,chunk1_indices]
  idx2 <- pairs_matrix[2,chunk1_indices]
  X_tilde[chunk1_indices, 1] <- 2 * (w_S[idx1] * w_S[idx2])
  for (j in 1:p) {
    col_idx <- j + 1
    X_tilde[chunk1_indices, col_idx] <- (w_S[idx1] * w_S[idx2]) * (X_S[idx1, j] + X_S[idx2, j])
  }
  Y_tilde[chunk1_indices] <- (w_S[idx1] * w_S[idx2]) * (y_S[idx1] + y_S[idx2])
  
  
  self_indices <- (n_pairs + 1):(n_pairs + n_obs)
  self_idx <- 1:n_obs 
  
  X_tilde[self_indices, 1] <- 2* w_S[self_idx]
  
  for (j in 1:p) {
    col_idx <- j + 1
    X_tilde[self_indices, col_idx] <- 2 * w_S[self_idx] * X_S[self_idx, j]
  }
  
  
  Y_tilde[self_indices] <- 2 * w_S[self_idx] * y_S[self_idx]
  
  
  tryCatch({
    fit <- quantreg::rq.fit(X_tilde, Y_tilde, tau = 0.5, method = "fnb")
    return(as.numeric(coef(fit)))
  }, error = function(e) {
    return(rep(NA, p + 1))
  })
}

# =============== Main Simulation Program ===============

data <- read.csv("data_train_split.csv", header = TRUE)
y <- data[, ncol(data)]  
X <- data[, 1:(ncol(data)-2)] 
X <- as.matrix(X)
y <- as.numeric(y)

data_1 <- read.csv("data_test_split.csv", header = TRUE)
y_test <- data_1[, ncol(data_1)]
X_test <- data_1[, 1:(ncol(data_1)-2)] 
y_test <- as.numeric(y_test)
X_test <- as.matrix(X_test)


H <- 1000     
n <- length(y)
p <- ncol(X)
r_values <- c(200, 250, 300, 350, 400)
methods <- c("gammaR", "poisR", "gammaW", "poisW")  
m <- 20     
cores <- 4

cl <- makeCluster(cores)
registerDoParallel(cl)


all_results <- foreach(h = 1:H, 
.packages = c("mvtnorm", "quantreg"),
.combine = "rbind",
.errorhandling = "remove") %dopar% {
  experiment_results <- matrix(NA, nrow = length(r_values), ncol = length(methods))
  for (method_idx in 1:length(methods)) {
    method <- methods[method_idx]
    for (r_idx in 1:length(r_values)) {
      r <- r_values[r_idx]
      if (method %in% c("gammaR", "poisR")) {
        all_estimates <- matrix(NA, nrow = m, ncol = p)
        for (iter in 1:m) {
          if(method == "gammaR"){
            method_1 <- "gamma"
          }else{
            method_1 <- "pois"
          }
          gw <- generate_weights(n, r, method_1)
          w <- gw$w
          theta_hat <- tryCatch({
            fit_subsample_rank(X, y, w)
          }, error = function(e) rep(NA, p))
          all_estimates[iter, ] <- theta_hat
        }
        theta_hat_mean <- colMeans(all_estimates, na.rm = TRUE)
      } else {
        all_estimates <- matrix(NA, nrow = m, ncol = p + 1)
        
        for (iter in 1:m) {
          if(method == "gammaW"){
            method_1 <- "gamma"
          }else{
            method_1 <- "pois"
          }
          gw <- generate_weights(n, r, method_1)
          w <- gw$w
          theta_hat <- tryCatch({
            fit_subsample_lad_optimized(X, y, w)
          }, error = function(e) rep(NA, p + 1))
          all_estimates[iter, ] <- theta_hat
        }
        theta_hat_mean <- colMeans(all_estimates, na.rm = TRUE)
      }                       
      if (all(!is.na(theta_hat_mean))) {
        if(method %in% c("gammaR", "poisR")){
          y_pred <- X_test %*% theta_hat_mean
          PE_value <- mean((y_test - y_pred)^2, na.rm = TRUE)
        }else{
          y_pred <- cbind(1, X_test) %*% theta_hat_mean
          PE_value <- mean((y_test - y_pred)^2, na.rm = TRUE)
        }
      } else {
        PE_value <- NA
      }
      experiment_results[r_idx, method_idx] <- PE_value
    }
  }
  experiment_results
}

stopCluster(cl)

if (is.null(dim(all_results))) {
  all_results <- matrix(all_results, nrow = 1)
}

results <- matrix(NA, nrow = length(r_values), ncol = length(methods))

for (i in 1:length(r_values)) {
  result_1 <- matrix(NA, nrow = H, ncol = length(methods))
  for (j in 1:H) {
    result_1[j, ] <- all_results[i+length(r_values)*(j-1), ] 
  }
  results[i, ] <- colMeans(result_1, na.rm = TRUE)
}
rownames(results) <- paste0("r=", r_values)

colnames(results) <- methods

print(results)

end_time <- Sys.time()
time_taken <- end_time - start_time
cat("\nTotal run time:", format(time_taken), "\n")

write.csv(results, file = "(rank)Real_PE_experiment_2_1000.csv", row.names = FALSE)
