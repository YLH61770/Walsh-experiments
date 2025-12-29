library(mvtnorm)    # For generating multivariate normal and t-distributed data
library(quantreg)   # For quantile regression (LAD estimation)
library(doParallel) # For parallel computing
library(foreach)    # For parallel loops
library(nor1mix)    # For generating mixed normal distributions

start_time <- Sys.time()

# =============== Data generation function ===============
generate <- function(theta, N, Xtype = c("Normal", "t3"), etype = c("Normal", "Mixed", "t3")) {
  Xtype <- match.arg(Xtype); etype <- match.arg(etype)
  p <- length(theta) - 1
  intercept <- theta[1]
  beta <- theta[-1]
  
  mu <- rep(0, p)
  Sigma <- outer(1:p, 1:p, function(i, j) 0.5^abs(i - j))
  
  if (Xtype == "Normal") {
    x <- mvtnorm::rmvnorm(N, mu, Sigma)
  } else {
    x <- mvtnorm::rmvt(N, Sigma, df = 3)
  }
  
  if (etype == "Normal") {
    err <- rnorm(N)
  } else if (etype == "Mixed") {
    err <- rnorMix(N,norMix(mu=c(-2,2),sigma=c(0.5,0.5),w=c(0.5,0.5)))
  } else {
    err <- rt(N, df = 3)
  }
  
  y <- intercept + as.vector(x %*% beta) + err
  list(y = y, x = x)
}

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
fit_lad_model <- function(X_tilde, Y_tilde, theta_true) {
  methods_to_try <- c("pfn", "fnc", "fnb")
  for (method in methods_to_try) {
    fit_full <- tryCatch({
      quantreg::rq.fit(X_tilde, Y_tilde, tau = 0.5, method = method)
    }, error = function(e) NULL)
    if (!is.null(fit_full)) {
      theta_hat_full <- as.numeric(coef(fit_full))
      if (length(theta_hat_full) == length(theta_true) &&
          all(is.finite(theta_hat_full)) && !any(is.na(theta_hat_full))) {
        return(mean((theta_hat_full - theta_true)^2))
      }
    }
  }
  return(NA)
}

# =============== Subsample LAD Fitting Function (Pairwise Method) ===============
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

H <- 1000   
n <- 10000

p <- 19
#p <- 39
#p <- 59

m_values <- c(1, 10, 20, 30, 40) 
methods <- c("gamma", "pois")  

r <- 100
#r <- 150

cores <- 4  # 并行核

intercept_true <- 1.0  # 真实截距
beta_0 <- rep(c(1, 1.5), length.out = p)
theta_true <- c(intercept_true, beta_0)  # 真实参数向量

cl <- makeCluster(cores)
registerDoParallel(cl)

# ===============Full-sample ===============
pairs_matrix <- combn(n, 2, simplify = TRUE)
storage.mode(pairs_matrix) <- "integer" 

all_mse_full <- numeric(H)

all_mse_full <- foreach(h = 1:H, 
.packages = c("mvtnorm", "quantreg"),
.combine = "c",
.errorhandling = "remove") %dopar% {
  G <- generate(theta_true, n, "Normal", "Normal")
  y <- G$y
  X <- G$x
  rm(G) 
  gc()

  matrices <- build_tilde_matrix(X, y, pairs_matrix)
  X_tilde_full <- matrices$X_tilde
  Y_tilde_full <- matrices$Y_tilde
  rm(matrices, X, y) 
  gc()
  mse_result <- fit_lad_model(X_tilde_full, Y_tilde_full, theta_true)
  rm(X_tilde_full, Y_tilde_full)
  gc()
  as.numeric(mse_result)
}

stopCluster(cl)

mean_mse_full <- mean(all_mse_full, na.rm = TRUE)

results <- data.frame(
  Method = "full",
  MSE = mean_mse_full,
  H = H,
  n = n,
  p = p
)

print(results)

write.csv(results, file = "full_sample_results.csv", row.names = TRUE)

# =============== Subsample===============
cl <- makeCluster(cores)
registerDoParallel(cl)

all_results <- foreach(h = 1:H, 
.packages = c("mvtnorm", "quantreg", "nor1mix"),
.combine = "rbind",
.errorhandling = "remove") %dopar% {
  
  G <- generate(theta_true, n, "Normal", "Normal")
  y <- G$y
  X <- G$x
  experiment_results <- matrix(NA, nrow = length(m_values), ncol = length(methods))
  for (method_idx in 1:length(methods)) {
    method <- methods[method_idx]
    for (m_idx in 1:length(m_values)) {
      m <- m_values[m_idx]
      all_estimates <- matrix(NA, nrow = m, ncol = p + 1)
      for (iter in 1:m) {
        gw <- generate_weights(n, r, method = method)
        w <- gw$w
        theta_hat <- tryCatch({
          fit_subsample_lad_optimized(X, y, w)
        }, error = function(e) rep(NA, p + 1))
        all_estimates[iter, ] <- theta_hat
      }
      theta_hat_mean <- colMeans(all_estimates, na.rm = TRUE)

      if (all(!is.na(theta_hat_mean)) && length(theta_hat_mean) == length(theta_true)) {
        mse_value <- mean((theta_hat_mean - theta_true)^2)
      } else {
        mse_value <- NA
      }
      experiment_results[m_idx, method_idx] <- mse_value
    }
  }
  experiment_results
}

stopCluster(cl)

if (is.null(dim(all_results))) {
  all_results <- matrix(all_results, nrow = 1)
}

results <- matrix(NA, nrow = length(m_values), ncol = length(methods))

for (i in 1:length(m_values)) {
  result_1 <- matrix(NA, nrow = H, ncol = length(methods))
  for (j in 1:H) {
    result_1[j, ] <- all_results[i+length(m_values)*(j-1), ] 
  }
  results[i, ] <- colMeans(result_1, na.rm = TRUE)
}
rownames(results) <- paste0("m=", m_values)

colnames(results) <- methods

print(results)

end_time <- Sys.time()
time_taken <- end_time - start_time
cat("\nTotal run time:", format(time_taken), "\n")

write.csv(results, file = "Experiment_2_p=20_1000.csv", row.names = TRUE)