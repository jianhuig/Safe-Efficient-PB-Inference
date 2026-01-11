# Main function to generate data and calculate delta var.
simple_data_gen <- function(n, p = 2, prop_lab = 0.10, 
                            total_var = 20, ratio_var = 1, rho = 0, 
                            sigma_y = 1, beta = c(1,1),
                            predset_ind = 1, 
                            pred_noise = FALSE, noise_diff = 0,
                            dist_shift = FALSE, dist_shift_diff = 0,
                            dist_shift_type = "x1-only",
                            family = "gaussian") {
  
  # Generate inference data.
  inference_data <- get_data(n = n , p = p,
                             total_var = total_var, ratio_var = ratio_var, 
                             rho = rho, sigma_y = sigma_y, beta = beta,
                             family = family)
  
  pred_covs <- as.matrix(inference_data[, paste0("x", predset_ind)])
  
  if (dist_shift) {
    
    if (dist_shift_type == "x1-only") {
      
      beta_shift <- c(unit_vec(p / 2), rep(0, p / 2))
      
    } else if (dist_shift_type == "x2-only") {
      
      beta_shift <- c(rep(0, p / 2), unit_vec(p / 2))
      
    } else if (dist_shift_type == "both") {
      
      beta_shift <- c(unit_vec(p / 2), unit_vec(p / 2))
      
    }
    
    pred_beta <- beta + dist_shift_diff * beta_shift
    pred <- pred_covs %*% pred_beta
    
    
  } else {
    
    pred_beta <- c(beta[predset_ind])
    
    if (pred_noise) {
      
      sigma_x2 <- total_var / 2
      sigma_x3 <- noise_diff * sigma_x2
      
      noise_covs <- MASS::mvrnorm(n, mu = rep(0, p / 2),
                                  Sigma = sigma_x3 * diag(p / 2))
      
      noise_beta <- unit_vec(p / 2)
      
      pred <- pred_covs %*% pred_beta + noise_covs %*% noise_beta 
      
    } else {
      
      pred <- pred_covs %*% pred_beta 
      
    }
    
  }
  
  if(length(pred_beta) < p & family == "poisson") {
    
    # Add intercept for poisson case.
    sigma_x <- total_var / (ratio_var + 1)
    pred <- pred  + (0.5 *  sigma_x)
    
  }
  
  if (family == "poisson") {
    
    pred <- exp(pred)
    
  }
  
  # Simulated dataset.
  analysis_data <- data.frame(inference_data, pred)
  samp_ind <- rbinom(n_tot, 1, prop_lab)
  analysis_data$set <- ifelse(samp_ind == 1, "testing", "validation")
  
  return(list(analysis_data = analysis_data))
  
}

# Function to generate (Y, X) data.
get_data <- function(n, p = 2,
                     total_var = 20, ratio_var = 1, rho = 0, 
                     sigma_y = 1, beta = c(1,1), family = "gaussian"){
  
  # Obtain covariate covariance matrix. 
  sigma_x <- total_var / (ratio_var + 1)
  sigma_w <- ratio_var * sigma_x
  
  cov_diag <- sqrt(c(rep(sigma_x, p / 2), rep(sigma_w, p / 2))) 
  
  cov_matrix <- (1 - rho) * diag(cov_diag^2) + rho * (cov_diag %o% cov_diag)
  
  # Generate covariates.
  x <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = cov_matrix)
  colnames(x) <- paste0('x', 1 : p)
  
  # Generate outcome.
  if (family == 'gaussian') {
    
    eps <- rnorm(n, sd = sigma_y) 
    lin_pred <- x %*% beta 
    y <- lin_pred + eps
    
  } else if (family == 'poisson') {
    
    y <- rpois(n, exp(x %*% beta))
    
  }
  
  return(data.frame(y = y, x))
  
}

# Function to generate from unit sphere S^{d-1}.
unit_vec <- function(d) {
  
  v <- rnorm(d)
  
  return(v / l2_norm(v))
  
}

# Function for l2 norm.
l2_norm <- function(x, squared = FALSE) {
  
  y <- sum(x^2)
  
  if (!squared) {
    
    y <- sqrt(y)
    
  }
  
  return(y)
}
