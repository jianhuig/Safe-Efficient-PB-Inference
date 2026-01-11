# ----------- Libraries  ----------- #
library(dplyr)
library(tidyr)
library(lmtest)
library(sandwich)

#----------------------------------------------------------------------------- #
# ---------------------------- Method Functions ------------------------------ #
#----------------------------------------------------------------------------- #

# ----------- Estimation helper function  ----------- #
compute_residuals <- function(Y, X_int, theta, family){
  
  eta_hat <- X_int %*% theta
  samp_size <- nrow(X_int)
    
  if (family == "gaussian") {
    
    residuals <- c(Y - eta_hat)
    hessian <- -crossprod(X_int) / samp_size
    
  } else if (family == "poisson") {
    
    mu_hat <- exp(eta_hat)
    residuals<- c(Y - mu_hat)
    exp_deriv <- c(exp(eta_hat))
    hessian <- -t(X_int) %*% (X_int * exp_deriv) / samp_size
    
  }
  
  return(list(residuals = residuals, hessian = hessian))
  
}

# ----------- Estimation functions  ----------- #
compute_augmentation_est <- function(Y_lab, pred_lab,
                                     X_lab_int, X_all_int,
                                     beta_lab, gamma_lab, gamma_all,
                                     est_type){
  
  n <- nrow(X_lab_int)
  N <- nrow(X_all_int)
  rho <- n / N
  
  # Compute residuals and hessian for model of Y ~ X.
  beta_res <- compute_residuals(Y_lab, X_lab_int, beta_lab, family)
  
  # Compute residuals and hessian for model of Yhat ~ X.
  gamma_res <- compute_residuals(pred_lab, X_lab_int, gamma_all, family)
  
  # Obtain residuals and Hessian.
  residuals_beta <- beta_res$residuals
  D1 <- beta_res$hessian
  
  residuals_gamma <- gamma_res$residuals
  D2 <- gamma_res$hessian
  
  # Compute sandwich pieces.
  score_beta <- X_lab_int * residuals_beta
  score_gamma <- X_lab_int * residuals_gamma
  
  C11 <- crossprod(score_beta) / n
  C12 <- crossprod(score_beta, score_gamma) / n
  C22 <- crossprod(score_gamma) / n
  D1_inv <- solve(D1)
  
  # Estimate and estimated variance-covariance matrix.
  if (est_type == "chen-chen") {
    
    C22_inv <- solve(C22)
    
    D1C12 <- D1_inv %*% C12
    D1C12C22 <- D1C12 %*% C22_inv
    
    theta_hat <- as.vector(beta_lab - D1C12C22 %*% D2 %*% (gamma_lab - gamma_all))
    
    theta_hat_cov <- (D1_inv %*% C11 %*% D1_inv -
                        (1 - rho) * D1C12C22 %*% t(D1C12)) / n
    
  } else if (est_type == "ppi") {
    
    C22_inv <- solve(C22)
    
    D1C12 <- D1_inv %*% C12
    D1C12C22 <- D1C12 %*% C22_inv
    
    theta_hat <- beta_lab - (gamma_lab - gamma_all)
    
    theta_hat_cov <- ((D1_inv %*% C11 %*% D1_inv) -
                         (1 - rho) * (D1_inv %*% (2*C12 - C22) %*% D1_inv))  / n
    
  }
  
  return(list(theta_hat = theta_hat, theta_hat_cov = theta_hat_cov))
}


compute_augmented_esteq_est <- function(Y_lab, pred_lab, pred_all,
                                        X_lab_int, X_all_int,
                                        beta_lab, est_type){
  
  n <- nrow(X_lab_int)
  N <- nrow(X_all_int)
  rho <- n / N
  
  # Compute residuals and hessian for model of Y ~ X.
  beta_res <- compute_residuals(Y_lab, X_lab_int, beta_lab, family)
  
  # Compute residuals and hessian for model of Yhat ~ X, evaluated at betahat.
  gamma_res <- compute_residuals(pred_lab, X_lab_int, beta_lab, family)
  gamma_res_all <- compute_residuals(pred_all, X_all_int, beta_lab, family)
  
  # Obtain residuals and Hessian.
  residuals_beta <- beta_res$residuals
  D1 <- beta_res$hessian
  D1_inv <- solve(D1)
  
  residuals_gamma <- gamma_res$residuals
  D2 <- gamma_res$hessian
  
  # Compute scores.
  score_beta <- X_lab_int * residuals_beta
  
  score_gamma_all <- X_all_int * gamma_res_all$residuals
  score_gamma_all_int <- cbind(1, score_gamma_all)
  
  score_gamma <- X_lab_int * residuals_gamma
  score_gamma_int <- cbind(1, score_gamma)
  
  # Sandwich pieces.
  C11 <- crossprod(score_beta) / n
  C12 <- crossprod(score_beta, score_gamma) / n
  C22 <- crossprod(score_gamma) / n
  
  C22_inv <- solve(C22)
  
  D1C12 <- D1_inv %*% C12
  D1C12C22 <- D1C12 %*% C22_inv
  
  # Calculate difference in score functions.
  diff <- matrix(apply(score_gamma_int, 2, mean) - 
                   apply(score_gamma_all_int, 2, mean), 
                 ncol(score_gamma_all_int), 1)
  

  if (est_type == "pdc") {
    
    # Calculate augmentation weights.
    weights <- t(apply(score_beta, 2, function(s) {
      lm(s ~ 0 + score_gamma_int)$coefficients}))
    
    # Calculate final estimator. 
    augmentation <- -D1_inv %*% weights %*% diff
    theta_hat <- beta_lab - augmentation
    
    # Standard error. 
    theta_hat_cov <-  (D1_inv %*% C11 %*% D1_inv -
                         (1 - rho) * D1C12C22 %*% t(D1C12)) / n
    
  } else if (est_type == "ppi") {
    
    # Calculate final estimator. 
    augmentation <- (-D1_inv %*% diff[-1])
    theta_hat <- beta_lab - augmentation
    
    # Standard error. 
    theta_hat_cov <- ((D1_inv %*% C11 %*% D1_inv) -
                        (1 - rho) * (D1_inv %*% (2*C12 - C22) %*% D1_inv))  / n
    
  }
  
  return(list(theta_hat = theta_hat, theta_hat_cov = theta_hat_cov))
}


# ----------- Wrapper functions for estimation ----------- #
# --------- Computes PPI, CC, and PDC estimators --------- #
pb_estimation <- function(dat_tv,
                          formula, 
                          family,
                          est_type = "chen-chen",
                          alpha = 0.05) {
  
  # Extract response and covariates.
  response <- all.vars(formula[[2]])[1]
  rhs_vars <- attr(terms(formula), "term.labels")
  
  # Prepare data.
  labeled_data <- dplyr::filter(dat_tv, set == "testing")
  all_data <- dat_tv
  
  X_lab <- as.matrix(labeled_data[, rhs_vars])
  X_all <- as.matrix(all_data[, rhs_vars])
  X_lab_int <- cbind(1, X_lab)
  X_all_int <- cbind(1, X_all)
  
  Y_lab <- labeled_data[[response]]
  pred_lab <- labeled_data$pred
  pred_all <- all_data$pred
  
  # Model 1: Y ~ X in the labeled set.
  model_beta_lab <- glm(Y_lab ~ X_lab, family = family)
  beta_lab <- coef(model_beta_lab)
  
  # Model 2: pred ~ X in the labeled set.
  model_gamma_lab <- glm(pred_lab ~ X_lab, family = family)
  gamma_lab <- coef(model_gamma_lab)
  
  # Model 3: pred ~ X on the full set. 
  model_gamma_all <- glm(pred_all ~ X_all, family = family)
  gamma_all <- coef(model_gamma_all)
  
  if (est_type == "chen-chen") {
    
    # Compute CC estimator.
    result <- compute_augmentation_est(Y_lab, pred_lab,
                                       X_lab_int, X_all_int,
                                       beta_lab, gamma_lab, gamma_all,
                                       est_type)
    
    beta_hat <- result$theta_hat
    se_beta_hat <- sqrt(diag(result$theta_hat_cov))
      
    } else if (est_type == "ppi") {
      
      # Compute PPI estimator.
      if (family == "gaussian") {
        
        result <- compute_augmentation_est(Y_lab, pred_lab,
                                           X_lab_int, X_all_int,
                                           beta_lab, gamma_lab, gamma_all,
                                           est_type)
        
        beta_hat <- result$theta_hat
        se_beta_hat <- sqrt(diag(result$theta_hat_cov))
        
      } else if (family == "poisson") {
        
        result <- compute_augmented_esteq_est(Y_lab, pred_lab, pred_all,
                                               X_lab_int, X_all_int,
                                               beta_lab, est_type)
        
        beta_hat <- result$theta_hat
        se_beta_hat <- sqrt(diag(result$theta_hat_cov))
        
      }
        
      } else if (est_type == "pdc"){
      
      # Compute PDC estimator.
      result <- compute_augmented_esteq_est(Y_lab, pred_lab, pred_all,
                                             X_lab_int, X_all_int,
                                             beta_lab, est_type)
      
      beta_hat <- result$theta_hat
      se_beta_hat <- sqrt(diag(result$theta_hat_cov))
      
    } 
  
  # Output.
  terms <- c("(Intercept)", rhs_vars)
  critical_val <- qnorm(1 - alpha / 2)
  
  df <- data.frame(
    Estimate = beta_hat,
    Std.Error = se_beta_hat,
    Method = est_type,
    term = terms
  ) %>%
    mutate(
      Lower.CI = Estimate - critical_val * Std.Error,
      Upper.CI = Estimate + critical_val * Std.Error
    ) %>%
    dplyr::select(Estimate, Std.Error, Lower.CI, Upper.CI, Method, term)
  
  return(df)
}

# ----------- Main function for computing classical est. ----------- #
# ----------- Computes True, Naive, Classical ests. ----------- #

# Method 4: naive beta
classical_estimation <- function(dat_tv, 
                                 formula, 
                                 family, 
                                 robust = TRUE,
                                 est_type = "classical",
                                 alpha = 0.05) {
  
  # Subset data if necessary.
  if(est_type == "classical"){
    
    dat_tv_fit <- filter(dat_tv, set == "testing")
    
  } else {
    
    dat_tv_fit <- dat_tv
  }
  
  # Obtain formula.
  if(est_type == "naive"){
    
    formula <- as.formula(call("~", as.name("pred"), formula[[3]]))
    
  } else {
    
    formula <- as.formula(call("~", as.name("y"), formula[[3]]))
    
  }
  
  # Fit the model.
  fit <- glm(formula, data = dat_tv_fit, family = family)
  
  if(robust){
    
    fit <- coeftest(fit, vcov. = vcovHC(fit, type = "HC"))
    
  }
  
  # Output.
  critical_val <- qnorm(1 - alpha / 2)
  
  df <- broom::tidy(fit) %>%
    mutate(
      Method = est_type,
      conf.low = estimate - critical_val * std.error,
      conf.high = estimate + critical_val * std.error
    ) %>%
    dplyr::select(Estimate = estimate, Std.Error = std.error,
           Lower.CI = conf.low, Upper.CI = conf.high,
           Method, term)
  
  return(df)
}
