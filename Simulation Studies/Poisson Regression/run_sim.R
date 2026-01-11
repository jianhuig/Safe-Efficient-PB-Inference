# Libraries for parallelization.
library(foreach)
library(doParallel)

# Simulation Parameters.
n_tot <- 10000 # Total data size.
prop_lab <- 0.10 # Proportion of labeled data.

# Inference model formula.
formula <- y - pred ~ x1 + x2 

# Parameters for outcome and covariate generation.
family <- "poisson" # GLM family.
sigma_y <- 1 # Outcome error variance.
rho <- 0 # Covariate correlation.
p <- 4 # Covariate dimension.
infset_ind <- 1 : (p / 2) # Indices for covariates used for inference.

# Set seed.
set.seed(1947)

# Simulation loop.
sim_result_list <- foreach(i = 1:n_sim, 
                           .combine = rbind, 
                           .packages = c("dplyr", "tidyr", "lmtest", "sandwich")) %dopar% 
  {
    
    # Simulate beta_1 and beta_2 from unit sphere S^{d-1}.
    beta <- c(unit_vec(p / 2), unit_vec(p / 2))
    
    # Simulate inference data and compute analytical form of delta variance.
    sim_dat_ob <- simple_data_gen(n = n_tot, p = p, prop_lab = prop_lab,
                                  total_var = total_var, ratio_var = ratio_var,
                                  rho = rho, 
                                  sigma_y = sigma_y, beta = beta, 
                                  predset_ind =  predset_ind, 
                                  pred_noise = pred_noise, noise_diff = noise_diff,
                                  dist_shift = dist_shift, 
                                  dist_shift_diff = dist_shift_diff,
                                  dist_shift_type = dist_shift_type,
                                  family = family)
    
    sim_dat <- sim_dat_ob$analysis_data
    
    # Run analysis.
    analysis_results <- rbind(
      classical_estimation(sim_dat, formula, family, est_type = "true"),
      classical_estimation(sim_dat, formula, family, est_type = "naive"),
      classical_estimation(sim_dat, formula, family, est_type = "classical"),
      pb_estimation(sim_dat, formula, family, est_type = "ppi"),
      pb_estimation(sim_dat, formula, family, est_type = "chen-chen"),
      pb_estimation(sim_dat, formula, family, est_type = "pdc"))
    
    # Calculate bias + coverage + add delta variance.
    beta_true <- tibble::tibble(
      term = c("(Intercept)", paste0("x", infset_ind)),
      beta_true = c(0, beta[infset_ind]))
    
    # Calculate bias + coverage.
    beta_true <- tibble::tibble(
      term = c("(Intercept)", paste0("x", infset_ind)),
      beta_true = c(0, beta[infset_ind]))
    
    method_dfs <- analysis_results %>%
      left_join(beta_true, by = "term") %>%
      mutate(bias = Estimate - beta_true,
             covered = beta_true >= Lower.CI & beta_true <= Upper.CI)
    
    print(method_dfs, n = 50)
    
    return(method_dfs)
    
  }


sim_result <- bind_rows(sim_result_list)

summary_df_all <- sim_result %>%
  group_by(Method, term)  %>%
  summarise(
    bias = mean(Estimate),
    mean_std_error = mean(Std.Error),
    mean_ci_width = 1.96 * 2 * mean(Std.Error),
    cp = mean(covered),
    .groups = "drop")  


summary_df <- summary_df_all %>%
  filter(term == "x1") %>%
  mutate(
    RE = mean_std_error /
      mean_std_error[Method == "classical"]) 

summary_df
