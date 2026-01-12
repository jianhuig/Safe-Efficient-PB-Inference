make_formula_ipd <- function(y, yhat, xs) {
  as.formula(paste0(y, " - ", yhat, " ~ ", paste(xs, collapse = " + ")))
}
make_formula_std <- function(y, xs) {
  as.formula(paste0(y, " ~ ", paste(xs, collapse = " + ")))
}

run_one_model <- function(xs, test) {
  
  formula <- make_formula_ipd(y, yhat, xs)
  
  # your methods
  ppi <- ipd(formula, method = "ppi", model = "ols", data = test, label = lab)
  
  fake_dat <- test %>% filter(.data[[lab]] == "labeled") %>% mutate(!!lab := "unlabeled")
  ppia <- ipd(formula, method = "ppi", model = "ols", data = rbind(test, fake_dat), label = lab)
  
  cc <- chen_chen(
    Y    = test %>% filter(.data[[lab]] == "labeled") %>% pull(.data[[y]]),
    Xlab = test %>% filter(.data[[lab]] == "labeled") %>% select(all_of(xs)),
    flab = test %>% filter(.data[[lab]] == "labeled") %>% pull(.data[[yhat]]),
    fall = test %>% pull(.data[[yhat]]),
    Xall = test %>% select(all_of(xs))
  )
  
  pspa <- ipd(formula, method = "pspa", model = "ols", data = test, label = lab)
  ppiplus <- ipd(formula, method = "ppi_plusplus", model = "ols", data = test, label = lab)
  
  std <- lm(make_formula_std(y, xs), data = test)
  
  sim_dat_tv <- test %>%
    rename(set = !!lab) %>%
    mutate(set = ifelse(set == "labeled", "testing", "validation")) %>%
    rename(pred = !!yhat)
  pdc_result <- pdc(sim_dat_tv, formula, family = "gaussian")
  
  naive <- lm(make_formula_std(yhat, xs), data = test)
  
  # ---------- build `result` EXACTLY like you did ----------
  result <- c()
  
  result <- rbind(result, t(summary(std)$coefficients[, c("Estimate", "Std. Error")]))
  result <- rbind(result, (result[2,] / result[nrow(result),]))
  
  result <- rbind(result, t(summary(ppi)$coefficients[, c("Estimate", "Std. Error")]))
  result <- rbind(result, (result[2,] / result[nrow(result),]))
  
  result <- rbind(result, t(summary(ppia)$coefficients[, c("Estimate", "Std. Error")]))
  result <- rbind(result, (result[2,] / result[nrow(result),]))
  
  result <- rbind(result, t(summary(pspa)$coefficients[, c("Estimate", "Std. Error")]))
  result <- rbind(result, (result[2,] / result[nrow(result),]))
  
  result <- rbind(result, t(summary(ppiplus)$coefficients[, c("Estimate", "Std. Error")]))
  result <- rbind(result, (result[2,] / result[nrow(result),]))
  
  temp <- t(pdc_result[, c("Estimate", "Std.Error")])
  colnames(temp) <- colnames(result)          # <-- YOUR ALIGNMENT TRICK
  result <- rbind(result, temp)
  result <- rbind(result, (result[2,] / result[nrow(result),]))
  
  tempcc <- t(cc)
  colnames(tempcc) <- colnames(result)        # <-- do the same for cc (prevents NAs from name mismatch)
  result <- rbind(result, tempcc)
  result <- rbind(result, (result[2,] / result[nrow(result),]))
  
  result <- rbind(result, t(summary(naive)$coefficients[, c("Estimate", "Std. Error")]))
  result <- rbind(result, (result[2,] / result[nrow(result),]))
  
  rownames(result) <- rep(c("Estimate", "SE", "RE"), 8)
  
  result <- as.data.frame(result) %>% tibble::rownames_to_column("rowname")
  result$rowname <- rep(c("Estimate", "SE", "RE"), 8)
  
  result <- result %>%
    mutate(
      across(
        where(is.numeric),
        ~ case_when(
          rowname == "Estimate" ~ formatC(., format = "f", digits = 3),
          rowname == "SE"       ~ formatC(., format = "f", digits = 4),
          rowname == "RE"       ~ formatC(., format = "f", digits = 2),
          TRUE ~ as.character(.)
        )
      )
    )
  
  result
}