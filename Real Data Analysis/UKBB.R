library(dplyr)
library(kableExtra)
library(ipd)
library(glmnet)
source("method_functions.R")

# DEXA
dexa <- read.table("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data/DEXA.tab", header = TRUE)
ancestry <- read.table("/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/Desktop/SyntheticSurrogateAnalysis/Data/Subpopulation.tab", header = TRUE)
dexa <- dexa %>% left_join(ancestry)
dexa$f.23248.2.0 <- dexa$f.23248.2.0 / 1000

dat <- readRDS('/Users/jianhuigao/Library/CloudStorage/OneDrive-UniversityofToronto/full_data.rds')
dat <- as.data.frame(dat)
colnames(dat) <- c("f.eid", "BMD", "PA", "age", "smoke", "drink", "driving", "computer", "tv", "sex")
dat <- dat[!dat$smoke == -3,]
dat <- dat[!dat$drink == -3,]
dat <- dat[!dat$driving == -3 & !dat$driving == -1,]
dat$driving[dat$driving == -10] <- 0.5
dat <- dat[!dat$computer == -3 & !dat$computer == -1,]
dat$computer[dat$computer == -10] <- 0.5
dat <- dat[!dat$tv == -3 & !dat$tv == -1,]
dat$tv[dat$tv == -10] <- 0.5
dat$SB <- dat$driving + dat$computer + dat$tv
dat$PA <- factor(dat$PA)

dexa <- dexa %>% left_join(dat, by = c("f.eid" = "f.eid"))

pheno_id <- "f.23248.2.0" 
cov_id <- c("f.21022.0.0", "f.22001.0.0", "f.23104.0.0", "PA", "smoke", "SB", "drink", "f.50.0.0", "f.23107.0.0", "f.23108.0.0")

###############################################################################
############################# Random Forest ###################################
###############################################################################

model <- c("strong")
if(model == "strong"){
  W <- NULL
} else {
  W <- c("f.21022.0.0", "f.22001.0.0", "f.23104.0.0") # age, sex, BMI
}

ancestry_id <- c(3, 3001, 3002, 3003)

# South Asian subpopulation
test <- dexa %>%
  filter(f.21000.0.0 %in% ancestry_id) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id))

train <- dexa %>%
  filter(f.21000.0.0 %in% c(1, 1001, 1002, 1003)) %>%
  filter(!is.na(f.23248.2.0)) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id)) %>%
  select(-W)

# Train the model
set.seed(1234)
model <- ranger::ranger(f.23248.2.0 ~ ., data = train, num.trees = 1000)
train$yhat <- predict(model, train)$predictions
test$yhat <- predict(model, test)$predictions


# Inverse normal transformation
test <- INT(test, "f.23248.2.0")
test <- INT(test, "yhat")

test$f.23248.2.0 <- test$f.23248.2.0_int
test$yhat <- test$yhat_int
test$set_label <- ifelse(is.na(test$f.23248.2.0), "unlabeled", "labeled")



run_one_model(c("f.21022.0.0","f.23104.0.0"), test) %>%
  slice(1:3, 7:9, 19:21, 16:18) %>%
  kable(
    col.names = c(
      "",
      "Intercept","Age","BMI"
    ),
    align = "c"
  )  %>%
  group_rows("classical",   1,  3) %>%
  group_rows("ppi", 4,  6) %>%
  group_rows("chen-chen",    7, 9) %>%
  group_rows("pdc",   10, 12)
  
model <- c("weak")
if(model == "strong"){
  W <- NULL
} else {
  W <- c("f.21022.0.0", "f.22001.0.0", "f.23104.0.0") # age, sex, BMI
}

# South Asian subpopulation
test <- dexa %>%
  filter(f.21000.0.0 %in% ancestry_id) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id))

train <- dexa %>%
  filter(f.21000.0.0 %in% c(1, 1001, 1002, 1003)) %>%
  filter(!is.na(f.23248.2.0)) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id)) %>%
  select(-W)

# Train the model
set.seed(1234)
model <- ranger::ranger(f.23248.2.0 ~ ., data = train, num.trees = 1000)

train$yhat <- predict(model, train)$predictions


test$yhat <- predict(model, test)$predictions


# Inverse normal transformation
test <- INT(test, "f.23248.2.0")
test <- INT(test, "yhat")

test$f.23248.2.0 <- test$f.23248.2.0_int
test$yhat <- test$yhat_int

test$set_label <- ifelse(is.na(test$f.23248.2.0), "unlabeled", "labeled")

# check R squared
r_squared <- cor(test$f.23248.2.0[test$set_label == "labeled"], test$yhat[test$set_label == "labeled"])^2
print(paste("R squared in the labeled set:", round(r_squared, 4)))


run_one_model(c("f.21022.0.0","f.23104.0.0"), test) %>%
  slice(1:3, 7:9, 19:21, 16:18) %>%
  kable(
    col.names = c(
      "",
      "Intercept","Age","BMI"
    ),
    align = "c"
  )  %>%
  group_rows("classical",   1,  3) %>%
  group_rows("ppi", 4,  6) %>%
  group_rows("chen-chen",    7, 9) %>%
  group_rows("pdc",   10, 12)

###########################################################################
########################### Ridge #########################################
###########################################################################

library(glmnet)
model <- c("strong")
if(model == "strong"){
  W <- NULL
} else {
  W <- c("f.21022.0.0", "f.22001.0.0", "f.23104.0.0") # age, sex, BMI
}

# South Asian subpopulation
test <- dexa %>%
  filter(f.21000.0.0 %in% ancestry_id) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id))

train <- dexa %>%
  filter(f.21000.0.0 %in% c(1, 1001, 1002, 1003)) %>%
  filter(!is.na(f.23248.2.0)) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id)) %>%
  select(-W)

# Lasso model
set.seed(1234)
model <- cv.glmnet(
  x       = as.matrix(train %>% select(-f.23248.2.0)),
  y       = train$f.23248.2.0,
  alpha   = 0
)
train$yhat <- predict(model, as.matrix(train %>% select(-f.23248.2.0)), s = "lambda.min") %>% as.numeric()
test$yhat <- predict(model, as.matrix(test %>% select(-f.23248.2.0)), s = "lambda.min") %>% as.numeric()


# Inverse normal transformation
test <- INT(test, "f.23248.2.0")
test <- INT(test, "yhat")

test$f.23248.2.0 <- test$f.23248.2.0_int
test$yhat <- test$yhat_int
test$set_label <- ifelse(is.na(test$f.23248.2.0), "unlabeled", "labeled")

# check R squared
r_squared <- cor(test$f.23248.2.0[test$set_label == "labeled"], test$yhat[test$set_label == "labeled"])^2
print(paste("R squared in the labeled set:", round(r_squared, 4)))


run_one_model(c("f.21022.0.0","f.23104.0.0"), test) %>%
  slice(1:3, 7:9, 19:21, 16:18) %>%
  kable(
    col.names = c(
      "",
      "Intercept","Age","BMI"
    ),
    align = "c"
  )  %>%
  group_rows("classical",   1,  3) %>%
  group_rows("ppi", 4,  6) %>%
  group_rows("chen-chen",    7, 9) %>%
  group_rows("pdc",   10, 12)

model <- c("weak")
if(model == "strong"){
  W <- NULL
} else {
  W <- c("f.21022.0.0", "f.22001.0.0", "f.23104.0.0") # age, sex, BMI
}


# South Asian subpopulation
test <- dexa %>%
  filter(f.21000.0.0 %in% ancestry_id) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id))

train <- dexa %>%
  filter(f.21000.0.0 %in% c(1, 1001, 1002, 1003)) %>%
  filter(!is.na(f.23248.2.0)) %>% 
  select(all_of(c(pheno_id, cov_id))) %>% 
  tidyr::drop_na(all_of(cov_id)) %>%
  select(-W)

# Train the model
set.seed(1234)
model <- cv.glmnet(
  x       = as.matrix(train %>% select(-f.23248.2.0)),
  y       = train$f.23248.2.0,
  alpha   = 0
)
train$yhat <- predict(model, as.matrix(train %>% select(-f.23248.2.0)), s = "lambda.min") %>% as.numeric()
test$yhat <- predict(model, as.matrix(test %>% select(all_of(cov_id))  %>%
                                        select(-W)), s = "lambda.min") %>% as.numeric()

# Inverse normal transformation
test <- INT(test, "f.23248.2.0")
test <- INT(test, "yhat")

test$f.23248.2.0 <- test$f.23248.2.0_int
test$yhat <- test$yhat_int

test$set_label <- ifelse(is.na(test$f.23248.2.0), "unlabeled", "labeled")

# check R squared
r_squared <- cor(test$f.23248.2.0[test$set_label == "labeled"], test$yhat[test$set_label == "labeled"])^2
print(paste("R squared in the labeled set:", round(r_squared, 4)))

run_one_model(c("f.21022.0.0","f.23104.0.0"), test) %>%
  slice(1:3, 7:9, 19:21, 16:18) %>%
  kable(
    col.names = c(
      "",
      "Intercept","Age","BMI"
    ),
    align = "c"
  )  %>%
  group_rows("classical",   1,  3) %>%
  group_rows("ppi", 4,  6) %>%
  group_rows("chen-chen",    7, 9) %>%
  group_rows("pdc",   10, 12)
