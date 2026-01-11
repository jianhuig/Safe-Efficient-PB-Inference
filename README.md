# Safe & Efficient Statistical Inference with ML-Augmented Outcome Data

This repository contains the code and results for the paper *"Safe & Efficient Statistical Inference with ML-Augmented Outcome Data"* by Jessica Gronsbell,  Jianhui Gao, Zachary R. McCaw, Yaqi Shi, and David Cheng. You can find the preprint [here](https://arxiv.org/abs/2411.19908).

## Overview

The repository includes:
- Implementation of the PB inference methods discussed in the paper, including the CC, PPI, and PDC methods.
- Code for the simulation studies.
- Code to reproduce our UK Biobank Analysis. [Note: access to UK Biobank is required as the data cannot be publicly released.]

## Repository Structure

Within the simulation folder:

* `run_sim.R`: Script for running simulations
* `simple_data_generation.R`: Contains functions for data generation
* `method_functions.R`: Contains functions for the PB inference methods methods 

## Requirements

Install the following R packages before running the methods code:

```r
install.packages(c("dplyr", "tidyr", "lmtest", "sandwich"))
```


## Example

Below is a simple demonstration of how to run the code.

```r
# Load analysis functions.
source('method_functions.R')

# Read in example data for linear regression.
analysis_data <- read.csv('example_data.csv', row.names = 1)

# Quick peak at the data.
head(analysis_data)

# Specify the model formula and GLM family.
formula <- y - pred ~ x1 + x2 + x3 + x4 + x5 
family <- "gaussian"

# Run the analysis. 
analysis_results <- rbind(
  classical_estimation(analysis_data, formula, family, est_type = "classical"),
  pb_estimation(analysis_data, formula, family, est_type = "ppi"),
  pb_estimation(analysis_data, formula, family, est_type = "chen-chen"),
  pb_estimation(analysis_data, formula, family, est_type = "pdc"))
```

You will obtain the following output. 

```r
# Review results for coefficient for x1. 
analysis_results %>% filter(term == "x1")
# A tibble: 4 × 6
  Estimate Std.Error Lower.CI Upper.CI Method    term 
*    <dbl>     <dbl>    <dbl>    <dbl> <chr>     <chr>
1   -0.170    0.0324   -0.233  -0.106  classical x1   
2   -0.169    0.0318   -0.231  -0.107  ppi       x1   
3   -0.169    0.0275   -0.223  -0.115  chen-chen x1   
4   -0.148    0.0301   -0.207  -0.0888 pdc       x1   
```

## Contact

For questions, please contact Jesse Gronsbell or open an issue on this repository.


