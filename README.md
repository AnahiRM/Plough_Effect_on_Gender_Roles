# Replication: The Effect of Plough Agriculture on Gender Roles

This repository contains the replication code for the analysis of the paper:

> Baiardi, Donatella and Naghi, Alireza (2024)  
> *The Effect of Plough Agriculture on Gender Roles*

The goal of this project is to reproduce the main empirical results using a customized Double Machine Learning (DML) framework implemented in R.

---

## 📁 Repository Structure

```
├── admin/          # Reference materials (paper, appendix, citations)
├── data/           # Input datasets (cross-country data)
├── functions/      # Custom R functions (ML + moment functions)
├── outputs/        # Replication results (tables and saved objects)
├── scripts/        # Main scripts to run the analysis
└── plough_replication.Rproj
```

### Folder details

- **admin/**  
  Contains the original paper, appendix, and citation files.

- **data/**  
  Includes the cross-country dataset used in the analysis (`.csv` and `.dta` formats).

- **functions/**  
  Core implementation of the DML procedure:
  - `ML_Functions_01.R`: Machine learning learners (Lasso, RF, Boosting, etc.)
  - `Moment_Functions_01.R`: Moment estimation functions
  - `Moment_functions_05_v01_adjusted2.R`: Adjusted moment functions for IV estimation

- **scripts/**  
  Main scripts to reproduce results:
  - `master.R`: Main entry point (reproduces main tables)
  - `plough_t2_*.R`: Table 2 estimations
  - `plough_t4_*.R`: Additional specifications

- **outputs/**  
  Contains estimation results:
  - `.txt` files: table outputs
  - `.RData` files: saved model objects and intermediate results

---

## ⚙️ Methodology

The replication implements a **Double Machine Learning (DML)** framework with the following components:

- Cross-fitting (2-fold, repeated 100 times)
- Estimation of nuisance functions using multiple learners:
  - Lasso
  - Regression Trees
  - Random Forest
  - Boosting
  - Neural Networks
- Aggregation strategies:
  - **Best**: selects the best-performing learner
  - **Ensemble**: weighted combination of learners

The final causal parameters are estimated using:
- OLS (Partially Linear DML)
- 2SLS (DML-IV)

---

## ▶️ How to Run the Code

1. Open the R project: plough_replication.Rproj
2.	Run the main script: source("scripts/master.R")


## 👤 Authors

Anahí Reyes Miguel and Cynthia Francis  
MScT Data and Economics for Public Policy  
École Polytechnique
