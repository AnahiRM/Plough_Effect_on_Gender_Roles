###########################################################################
# MASTER SCRIPT: On the origins of gender roles: Women and the plough
# Replication of Baiardi & Naghi (2024) using Double Machine Learning

# Replicators: Cynthia Francis & Anahí Reyes Miguel
###########################################################################

#  As a benchmark in the results tables: There are "OLS" and "2SLS" listed as the final column in several of the tables. However, 
# the table notes clarify that the authors are not running new, standalone OLS or 2SLS models for these columns. Instead, they are reporting 
# the original paper's (Alesina et al., 2013) baseline OLS and 2SLS estimates so that you can directly compare them side-by-side against 
# the new Double Machine Learning estimates.


cat("--- Starting Tables Replication ---\n")

# --- TABLE 1: MAIN RESULTS (Baseline controls, 3 different outcomes) ---
cat("Processing Table 1 Panel A...\n")
# Generates Table 1 Panel A: Outcome is Female labor force participation (2000).
source("scripts/plough_t4_1.R") 

cat("Processing Table 1 Panel B...\n")
# Generates Table 1 Panel B: Outcome is Share of firms with female ownership.
source("scripts/plough_t4_3.R") 

cat("Processing Table 1 Panel C...\n")
# Generates Table 1 Panel C: Outcome is Share of political positions held by women.
source("scripts/plough_t4_5.R") 


# --- TABLE 2: FULL CONTROLS & DML-IV (Outcome: Female labor force participation) ---
cat("Processing Table 2 Panel A...\n")
# Generates Table 2 Panel A: The Full Model specification (controls for 36 raw covariates).
source("scripts/plough_t2_pA.R")

cat("Processing Table 2 Panel B...\n")
# Generates Table 2 Panel B: The DML-IV model using geo-climatic controls and instruments.
source("scripts/plough_t2_2.R")



cat("\n--- Replication Complete ---\n")