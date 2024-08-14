# Objective: Use Matilda to generate an ensemble of Hector results to be used
# in the Hector-Facts workflow.

# 0. Set Up --------------------------------------------------------------------
library(matilda)
library(dplyr)
remotes::install_github("jgcri/hector@v3.2.0")
devtools::load_all("/Users/dorh012/Documents/Hector-WD/matilda")
stopifnot(packageVersion("hector") == "3.2.0")

# Set the seed to be reproducible
set.seed(42)

# Set up locations for output
BASEDIR <- here::here()
TMP_OUTDIR <- file.path(BASEDIR, "TEMP")
DATADIR <- file.path(BASEDIR, "data")

dir.create(TMP_OUTDIR, showWarnings = FALSE)
dir.create(DATADIR, showWarnings = FALSE)

# 1. Generate parameters, run Hector, & score results  -------------------------

ini <- system.file("input/hector_ssp245.ini", package = "hector")
hc <- newcore(ini)

param_values <- generate_params(hc, draws = 5000)

# Run Hector repeatedly over all parameter values
rslts <- iterate_model(core = hc, params = param_values)
write.csv(rslts, file = file.path(TMP_OUTDIR, "all_hector_rslts.csv"),
          row.names = FALSE)

# Score Hector runs with observed CO2 data
scores <- score_runs(rslts, criterion_co2_obs(), score_ramp, w1 = 2, w2 = 20)

# TODO address this!
# This cut off was honestly chosen some what arbitrarily we should revisit this!
cut_off <- 0.00020
to_keep <- which(scores$weights >= cut_off)

final_params <- param_values[to_keep, ]

# 2. Save Results  -------------------------------------------------------------

write.csv(final_params, file = file.path(DATADIR, "hector_params.csv"),
          row.names = FALSE)
