###############################################################################
# make_example_outputs.R
#
# Runs the BEACON demo with short chains to generate example output files
# for the vignette. Suitable for illustrative purposes only — estimates
# will be noisier than the full-chain run. Runtime: ~3-5 minutes.
#
# To generate publication-quality outputs, run demo_aian_gdm_oklahoma.R
# (3 chains x 10,000 iterations, ~20-30 minutes).
#
# Outputs saved to output/:
#   demo_results_table.csv
#   demo_county_posteriors.csv
#   demo_map_rho_true.pdf
#   demo_map_phi.pdf
#   demo_map_kappa.pdf
###############################################################################

# Override chain settings before sourcing the demo
n_chains <- 1
n_iter   <- 1500
n_burn   <- 500
n_thin   <- 1

# Source the full demo — it reads n_chains/n_iter/n_burn/n_thin from environment
source("demo_aian_gdm_oklahoma.R")
