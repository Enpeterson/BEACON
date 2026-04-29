###############################################################################
# BEACON Demo: Synthetic Data Generation and Model Fitting
#
# Peterson et al. (2025) -- Journal of the Royal Statistical Society Series C
#
# Purpose:
#   Demonstrates the full BEACON pipeline on synthetic data that mirrors the
#   structure of the Oklahoma AI/AN GDM application. Cherokee Nation clinical
#   data cannot be shared publicly due to Tribal data governance; synthetic
#   equivalents are generated here using parameters fixed at the posterior
#   means from the fitted model.
#
#   To obtain the Cherokee Nation data used in the paper, contact:
#   Cherokee Nation Public Health (cherokee.org/services/health)
#   A data use agreement with CNPH is required.
#
# Runtime: ~20-30 minutes (3 chains x 10,000 iterations)
#          Adjust n_chains / n_iter in Section 7 to shorten.
#
# Outputs (saved to output/):
#   demo_results_table.csv   -- posterior summaries + convergence diagnostics
#   demo_map_rho_true.pdf    -- adjusted AI/AN GDM prevalence
#   demo_map_phi.pdf         -- OSDH sensitivity surface
#   demo_map_kappa.pdf       -- CN detection/coverage probability
###############################################################################

rm(list = ls())
set.seed(2025)

###############################################################################
# 1. Libraries
###############################################################################

library(nimble)
library(coda)
library(tidyverse)
library(tigris)
library(sf)
library(spdep)
library(MASS)

options(tigris_use_cache = TRUE)

###############################################################################
# 2. Oklahoma County Spatial Structure
###############################################################################

cat("Loading Oklahoma county spatial structure...\n")

OKcounties <- tigris::counties("Oklahoma", cb = FALSE)
OKcounties <- OKcounties[order(OKcounties$NAME), ]

county_names <- OKcounties$NAME
C            <- length(county_names)   # 77 counties

Co_nb  <- poly2nb(OKcounties)
nbWB_A <- nb2WB(nb = Co_nb)

ok_map <- tigris::counties(state = "OK", cb = TRUE, class = "sf")
ok_map <- ok_map[order(ok_map$NAME), ]

###############################################################################
# 3. Synthetic Covariates
#
#   x_c: covariates for true prevalence model (SNAP, food access, prenatal care)
#   w_c: covariate for sensitivity model (log proportion AIAN)
#   d_c: covariates for CN coverage model (CN clinic count, distance-decay spillover)
###############################################################################

cat("Generating synthetic covariates...\n")

# x_c -- calibrated to Oklahoma county-level distributions
propSNAP       <- rbeta(C, 2,  8)   # proportion receiving SNAP (~mean 0.20)
prop_lowaccess <- rbeta(C, 2, 10)   # proportion with low food access (~mean 0.17)
adequate_care  <- rbeta(C, 7,  2)   # proportion with adequate prenatal care (~mean 0.78)

# w_c -- log(proportion AIAN); bound away from zero for log transform
prop_aian <- pmax(rbeta(C, 1, 18), 0.001)   # ~mean 0.05, right-skewed

# d_c -- Cherokee Nation clinic count + distance-decay spillover
# CN service counties (7 counties with CN clinical data, publicly known)
cn_county_names <- c("Adair", "Cherokee", "Delaware", "Mayes",
                     "Muskogee", "Rogers", "Sequoyah")
cn_idx    <- which(county_names %in% cn_county_names)
n_clinics <- rep(0L, C)
n_clinics[cn_idx] <- sample(1:2, length(cn_idx), replace = TRUE)

# Distance-decay CN spillover covariate.
# Publicly known Cherokee Nation Health Services facility coordinates
# (cherokee.org/services/health; no data use agreement required for locations).
# Weight = exp(-distance_km / 75); own-county clinics excluded (captured by n_clinics).
cn_clinic_coords <- data.frame(
  lat    = c(35.9106, 35.9098, 35.7297, 36.4212, 36.7480,
             35.4587, 35.8144, 36.2929, 36.7006, 36.6408),
  long   = c(-94.9710, -94.9478, -95.3104, -94.8043, -95.9724,
             -94.7879, -94.6286, -95.1533, -95.6380, -95.1564),
  county = c("Cherokee", "Cherokee", "Muskogee", "Delaware", "Washington",
             "Sequoyah", "Adair", "Mayes", "Nowata", "Craig")
)

lambda_km           <- 75
cn_clinics_sf       <- st_as_sf(cn_clinic_coords, coords = c("long", "lat"), crs = 4326) %>%
  st_transform(crs = 5070)
county_centroids_sf <- ok_map %>%
  st_transform(crs = 5070) %>%
  st_centroid() %>%
  arrange(NAME)

dist_km          <- matrix(as.numeric(st_distance(county_centroids_sf, cn_clinics_sf)),
                           nrow = C) / 1000
same_county_mask <- outer(county_centroids_sf$NAME, cn_clinic_coords$county, FUN = "==")
decay_matrix     <- exp(-dist_km / lambda_km)
decay_matrix[same_county_mask] <- 0
cn_spillover     <- rowSums(decay_matrix)

###############################################################################
# 4. Simulate Synthetic Data
#
#   True parameters fixed at posterior means from Peterson et al. (2025).
#   Spatial random effects redrawn from their respective priors.
###############################################################################

cat("Simulating synthetic data...\n")

# --- True parameter values ---
beta0  <- -1.56;  beta1  <-  3.50;  beta2  <- -0.36;  beta3  <- -0.36
gamma0 <- -0.46;  gamma1 <- -0.36
alpha0 <- -1.48;  alpha1 <-  0.97;  alpha3 <-  1.00

sigma_p_v   <- 0.345
sigma_p_u   <- 0.321
sigma_phi_v <- 0.540
sigma_phi_u <- 0.288

# --- Helper: draw from ICAR prior ---
sim_icar <- function(nb, sigma, C) {
  W     <- nb2mat(nb, style = "B", zero.policy = TRUE)
  Q     <- diag(rowSums(W)) - W + diag(C) * 1e-5   # small ridge for invertibility
  Sigma <- solve(Q) * sigma^2
  v     <- drop(MASS::mvrnorm(1, mu = rep(0, C), Sigma = Sigma))
  v - mean(v)   # sum-to-zero identifiability constraint
}

# --- Draw spatial random effects ---
v_p_c     <- sim_icar(Co_nb, sigma_p_v,     C)
u_p_c     <- rnorm(C, 0, sigma_p_u)
v_phi_c   <- sim_icar(Co_nb, sigma_phi_v,   C)
u_phi_c   <- rnorm(C, 0, sigma_phi_u)
# --- Compute latent surfaces ---
rho_true_c <- plogis(beta0 + beta1 * propSNAP + beta2 * prop_lowaccess +
                       beta3 * adequate_care + v_p_c + u_p_c)

phi_c      <- plogis(gamma0 + gamma1 * log(prop_aian) + v_phi_c + u_phi_c)

kappa_c    <- plogis(alpha0 + alpha1 * n_clinics + alpha3 * cn_spillover)

# --- Synthetic OSDH data (y^{OSDH}_c, n^{OSDH}_c) ---
# AI/AN births per county proportional to prop_aian; counties with < 5 suppressed
n_OSDH_all  <- pmax(0, round(prop_aian * runif(C, 3000, 8000)))
observed_idx <- which(n_OSDH_all >= 5)

n_OSDH_i <- n_OSDH_all[observed_idx]
y_OSDH_i <- rbinom(length(observed_idx), n_OSDH_i,
                   rho_true_c[observed_idx] * phi_c[observed_idx])
getc_i   <- observed_idx
N_OSDH   <- length(observed_idx)

# --- Synthetic CN data (y^{CN}_c, n^{CN}_c) for 7 CN counties ---
n_CN_k <- round(runif(length(cn_idx), 80, 400))
y_CN_k <- rbinom(length(cn_idx), n_CN_k,
                 rho_true_c[cn_idx] * kappa_c[cn_idx])
getc_k <- cn_idx
N_CN   <- length(cn_idx)

cat(sprintf("  Counties with OSDH data : %d / %d\n", N_OSDH, C))
cat(sprintf("  CN-covered counties     : %d\n", N_CN))
cat(sprintf("  Mean true prevalence    : %.3f\n", mean(rho_true_c)))
cat(sprintf("  Mean OSDH sensitivity   : %.3f\n", mean(phi_c)))
cat(sprintf("  Mean CN coverage        : %.3f\n", mean(kappa_c[cn_idx])))

###############################################################################
# 5. NIMBLE Model Definition
#    Notation matches manuscript (Peterson et al. 2025, Eqs. 1--5)
###############################################################################

myCode <- nimbleCode({

  # --- Likelihood: OSDH AI/AN GDM counts [Eq. 1] ---
  # y^{OSDH}_c ~ Binomial(n^{OSDH}_c, rho_error_c)
  for (i in 1:N_OSDH) {
    y_OSDH_i[i] ~ dbin(rho_error_c[getc_i[i]], n_OSDH_i[i])
  }

  # --- Likelihood: Cherokee Nation GDM counts [Eq. 2] ---
  # y^{CN}_c ~ Binomial(n^{CN}_c, rho_true_c * kappa_c)
  for (k in 1:N_CN) {
    y_CN_k[k] ~ dbin(rho_true_c[getc_k[k]] * kappa_c[getc_k[k]], n_CN_k[k])
  }

  # --- County-level models ---
  for (c in 1:Ncounties) {

    # Attenuated OSDH rate [Eq. 1]
    rho_error_c[c] <- rho_true_c[c] * phi_c[c]

    # Latent true prevalence [Eq. 3]:
    # logit(rho_true_c) = beta_0 + x_c'*beta + v_p_c + u_p_c
    logit(rho_true_c[c]) <- beta0 +
                            beta1 * propSNAP[c] +        # x_c: SNAP
                            beta2 * prop_lowaccess[c] +  # x_c: low food access
                            beta3 * adequate_care[c] +   # x_c: adequate prenatal care
                            v_p_c[c] + u_p_c[c]

    # OSDH sensitivity [Eq. 4]:
    # logit(phi_c) = gamma_0 + w_c'*gamma + v_phi_c + u_phi_c
    logit(phi_c[c]) <- gamma0 +
                       gamma1 * log(prop_aian[c]) +  # w_c: log(proportion AIAN)
                       v_phi_c[c] + u_phi_c[c]

    # CN detection/coverage [Eq. 5]:
    # logit(kappa_c) = alpha_0 + d_c'*alpha
    logit(kappa_c[c]) <- alpha0 +
                         alpha1 * n_clinics[c] +    # d_c: CN clinic count
                         alpha3 * cn_spillover[c]   # d_c: distance-decay spillover

    # iid random effects
    u_p_c[c]   ~ dnorm(0, sd = sigma_p_u)    # u_c^{(p)}   ~ N(0, sigma_{p,u}^2)
    u_phi_c[c] ~ dnorm(0, sd = sigma_phi_u)  # u_c^{(phi)} ~ N(0, sigma_{phi,u}^2)
  }

  # --- ICAR structured spatial random effects ---
  v_p_c[1:Ncounties]     ~ dcar_normal(adj[1:L], weights[1:L], num[1:Ncounties],  # v_c^{(p)}
                                       tau = 1 / (sigma_p_v^2),     zero_mean = 1)
  v_phi_c[1:Ncounties]   ~ dcar_normal(adj[1:L], weights[1:L], num[1:Ncounties],  # v_c^{(phi)}
                                       tau = 1 / (sigma_phi_v^2),   zero_mean = 1)
  # --- Priors: regression coefficients [Eq. 6] ---
  beta0 ~ dnorm(-2,  sd = 1)   # centered at logit(0.12), ~10-15% AIAN GDM prevalence
  beta1 ~ dnorm(0.5, sd = 5)   # beta: SNAP
  beta2 ~ dnorm(0.5, sd = 5)   # beta: low food access
  beta3 ~ dnorm(0.5, sd = 5)   # beta: adequate prenatal care

  gamma0 ~ dnorm(0, sd = 1)    # centered at logit(0.5) = 50% sensitivity
  gamma1 ~ dnorm(0, sd = 1)    # gamma: log(proportion AIAN)

  alpha0 ~ dnorm(-1,  sd = 1)
  alpha1 ~ dnorm(0,   sd = 1)  # alpha: CN clinic count
  alpha3 ~ dnorm(0,   sd = 1)  # alpha: distance-decay spillover

  # --- Priors: variance components, Half-Normal(0, 2.5) [Eq. 6] ---
  sigma_p_u   ~ T(dnorm(0, sd = 2.5), 0.0001, )
  sigma_p_v   ~ T(dnorm(0, sd = 2.5), 0.0001, )
  sigma_phi_u ~ T(dnorm(0, sd = 2.5), 0.0001, )
  sigma_phi_v ~ T(dnorm(0, sd = 2.5), 0.0001, )

})

###############################################################################
# 6. Model Inputs and Initial Values
###############################################################################

myData <- list(
  y_OSDH_i      = y_OSDH_i,
  n_OSDH_i      = n_OSDH_i,
  y_CN_k        = y_CN_k,
  n_CN_k        = n_CN_k,
  propSNAP      = propSNAP,
  prop_lowaccess = prop_lowaccess,
  adequate_care  = adequate_care,
  prop_aian     = prop_aian,
  n_clinics    = n_clinics,
  cn_spillover = cn_spillover
)

myConstants <- list(
  Ncounties = C,
  N_OSDH    = N_OSDH,
  N_CN      = N_CN,
  getc_i    = getc_i,
  getc_k    = getc_k,
  adj       = nbWB_A$adj,
  weights   = nbWB_A$weights,
  num       = nbWB_A$num,
  L         = length(nbWB_A$weights)
)

inits_fn <- function() {
  list(
    beta0         = rnorm(1, -2, 0.5),
    beta1         = rnorm(1,  0.5, 1),
    beta2         = rnorm(1,  0.5, 1),
    beta3         = rnorm(1,  0.5, 1),
    gamma0        = rnorm(1,  0, 0.5),
    gamma1        = rnorm(1,  0, 0.5),
    alpha0      = rnorm(1, -1, 0.5),
    alpha1      = rnorm(1,  0, 0.5),
    alpha3      = rnorm(1,  0, 0.5),
    sigma_p_u   = runif(1, 0.1, 1),
    sigma_p_v   = runif(1, 0.1, 1),
    sigma_phi_u = runif(1, 0.1, 1),
    sigma_phi_v = runif(1, 0.1, 1),
    u_p_c       = rnorm(C, 0, 0.1),
    u_phi_c     = rnorm(C, 0, 0.1),
    v_p_c       = rep(0, C),
    v_phi_c     = rep(0, C)
  )
}

###############################################################################
# 7. Fit Model
#    These are shortened chains for demonstration purposes.
#    The full analysis used 5 chains x 60,000 iterations (30,000 burn-in).
###############################################################################

if (!exists("n_chains")) n_chains <- 3
if (!exists("n_iter"))   n_iter   <- 10000
if (!exists("n_burn"))   n_burn   <- 5000
if (!exists("n_thin"))   n_thin   <- 5

params_monitor <- c(
  "beta0", "beta1", "beta2", "beta3",
  "gamma0", "gamma1",
  "alpha0", "alpha1", "alpha3",
  "sigma_p_u", "sigma_p_v",
  "sigma_phi_u", "sigma_phi_v",
  "rho_true_c", "phi_c", "kappa_c"
)

cat("\nFitting BEACON model (this takes ~20-30 minutes)...\n")

fit <- nimbleMCMC(
  code              = myCode,
  data              = myData,
  constants         = myConstants,
  inits             = lapply(1:n_chains, function(i) inits_fn()),
  monitors          = params_monitor,
  nchains           = n_chains,
  niter             = n_iter,
  nburnin           = n_burn,
  thin              = n_thin,
  progressBar       = TRUE,
  samplesAsCodaMCMC = TRUE,
  summary           = TRUE,
  WAIC              = TRUE
)

###############################################################################
# 8. Posterior Summary Table and Convergence Diagnostics
###############################################################################

cat("\nExtracting results...\n")

scalar_params <- c(
  "beta0", "beta1", "beta2", "beta3",
  "gamma0", "gamma1",
  "alpha0", "alpha1", "alpha3",
  "sigma_p_u", "sigma_p_v",
  "sigma_phi_u", "sigma_phi_v"
)

rhat_vals <- coda::gelman.diag(fit$samples[, scalar_params],
                                multivariate = FALSE)$psrf[, "Point est."]
ess_vals  <- coda::effectiveSize(fit$samples[, scalar_params])

results_table <- data.frame(
  Parameter = c(
    "beta0  (GDM intercept)",
    "beta1  (SNAP)",
    "beta2  (Low food access)",
    "beta3  (Adequate prenatal care)",
    "gamma0 (Sensitivity intercept)",
    "gamma1 (log proportion AIAN)",
    "alpha0 (CN coverage intercept)",
    "alpha1 (CN clinic count)",
    "alpha3 (distance-decay spillover)",
    "sigma_p_u   (iid SD, prevalence)",
    "sigma_p_v   (ICAR SD, prevalence)",
    "sigma_phi_u (iid SD, sensitivity)",
    "sigma_phi_v (ICAR SD, sensitivity)"
  ),
  True_Value = c(beta0, beta1, beta2, beta3,
                 gamma0, gamma1,
                 alpha0, alpha1, alpha3,
                 sigma_p_u, sigma_p_v,
                 sigma_phi_u, sigma_phi_v),
  Post_Mean  = round(fit$summary$all.chains[scalar_params, "Mean"],      3),
  Post_SD    = round(fit$summary$all.chains[scalar_params, "St.Dev."],   3),
  LCI_95     = round(fit$summary$all.chains[scalar_params, "95%CI_low"], 3),
  UCI_95     = round(fit$summary$all.chains[scalar_params, "95%CI_upp"], 3),
  Rhat       = round(rhat_vals, 3),
  ESS        = round(ess_vals,  0),
  stringsAsFactors = FALSE
)

print(results_table, row.names = FALSE)
cat(sprintf("\nWAIC: %.2f\n", fit$WAIC$WAIC))

poor_conv <- results_table[results_table$Rhat > 1.1, "Parameter"]
if (length(poor_conv) > 0) {
  warning("Rhat > 1.1 for: ", paste(poor_conv, collapse = ", "),
          "\nConsider longer chains for final analysis.")
} else {
  cat("All scalar parameters have Rhat <= 1.1.\n")
}

dir.create("output", showWarnings = FALSE)
write.csv(results_table, "output/demo_results_table.csv", row.names = FALSE)

###############################################################################
# 9. Extract County-Level Posterior Summaries
###############################################################################

get_county_post <- function(fit, prefix, C, county_names) {
  nms <- paste0(prefix, "[", 1:C, "]")
  data.frame(
    county = county_names,
    mean   = fit$summary$all.chains[nms, "Mean"],
    lci    = fit$summary$all.chains[nms, "95%CI_low"],
    uci    = fit$summary$all.chains[nms, "95%CI_upp"],
    stringsAsFactors = FALSE
  )
}

post_rho   <- get_county_post(fit, "rho_true_c", C, county_names)
post_phi   <- get_county_post(fit, "phi_c",      C, county_names)
post_kappa <- get_county_post(fit, "kappa_c",    C, county_names)

# Attach true values (available because we generated the data)
post_rho$true_value   <- rho_true_c
post_phi$true_value   <- phi_c
post_kappa$true_value <- kappa_c

# Observed (attenuated) OSDH rate for comparison
obs_df <- data.frame(county   = county_names,
                     obs_rate = NA_real_)
obs_df$obs_rate[observed_idx] <- y_OSDH_i / n_OSDH_i

# Save combined county posteriors for vignette rendering
county_post_df <- post_rho %>%
  rename(rho_mean = mean, rho_lci = lci, rho_uci = uci, rho_true = true_value) %>%
  left_join(post_phi   %>% rename(phi_mean   = mean, phi_lci   = lci, phi_uci   = uci,
                                  phi_true   = true_value), by = "county") %>%
  left_join(post_kappa %>% rename(kappa_mean = mean, kappa_lci = lci, kappa_uci = uci,
                                  kappa_true = true_value), by = "county") %>%
  left_join(obs_df, by = "county")

write.csv(county_post_df, "output/demo_county_posteriors.csv", row.names = FALSE)

###############################################################################
# 10. Maps
###############################################################################

cat("Producing maps...\n")

map_df <- ok_map %>%
  left_join(post_rho,   by = c("NAME" = "county")) %>%
  rename(rho_mean = mean, rho_lci = lci, rho_uci = uci, rho_true = true_value) %>%
  left_join(post_phi   %>% select(county, mean, true_value),
            by = c("NAME" = "county")) %>%
  rename(phi_mean = mean, phi_true = true_value) %>%
  left_join(post_kappa %>% select(county, mean, true_value),
            by = c("NAME" = "county")) %>%
  rename(kappa_mean = mean, kappa_true = true_value) %>%
  left_join(obs_df, by = c("NAME" = "county"))

# Map 1: Posterior mean adjusted prevalence (rho^{true}_c)
p_rho <- ggplot(map_df) +
  geom_sf(aes(fill = rho_mean), color = "white", linewidth = 0.2) +
  scale_fill_distiller(palette = "Purples", direction = 1,
                       name = expression(hat(rho)[c]^{true}),
                       limits = c(0, 0.5), oob = scales::squish,
                       na.value = "grey80") +
  theme_minimal(base_size = 11) +
  labs(title    = "Posterior Mean Adjusted AI/AN GDM Prevalence",
       subtitle = "BEACON model — synthetic data demonstration",
       caption  = "Grey counties: OSDH data suppressed")

# Map 2: OSDH sensitivity (phi_c)
p_phi <- ggplot(map_df) +
  geom_sf(aes(fill = phi_mean), color = "white", linewidth = 0.2) +
  scale_fill_distiller(palette = "RdYlGn", direction = 1,
                       name = expression(hat(phi)[c]),
                       limits = c(0, 1)) +
  theme_minimal(base_size = 11) +
  labs(title    = "Posterior Mean OSDH Sensitivity",
       subtitle = "BEACON model — synthetic data demonstration")

# Map 3: CN detection/coverage probability (kappa_c) for CN counties only
cn_map <- map_df %>% filter(NAME %in% cn_county_names)

p_kappa <- ggplot() +
  geom_sf(data = map_df, fill = "grey92", color = "white", linewidth = 0.2) +
  geom_sf(data = cn_map, aes(fill = kappa_mean), color = "white", linewidth = 0.2) +
  scale_fill_distiller(palette = "Blues", direction = 1,
                       name = expression(hat(kappa)[c]),
                       limits = c(0, 1)) +
  theme_minimal(base_size = 11) +
  labs(title    = "Posterior Mean CN Detection/Coverage Probability",
       subtitle = "CN-covered counties only — synthetic data demonstration")

ggsave("output/demo_map_rho_true.pdf", p_rho,   width = 8, height = 5)
ggsave("output/demo_map_phi.pdf",      p_phi,   width = 8, height = 5)
ggsave("output/demo_map_kappa.pdf",    p_kappa, width = 8, height = 5)

cat("\nDone. Outputs saved to output/:\n")
cat("  demo_results_table.csv\n")
cat("  demo_map_rho_true.pdf\n")
cat("  demo_map_phi.pdf\n")
cat("  demo_map_kappa.pdf\n")
