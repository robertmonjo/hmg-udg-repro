# ===============================
# GLOBAL single d0 fit (vectorized)
# ===============================
factor_error = 1

# -- anchors (keep your choices)
gamma_U_deg   <- 60
gamma_cen_deg <- 90
gammaU_s2     <- sin(pi*gamma_U_deg/180)^2
gammacen_s2   <- sin(pi*gamma_cen_deg/180)^2


# -----------------------
# Physical constants (SI)
# -----------------------
G_SI       <- 6.67430e-11           # m^3 kg^-1 s^-2
Msun_SI    <- 1.98847e30            # kg
c_SI       <- 2.99792458e8          # m s^-1
yr_SI      <- 365.25*24*3600        # s
Gyr_SI     <- 1e9*yr_SI
kpc_SI     <- 3.085677581e19        # m
kmps       <- 1e3                   # m/s per (km/s)

# -----------------------
# Cosmology / timescale
# -----------------------
t0_Gyr <- 13.8             # you can change this if you want
t0_SI  <- t0_Gyr * Gyr_SI
a_ct   <- c_SI / t0_SI     # c / t0 (m s^-2)


suppressPackageStartupMessages({
  library(tibble); library(dplyr); library(purrr); library(stringr)
})

udg <- tribble(
  ~AGC,   ~Distance_Mpc, ~Incl_deg, ~Rd_kpc, ~Rd_err_kpc, ~log10Mbar, ~log10Mbar_err, ~log10Mstar, ~log10Mstar_err, ~Mgas_over_Mstar, ~MHI_over_Mstar_plus, ~Mgas_over_Mstar_minus, ~Vcirc_kms, ~Vcirc_plus, ~Vcirc_minus, ~Rsys_kpc, ~vN_Rsys_kms, ~vN_Rsys_plus, ~vN_Rsys_minus,     
  114905, 76,            33,        1.79,    0.04,        9.21,       0.19,           8.30,        0.17,            7.1,               2.3,                   1.9,                    23,          4,           6,            8.02,           29,      8,     6,
  122966, 90,            34,        4.15,    0.19,        9.21,       0.13,           7.73,        0.12,            29.1,              7.0,                   11.9,                   37,          5,           6,            10.8,           25,     4,     4,
  219533, 96,            42,        2.35,    0.20,        9.36,       0.22,           8.04,        0.12,            19.7,              8.8,                   12.2,                   37,          6,           5,            9.78,           32,    12,     8,
  248945, 84,            66,        2.08,    0.07,        9.05,       0.19,           8.52,        0.17,            2.4,               0.8,                   1.6,                    27,          3,           3,            8.55,           24,     6,     5,
  334315, 73,            52,        3.76,    0.14,        9.25,       0.16,           7.93,        0.12,            23.7,              5.9,                   9.8,                    25,          5,           5,            8.49,           33,     6,     5,
  749290, 97,            39,        2.38,    0.14,        9.17,       0.15,           8.32,        0.13,            6.1,               1.7,                   2.9,                    26,          6,           6,            8.47,            27,     6,     5
)

# -- rho_vac (as you had)
rho_vac_SI <- (function(H0_km_s_Mpc = 70, Omega_L = 0.69) {
  MPC_TO_M <- 3.085677581e22
  H0_SI    <- (H0_km_s_Mpc*1000) / MPC_TO_M
  rho_crit <- 3*H0_SI^2/(8*pi*G_SI)
  Omega_L * rho_crit
})()

# -------------- Vectorized helpers --------------
vN_vec <- function(Mbar_Msun, R_kpc){
  M_kg <- Mbar_Msun * Msun_SI
  r_m  <- R_kpc * kpc_SI
  sqrt(G_SI * M_kg / r_m) / kmps
}
vH_vec <- function(R_kpc){ (R_kpc*kpc_SI / t0_SI) / kmps }

epsilon_from_dist0_vec <- function(log10Mbar, R_kpc, dist0){
  Mbar_kg <- 10^log10Mbar * Msun_SI
  R_m     <- dist0 * R_kpc * kpc_SI
  rho_nei <- (3 * Mbar_kg) / (4*pi*R_m^3)
  sqrt(pmax(rho_nei/rho_vac_SI + 1/6, 0))
}

gamma_sys_from_eps_vec <- function(eps_H, vN_kms, vH_kms){
  num  <- 2*vN_kms^2 - (eps_H^2)*vH_kms^2
  den  <- 2*vN_kms^2 + (eps_H^2)*vH_kms^2
  frac <- abs(num/den)
  sin2 <- pmin(pmax(gammaU_s2 + (gammacen_s2 - gammaU_s2)*frac, 0), 1)
  asin(sqrt(sin2))  # radians
}

vC_from_gamma0_vec <- function(vN_kms, R_kpc, gamma_sys_rad){
  # gamma0 = 1 / (cos(gamma_sys)/gamma_sys) = gamma_sys / cos(gamma_sys)
  g0    <- pmax(gamma_sys_rad, 1e-12) / pmax(cos(gamma_sys_rad), 1e-30)
  aN    <- ((vN_kms*kmps)^2) / (R_kpc*kpc_SI)            # m/s^2
  aC    <- sqrt(pmax(0, aN^2 + abs(aN) * (2*c_SI)/(g0*t0_SI)))  # m/s^2
  sqrt(aC * (R_kpc*kpc_SI)) / kmps                       # km/s
}

sigma_eff_2x <- function(residual, e_plus, e_minus){
  factor_error * ifelse(residual >= 0, pmax(e_plus, 1e-6), pmax(e_minus, 1e-6))
}

# -------------- Data vectors --------------
Rsys  <- udg$Rsys_kpc
Mbar  <- 10^udg$log10Mbar
logM  <- udg$log10Mbar
sigM  <- udg$log10Mbar_err
Vobs  <- udg$Vcirc_kms
eplus <- udg$Vcirc_plus
eminus<- udg$Vcirc_minus

# Precompute vH (independent of d0)
vH_R  <- vH_vec(Rsys)

# -------------- chi2_total(d0) --------------
chi2_total <- function(d0){
  vN   <- vN_vec(Mbar, Rsys)
  eps  <- epsilon_from_dist0_vec(logM, Rsys, d0)
  gsys <- gamma_sys_from_eps_vec(eps, vN, vH_R)
  vC   <- vC_from_gamma0_vec(vN, Rsys, gsys)
  r    <- vC - Vobs
  se   <- sigma_eff_2x(r, eplus, eminus)
  sum((r/se)^2)
}

# -------------- Optimize global d0 --------------
bounds <- c(0.1, 5e2)   # start broad
opt    <- optimize(chi2_total, interval = bounds, tol = 1e-6)
d0_hat <- opt$minimum
chi2min<- opt$objective
Ndof   <- length(Vobs) - 1
chi2_min = chi2min

# -------------- Find Δχ²=+1 bounds with geometric expansion --------------
#### pick one: #####
# DELTA_CHI2 <- 1.00      # (your current choice: 68% two-sided / ~84% one-sided)
DELTA_CHI2 <- 2.71    # 95% one-sided lower limit
# DELTA_CHI2 <- 3.84    # 95% two-sided CI

target <- chi2min + DELTA_CHI2

expand_find <- function(side = c("lo","hi"), max_expand = 60, factor = 2){
  side <- match.arg(side)
  if (side == "lo"){
    a <- bounds[1]; b <- d0_hat
    fa <- chi2_total(a) - target; fb <- chi2_total(b) - target
    # if no sign change, try expanding the left bound downward geometrically
    step <- 1
    while (fa*fb > 0 && step <= max_expand){
      a <- max(1e-6, a/factor)
      fa <- chi2_total(a) - target
      step <- step + 1
      if (a <= 1e-6) break
    }
    if (fa*fb > 0) return(NA_real_)
    uniroot(function(x) chi2_total(x) - target, interval = c(a,b))$root
  } else {
    a <- d0_hat; b <- bounds[2]
    fa <- chi2_total(a) - target; fb <- chi2_total(b) - target
    # expand right bound upward until we cross
    step <- 1
    while (fa*fb > 0 && step <= max_expand){
      b <- b*factor
      fb <- chi2_total(b) - target
      step <- step + 1
      if (b > 1e9) break
    }
    if (fa*fb > 0) return(NA_real_)
    uniroot(function(x) chi2_total(x) - target, interval = c(a,b))$root
  }
}

d0_lo <- expand_find("lo")
d0_hi <- expand_find("hi")

# -------------- Predicted vC at d0_hat (central) --------------
vN_hat   <- vN_vec(Mbar, Rsys)
eps_hat  <- epsilon_from_dist0_vec(logM, Rsys, d0_hat)
# Keep the asymptotic weak-coupling branch, but preserve vector length so the
# Monte Carlo loop remains valid even if the sample size changes.
eps_hat <- rep(sqrt(1/6), length(Vobs))
gsys_hat <- gamma_sys_from_eps_vec(eps_hat, vN_hat, vH_R)
vC_hat   <- vC_from_gamma0_vec(vN_hat, Rsys, gsys_hat)

# -------------- Predicted vC uncertainty via MC over log10 Mbar --------------
set.seed(123)
Nsamp <- 30000
vC_med  <- numeric(length(Vobs))
vC_p    <- numeric(length(Vobs))
vC_m    <- numeric(length(Vobs))

for (i in seq_along(Vobs)){
  logM_s <- rnorm(Nsamp, mean = logM[i], sd = sigM[i])
  M_s    <- 10^logM_s
  vN_s   <- vN_vec(M_s, Rsys[i])
  # (keep eps fixed at central logM for stability; if you want, recompute eps per draw)
  gsys_s <- gamma_sys_from_eps_vec(rep(eps_hat[i], Nsamp), vN_s, rep(vH_R[i], Nsamp))
  vC_s   <- vC_from_gamma0_vec(vN_s, rep(Rsys[i], Nsamp), gsys_s)
  qs     <- quantile(vC_s, c(0.16, 0.5, 0.84), na.rm = TRUE)
  vC_med[i] <- qs[2]; vC_p[i] <- qs[3]-qs[2]; vC_m[i] <- qs[2]-qs[1]
}

# -------------- Final comparison table --------------
compare_vc <- tibble(
  AGC = udg$AGC,
  Rsys_kpc = Rsys,
  Vcirc_obs_kms   = Vobs,
  Vcirc_plus_kms  = eplus,
  Vcirc_minus_kms = eminus,
  dist0_hat = d0_hat,
  vC_pred_kms      = vC_hat,
  vC_pred_plus_kms = vC_p,
  vC_pred_minus_kms= vC_m,
  residual_kms     = Vcirc_obs_kms - vC_pred_kms,
  # If you want a symmetric 2×σ for display only:
  sigma_obs_2x_sym = factor_error * 0.5 * (Vcirc_plus_kms + Vcirc_minus_kms),
  pull_sym         = residual_kms / pmax(sigma_obs_2x_sym, 1e-6)
)

# -------------- Per-galaxy chi2 contributions at d0_hat --------------
chi2_terms <- {
  r  <- vC_hat - Vobs
  se <- sigma_eff_2x(r, eplus, eminus)
  tibble(AGC = udg$AGC, residual = r, sigma_eff_2x = se, chi2_i = (r/se)^2)
}

# -------------- Summary --------------
summary_global <- data.frame(
  dist0_hat = d0_hat,
  dist0_lo  = d0_lo,
  dist0_hi  = d0_hi,
  chi2_min  = chi2min,
  Ndof      = Ndof,
  chi2_red  = chi2min / pmax(Ndof,1)
)

summary_global
compare_vc
chi2_terms
# write.csv(compare_vc, "global_dist0_compare_vc.csv", row.names = FALSE)

# Reduced chi-square and p-value
nu <- nrow(udg) - 1
chi2_red <- chi2_min / nu
pval <- 1 - pchisq(chi2_min, df = nu)
c(chi2_min = chi2_min, dof = nu, chi2_red = chi2_red, p_value = pval)

# Standardized residuals at d0_hat (sign-aware 2× sigma)
r   <- compare_vc$Vcirc_obs_kms - compare_vc$vC_pred_kms
sig <- factor_error * 0.5 * (compare_vc$Vcirc_plus_kms + compare_vc$Vcirc_minus_kms) # display-only symmetric alt.
z   <- r / pmax(sig, 1e-6)
z






# ============================================================
### MOND comparison at R_sys: v_circ^MOND vs v_N, v_C(HMG), Vcirc(obs) ####
# ============================================================

# ----- User choices -----
mond_mu <- "mls"   # options: "simple", "standard", "mls
a0_SI   <- 1.2e-10    # MOND acceleration scale (m/s^2); change if desired
Nsamp_mond <- 30000   # MC samples to propagate log10Mbar_err into vMOND

# ----- Vectorized helpers (reuse your constants) -----
vN_vec <- function(Mbar_Msun, R_kpc){
  M_kg <- Mbar_Msun * Msun_SI
  r_m  <- R_kpc * kpc_SI
  sqrt(G_SI * M_kg / r_m) / kmps
}

# --- MOND acceleration from a_N with chosen interpolating function ---
# mu options: "simple", "standard", "mls" (McGaugh–Lelli–Schombert / RAR)
mond_acc_from_aN <- function(aN, a0 = 1.2e-10, mu = c("simple","standard","mls")){
  mu <- match.arg(mu)
  aN <- pmax(aN, 0)  # safety
  
  if (mu == "simple"){
    # mu(x) = x/(1+x)  ->  a = 0.5*(aN + sqrt(aN^2 + 4 aN a0))
    0.5*(aN + sqrt(aN^2 + 4*aN*a0))
    
  } else if (mu == "standard"){
    # mu(x) = x/sqrt(1+x^2)  ->  solve for a; let A=a^2:
    # A^2 - aN^2 A - aN^2 a0^2 = 0  ->  A = 0.5*(aN^2 + sqrt(aN^4 + 4 aN^2 a0^2))
    A <- 0.5*( aN^2 + sqrt(aN^4 + 4*aN^2*a0^2) )
    sqrt(A)
    
  } else {  # "mls" (RAR-inspired): g = g_N / (1 - exp(-sqrt(g_N/a0)))
    denom <- 1 - exp(-sqrt(pmax(aN, 0)/a0))
    # guard extremely small denom to avoid blow-ups at machine precision
    denom <- pmax(denom, 1e-12)
    aN / denom
  }
}

# --- Vectorized MOND circular speed at R_sys (km/s) ---
vMOND_vec <- function(Mbar_Msun, R_kpc, mu = "mls", a0 = 1.2e-10){
  r_m <- R_kpc * kpc_SI
  # Newtonian a_N from baryons
  vN  <- {
    M_kg <- Mbar_Msun * Msun_SI
    r_mx <- R_kpc * kpc_SI
    sqrt(G_SI * M_kg / r_mx) / kmps
  }
  aN  <- (vN*kmps)^2 / r_m
  aM  <- mond_acc_from_aN(aN, a0 = a0, mu = mu)
  sqrt(aM * r_m) / kmps
}


# ----- Central MOND prediction (no MC) -----
Rsys  <- udg$Rsys_kpc
Mbar  <- 10^udg$log10Mbar
Vobs  <- udg$Vcirc_kms
eplus <- udg$Vcirc_plus
eminus<- udg$Vcirc_minus

vN_Rsys_kms   <- vN_vec(Mbar, Rsys)
vMOND_kms     <- vMOND_vec(Mbar, Rsys, mu = mond_mu, a0 = a0_SI)

#  chi2--------------

se <- ifelse(vMOND_kms - Vobs>=0, eplus, eminus)
sum(((vMOND_kms - Vobs)/se)^2)
sum(((vMOND_kms - Vobs)/(2*se))^2)

se <- ifelse(vN_Rsys_kms - Vobs>=0, eplus, eminus)
sum(((vN_Rsys_kms - Vobs)/se)^2)
sum(((vN_Rsys_kms - Vobs)/(2*se))^2)

se <- ifelse(vC_hat - Vobs>=0, eplus, eminus)
sum(((vC_hat - (Vobs+eplus))/se)^2)
sum(((vC_hat - (Vobs))/(2*se))^2)

# If you already computed HMG global prediction as 'compare_vc$vC_pred_kms', reuse it.
# Otherwise, you can omit it or set to NA here:
vHMG_kms <- if (exists("compare_vc")) compare_vc$vC_pred_kms else rep(NA_real_, length(Vobs))

# ----- Uncertainty on vMOND from log10Mbar_err (Monte Carlo) -----
set.seed(123)
logM  <- udg$log10Mbar
sigM  <- udg$log10Mbar_err

vMOND_med  <- numeric(length(Vobs))
vMOND_p    <- numeric(length(Vobs))
vMOND_m    <- numeric(length(Vobs))

for (i in seq_along(Vobs)){
  logM_s <- rnorm(Nsamp_mond, mean = logM[i], sd = sigM[i])
  M_s    <- 10^logM_s
  v_s    <- vMOND_vec(M_s, Rsys[i], mu = mond_mu, a0 = a0_SI)
  qs     <- quantile(v_s, c(0.16, 0.5, 0.84), na.rm = TRUE)
  vMOND_med[i] <- qs[2]; vMOND_p[i] <- qs[3]-qs[2]; vMOND_m[i] <- qs[2]-qs[1]
}

# ----- Build comparison table -----
compare_models <- tibble::tibble(
  AGC = udg$AGC,
  Rsys_kpc = Rsys,
  # Observed
  Vcirc_obs_kms   = Vobs,
  Vcirc_plus_kms  = eplus,
  Vcirc_minus_kms = eminus,
  # Newtonian
  vN_Rsys_kms = vN_Rsys_kms,
  # MOND (central and 1-sigma from Mbar)
  vMOND_med_kms    = vMOND_med,
  vMOND_plus_kms   = vMOND_p,
  vMOND_minus_kms  = vMOND_m,
  # HMG (if available)
  vHMG_kms = vHMG_kms
) %>%
  mutate(
    res_N   = Vcirc_obs_kms - vN_Rsys_kms,
    res_MOND= Vcirc_obs_kms - vMOND_med_kms,
    res_HMG = ifelse(is.na(vHMG_kms), NA_real_, Vcirc_obs_kms - vHMG_kms)
  )

compare_models
# write.csv(compare_models, "compare_vN_vMOND_vHMG.csv", row.names = FALSE)




# ============================================================
##### Barplot: compare v_circ (Obs, Newtonian, MOND, HMG) with errors #####
# Requires: compare_vc, compare_models, udg (as built earlier)
# ============================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(forcats)
})

# ---- Assemble a tidy table with central values + asym. errors ----
obs_df <- compare_models %>%
  transmute(
    AGC,
    model = "Observed",
    value = Vcirc_obs_kms,
    err_plus = Vcirc_plus_kms,
    err_minus = Vcirc_minus_kms
  )

newt_df <- udg %>%
  transmute(
    AGC,
    model = "Newtonian",
    value = vN_Rsys_kms,
    err_plus = coalesce(vN_Rsys_plus, 0),
    err_minus = coalesce(vN_Rsys_minus, 0)
  )

# Use MOND MC median ± (84-50, 50-16); switch to vMOND_kms if you prefer formula central
mond_df <- compare_models %>%
  transmute(
    AGC,
    model = "MOND",
    value = vMOND_med_kms,
    err_plus = vMOND_plus_kms,
    err_minus = vMOND_minus_kms
  )

hmg_df <- compare_vc %>%
  transmute(
    AGC,
    model = "HMG",
    value = vC_pred_kms,
    err_plus = vC_pred_plus_kms,
    err_minus = vC_pred_minus_kms
  )

beard_df <- bind_rows(obs_df, newt_df, mond_df, hmg_df) %>%
  filter(is.finite(value)) %>%
  mutate(
    AGC   = factor(AGC, levels = sort(unique(AGC))),
    model = factor(model, levels = c("Observed","Newtonian","MOND","HMG")),
    ymin  = pmax(0, value - err_minus),
    ymax  = value + err_plus
  )

# ---- Beardplot: 0-based bars + whiskers + center marker ----
pd <- position_dodge(width = 0.35)   # ← barbas más juntas

default_cols <- scales::hue_pal()(length(unique(beard_df$model))-1)

cols <- default_cols
cols <- c("gray30", cols)
names(cols) <- sort(unique(beard_df$model))

p_beard <- ggplot(beard_df, aes(x = AGC, y = value, fill = model, color = model)) +
  geom_linerange(aes(ymin = ymin, ymax = ymax), position = pd, linewidth = 0.9) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), position = pd, width = 0.22, linewidth = 0.6) +
  geom_point(position = pd, size = 1.7, stroke = 0.4) +
  labs(
    x = "UDG (AGC)",
    y = expression(paste("Circular speed at ", R[sys], "  (km s"^{-1},")")),
    fill = "Legend:",
    color = "Legend:",
    title = NULL
  ) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(expand = expansion(mult = c(0.14, 0.14))) +  # ← bloques más separados
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(p_beard)

ruta = "C:/Users/rober/OneDrive/Documents/Codex_Portatil/hmg_udg_article/"

# ---- Optional: export ----
ggsave(paste0(ruta,"vcirc_beardplot.png"), plot = p_beard, width = 9, height = 5, dpi = 300)
