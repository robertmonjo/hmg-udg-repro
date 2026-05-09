# Reproducible HMG / MOND / Newton comparison for gas-rich UDGs
# - base R only
# - notation aligned with the final manuscript
# - one-sided 95% lower limit on s uses Delta chi^2 = 2.71

script_path <- normalizePath(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]))
out_dir <- dirname(script_path)
if(length(out_dir)==0)
  out_dir = getwd()

gamma_U_deg   <- 60
gamma_cen_deg <- 90
gammaU_s2     <- sin(pi * gamma_U_deg / 180)^2
gammacen_s2   <- sin(pi * gamma_cen_deg / 180)^2

G_SI    <- 6.67430e-11
Msun_SI <- 1.98847e30
c_SI    <- 2.99792458e8
yr_SI   <- 365.25 * 24 * 3600
Gyr_SI  <- 1e9 * yr_SI
kpc_SI  <- 3.085677581e19
kmps    <- 1e3

t0_Gyr <- 13.8
t0_SI  <- t0_Gyr * Gyr_SI

MPC_TO_M <- 3.085677581e22
H0_km_s_Mpc <- 70
Omega_vac <- 0.70
H0_SI <- (H0_km_s_Mpc * 1000) / MPC_TO_M
rho_vac_SI <- Omega_vac * 3 * H0_SI^2 / (8 * pi * G_SI)

udg <- data.frame(
  AGC = c("114905", "122966", "219533", "248945", "334315", "749290"),
  log10Mbar = c(9.21, 9.21, 9.36, 9.05, 9.25, 9.17),
  Vobs = c(23, 37, 37, 27, 25, 26),
  eplus = c(4, 5, 6, 3, 5, 6),
  eminus = c(6, 6, 5, 3, 5, 6),
  Rsys = c(8.02, 10.80, 9.78, 8.55, 8.49, 8.47),
  vN_tab = c(29, 25, 32, 24, 30, 27),
  vN_plus = c(7, 4, 9, 6, 6, 5),
  vN_minus = c(6, 4, 7, 5, 5, 4),
  stringsAsFactors = FALSE
)

hmg_model_errors <- data.frame(
  AGC = udg$AGC,
  plus = c(5.72, 2.64, 7.05, 3.96, 4.80, 3.75),
  minus = c(4.03, 2.17, 4.79, 2.88, 3.48, 2.96)
)

mond_model_errors <- data.frame(
  AGC = udg$AGC,
  plus = c(9.50, 6.07, 12.0, 8.28, 8.21, 6.97),
  minus = c(8.15, 5.68, 10.2, 7.30, 7.04, 6.44)
)

format_pm <- function(value, plus, minus, digits = 1) {
  fmt <- paste0("%.", digits, "f")
  sprintf("$%s^{+%s}_{-%s}$", sprintf(fmt, value), sprintf(fmt, plus), sprintf(fmt, minus))
}

vN_vec <- function(Mbar_Msun, R_kpc) {
  M_kg <- Mbar_Msun * Msun_SI
  r_m  <- R_kpc * kpc_SI
  sqrt(G_SI * M_kg / r_m) / kmps
}

vH_vec <- function(R_kpc) {
  H0_km_s_Mpc * R_kpc / 1000
}

epsilon_from_s_vec <- function(log10Mbar, R_kpc, s) {
  Mbar_kg <- 10^log10Mbar * Msun_SI
  R_m     <- s * R_kpc * kpc_SI
  rho_nei <- (3 * Mbar_kg) / (4 * pi * R_m^3)
  rho_nei / rho_vac_SI + 1 / 6
}

gamma_sys_from_eps_vec <- function(eps2_H, vN_kms, vH_kms) {
  num  <- 2 * vN_kms^2 - eps2_H * vH_kms^2
  den  <- 2 * vN_kms^2 + eps2_H * vH_kms^2
  frac <- abs(num / den)
  sin2 <- pmin(pmax(gammaU_s2 + (gammacen_s2 - gammaU_s2) * frac, 0), 1)
  asin(sqrt(sin2))
}

vC_from_gamma0_vec <- function(vN_kms, R_kpc, gamma_sys_rad) {
  g0 <- pmax(gamma_sys_rad, 1e-12) / pmax(cos(gamma_sys_rad), 1e-30)
  aN <- ((vN_kms * kmps)^2) / (R_kpc * kpc_SI)
  aC <- sqrt(pmax(0, aN^2 + abs(aN) * 2 * c_SI * H0_SI / g0))
  sqrt(aC * (R_kpc * kpc_SI)) / kmps
}

vMOND_vec <- function(vN_kms, R_kpc, a0_SI = 1.2e-10) {
  r_m <- R_kpc * kpc_SI
  aN <- ((vN_kms * kmps)^2) / r_m
  x <- aN / a0_SI
  nu_hat <- 1 / (1 - exp(-sqrt(x)))
  g_mond <- nu_hat * aN
  sqrt(g_mond * r_m) / kmps
}

sigma_eff <- function(residual, e_plus, e_minus) {
  ifelse(residual >= 0, pmax(e_plus, 1e-6), pmax(e_minus, 1e-6))
}

combined_sigma_z <- function(v_model, sigma_model_plus, sigma_model_minus, v_obs, sigma_obs_plus, sigma_obs_minus) {
  if (v_model >= v_obs) {
    sigma_comb <- sqrt(sigma_obs_plus^2 + sigma_model_minus^2)
  } else {
    sigma_comb <- sqrt(sigma_obs_minus^2 + sigma_model_plus^2)
  }
  (v_model - v_obs) / sigma_comb
}

Mbar <- 10^udg$log10Mbar
Rsys <- udg$Rsys
vH_R <- vH_vec(Rsys)

chi2_hmg <- function(s) {
  vN <- udg$vN_tab
  eps2 <- epsilon_from_s_vec(udg$log10Mbar, Rsys, s)
  gsys <- gamma_sys_from_eps_vec(eps2, vN, vH_R)
  vC <- vC_from_gamma0_vec(vN, Rsys, gsys)
  r <- vC - udg$Vobs
  se <- sigma_eff(r, udg$eplus, udg$eminus)
  sum((r / se)^2)
}

chi2_fixed <- function(model_name) {
  if (model_name == "newton") {
    v_model <- udg$vN_tab
  } else if (model_name == "mond") {
    vN <- vN_vec(Mbar, Rsys)
    v_model <- vMOND_vec(vN, Rsys)
  } else {
    stop("Unknown model")
  }
  r <- v_model - udg$Vobs
  se <- sigma_eff(r, udg$eplus, udg$eminus)
  sum((r / se)^2)
}

opt <- optimize(chi2_hmg, interval = c(1, 1e6), tol = 1e-6)
s_hat <- opt$minimum
chi2_min <- opt$objective
target_95 <- chi2_min + 2.71

find_lo <- function() {
  a <- 1
  b <- s_hat
  fa <- chi2_hmg(a) - target_95
  fb <- chi2_hmg(b) - target_95
  step <- 1
  while (fa * fb > 0 && step <= 60) {
    a <- max(1e-6, a / 2)
    fa <- chi2_hmg(a) - target_95
    step <- step + 1
    if (a <= 1e-6) break
  }
  if (fa * fb > 0) return(NA_real_)
  uniroot(function(x) chi2_hmg(x) - target_95, interval = c(a, b))$root
}

s95_lo <- find_lo()

vN_asym <- udg$vN_tab
eps2_asym <- rep(1 / 6, nrow(udg))
gamma_asym <- gamma_sys_from_eps_vec(eps2_asym, vN_asym, vH_R)
vHMG_asym <- vC_from_gamma0_vec(vN_asym, Rsys, gamma_asym)
vMOND <- vMOND_vec(vN_vec(Mbar, Rsys), Rsys)

z_newton <- mapply(
  combined_sigma_z,
  v_model = udg$vN_tab,
  sigma_model_plus = udg$vN_plus,
  sigma_model_minus = udg$vN_minus,
  v_obs = udg$Vobs,
  sigma_obs_plus = udg$eplus,
  sigma_obs_minus = udg$eminus
)

z_mond <- mapply(
  combined_sigma_z,
  v_model = vMOND,
  sigma_model_plus = mond_model_errors$plus,
  sigma_model_minus = mond_model_errors$minus,
  v_obs = udg$Vobs,
  sigma_obs_plus = udg$eplus,
  sigma_obs_minus = udg$eminus
)

z_hmg <- mapply(
  combined_sigma_z,
  v_model = vHMG_asym,
  sigma_model_plus = hmg_model_errors$plus,
  sigma_model_minus = hmg_model_errors$minus,
  v_obs = udg$Vobs,
  sigma_obs_plus = udg$eplus,
  sigma_obs_minus = udg$eminus
)

cat("# Global fit summary\n")
cat(sprintf("s_hat_at_scan_edge = %.6f\n", s_hat))
cat(sprintf("chi2_hmg_asymptotic = %.6f\n", chi2_min))
cat(sprintf("s_95_one_sided_lower = %.6f\n", s95_lo))
cat(sprintf("chi2_newton = %.6f\n", chi2_fixed("newton")))
cat(sprintf("chi2_mond_mls = %.6f\n\n", chi2_fixed("mond")))

cat("# Per-galaxy velocities at the asymptotic branch\n")
cat("AGC\tv_obs\tv_newton\tv_mond\tv_hmg_inf\n")
for (i in seq_len(nrow(udg))) {
  cat(sprintf(
    "%s\t%.3f\t%.3f\t%.3f\t%.3f\n",
    udg$AGC[i], udg$Vobs[i], vN_asym[i], vMOND[i], vHMG_asym[i]
  ))
}
cat("\n")

cat("# Combined tensions z = (model - obs) / sigma_comb\n")
cat("AGC\tz_newton\tz_mond\tz_hmg\n")
for (i in seq_len(nrow(udg))) {
  cat(sprintf("%s\t%.2f\t%.2f\t%.2f\n", udg$AGC[i], z_newton[i], z_mond[i], z_hmg[i]))
}
cat("\n")

cat("# Combined squared tensions\n")
cat(sprintf("chi2_comb_newton = %.6f\n", sum(z_newton^2)))
cat(sprintf("chi2_comb_hmg = %.6f\n", sum(z_hmg^2)))
cat(sprintf("chi2_comb_mond = %.6f\n", sum(z_mond^2)))

table2 <- data.frame(
  AGC = udg$AGC,
  v_obs = udg$Vobs,
  v_obs_plus = udg$eplus,
  v_obs_minus = udg$eminus,
  v_newton = udg$vN_tab,
  v_newton_plus = udg$vN_plus,
  v_newton_minus = udg$vN_minus,
  z_newton = z_newton,
  v_mond = vMOND,
  v_mond_plus = mond_model_errors$plus,
  v_mond_minus = mond_model_errors$minus,
  z_mond = z_mond,
  v_hmg = vHMG_asym,
  v_hmg_plus = hmg_model_errors$plus,
  v_hmg_minus = hmg_model_errors$minus,
  z_hmg = z_hmg,
  stringsAsFactors = FALSE
)

write.csv(table2, file.path(out_dir, "table2_repro.csv"), row.names = FALSE)

table_lines <- c(
  "% Reproducible Table 2 generated by fit_hmg_udg.R",
  "\\begin{tabular}{lccccccc}",
  "\\toprule",
  "Galaxy & $v_{\\rm obs}$ & $v_\\text{\\tiny N }$ & $z_\\text{\\tiny N }$ & $v_{\\rm MOND}$ & $z_{\\rm MOND}$ & $v_{\\rm HMG,\\infty}$ & $z_{\\rm HMG,\\infty}$ \\\\",
  " & (km s$^{-1}$) & (km s$^{-1}$) & ($\\sigma$) & (km s$^{-1}$) & ($\\sigma$) & (km s$^{-1}$) & ($\\sigma$) \\\\",
  "\\midrule"
)
for (i in seq_len(nrow(table2))) {
  table_lines <- c(
    table_lines,
    sprintf(
      "%s & %s & %s & %.2f & %s & %.2f & %s & %.2f \\\\",
      table2$AGC[i],
      format_pm(table2$v_obs[i], table2$v_obs_plus[i], table2$v_obs_minus[i], 0),
      format_pm(table2$v_newton[i], table2$v_newton_plus[i], table2$v_newton_minus[i], 0),
      table2$z_newton[i],
      format_pm(table2$v_mond[i], table2$v_mond_plus[i], table2$v_mond_minus[i], 1),
      table2$z_mond[i],
      format_pm(table2$v_hmg[i], table2$v_hmg_plus[i], table2$v_hmg_minus[i], 1),
      table2$z_hmg[i]
    )
  )
}
table_lines <- c(table_lines, "\\bottomrule", "\\end{tabular}", "")
writeLines(table_lines, file.path(out_dir, "table2_repro.tex"))

svg(file.path(out_dir, "vcirc_beardplot_repro.svg"), width = 11, height = 6.2)
par(mar = c(5, 5, 3, 1))
plot(
  NA,
  xlim = c(0.5, 6.5),
  ylim = c(10, 95),
  xaxt = "n",
  xlab = "UDG (AGC)",
  ylab = expression("Circular speed at " * R[sys] * " (km s"^{-1}*")"),
  main = "Beardplot: Observed vs Newtonian, MOND, and HMG circular speeds"
)
axis(1, at = 1:6, labels = table2$AGC)
grid(nx = NA, ny = NULL, col = "gray90", lty = 1)

series <- list(
  list(name = "Observed", x = 1:6 - 0.27, y = table2$v_obs, plus = table2$v_obs_plus, minus = table2$v_obs_minus, col = "#F8766D"),
  list(name = "Newtonian", x = 1:6 - 0.09, y = table2$v_newton, plus = table2$v_newton_plus, minus = table2$v_newton_minus, col = "#7CAE00"),
  list(name = "MOND", x = 1:6 + 0.09, y = table2$v_mond, plus = table2$v_mond_plus, minus = table2$v_mond_minus, col = "#00BFC4"),
  list(name = "HMG", x = 1:6 + 0.27, y = table2$v_hmg, plus = table2$v_hmg_plus, minus = table2$v_hmg_minus, col = "#C77CFF")
)

for (s in series) {
  arrows(s$x, s$y - s$minus, s$x, s$y + s$plus, angle = 90, code = 3, length = 0.05, col = s$col, lwd = 2)
  points(s$x, s$y, pch = 16, cex = 1.1, col = s$col)
}
legend("top", legend = vapply(series, `[[`, character(1), "name"), col = vapply(series, `[[`, character(1), "col"), pch = 16, horiz = TRUE, bty = "n")
dev.off()

cat("\n# Exported files\n")
cat(file.path(out_dir, "table2_repro.csv"), "\n")
cat(file.path(out_dir, "table2_repro.tex"), "\n")
cat(file.path(out_dir, "vcirc_beardplot_repro.svg"), "\n")
