import math
import os
import csv
from dataclasses import dataclass


# ======================
# Reproducible HMG / MOND / Newton comparison for gas-rich UDGs
# - no external dependencies
# - notation aligned with the final manuscript
# - one-sided 95% lower limit on s uses Delta chi^2 = 2.71
# ======================

G_SI = 6.67430e-11
MSUN_SI = 1.98847e30
C_SI = 2.99792458e8
YR_SI = 365.25 * 24 * 3600
GYR_SI = 1e9 * YR_SI
KPC_SI = 3.085677581e19
KMPS = 1e3

T0_GYR = 13.8
T0_SI = T0_GYR * GYR_SI

GAMMA_U_DEG = 60.0
GAMMA_CEN_DEG = 90.0
GAMMA_U_S2 = math.sin(math.pi * GAMMA_U_DEG / 180.0) ** 2
GAMMA_CEN_S2 = math.sin(math.pi * GAMMA_CEN_DEG / 180.0) ** 2

A0_MOND_SI = 1.2e-10

MPC_TO_M = 3.085677581e22
H0_KM_S_MPC = 70.0
H0_SI = (H0_KM_S_MPC * 1000.0) / MPC_TO_M
RHO_VAC_SI = 3.0 * H0_SI**2 / (8.0 * math.pi * G_SI)

DELTA_CHI2_95_ONE_SIDED = 2.71


@dataclass(frozen=True)
class Galaxy:
    agc: str
    log10_mbar: float
    r_sys_kpc: float
    v_obs_kms: float
    sigma_plus_kms: float
    sigma_minus_kms: float
    v_newton_tab_kms: float
    v_newton_plus_kms: float
    v_newton_minus_kms: float


DATA = [
    Galaxy("114905", 9.21, 8.02, 23.0, 4.0, 6.0, 29.0, 8.0, 6.0),
    Galaxy("122966", 9.21, 10.80, 37.0, 5.0, 6.0, 25.0, 4.0, 4.0),
    Galaxy("219533", 9.36, 9.78, 37.0, 6.0, 5.0, 32.0, 9.0, 7.0),
    Galaxy("248945", 9.05, 8.55, 27.0, 3.0, 3.0, 24.0, 6.0, 5.0),
    Galaxy("334315", 9.25, 8.49, 25.0, 5.0, 5.0, 30.0, 6.0, 5.0),
    Galaxy("749290", 9.17, 8.47, 26.0, 6.0, 6.0, 27.0, 6.0, 5.0),
]


def baryonic_mass_si(log10_mbar: float) -> float:
    return (10.0 ** log10_mbar) * MSUN_SI


def v_newton_kms(galaxy: Galaxy) -> float:
    return galaxy.v_newton_tab_kms


def v_newton_spherical_kms(galaxy: Galaxy) -> float:
    mbar_si = baryonic_mass_si(galaxy.log10_mbar)
    r_si = galaxy.r_sys_kpc * KPC_SI
    return math.sqrt(G_SI * mbar_si / r_si) / KMPS


def v_hubble_kms(galaxy: Galaxy) -> float:
    return H0_KM_S_MPC * galaxy.r_sys_kpc / 1000.0


def epsilon_s_sq(galaxy: Galaxy, s: float) -> float:
    mbar_si = baryonic_mass_si(galaxy.log10_mbar)
    r_nei_si = s * galaxy.r_sys_kpc * KPC_SI
    rho_nei_si = 3.0 * mbar_si / (4.0 * math.pi * r_nei_si**3)
    return rho_nei_si / RHO_VAC_SI + 1.0 / 6.0


def gamma_s(galaxy: Galaxy, s: float) -> float:
    v_n = v_newton_kms(galaxy)
    v_h = v_hubble_kms(galaxy)
    eps2 = epsilon_s_sq(galaxy, s)
    ratio = abs((2.0 * v_n**2 - eps2 * v_h**2) / (2.0 * v_n**2 + eps2 * v_h**2))
    sin2_gamma = GAMMA_U_S2 + (GAMMA_CEN_S2 - GAMMA_U_S2) * ratio
    sin2_gamma = min(max(sin2_gamma, 0.0), 1.0)
    return math.asin(math.sqrt(sin2_gamma))


def g_hmg_si(galaxy: Galaxy, s: float) -> float:
    v_n_si = v_newton_kms(galaxy) * KMPS
    r_si = galaxy.r_sys_kpc * KPC_SI
    a_n_si = v_n_si**2 / r_si
    gamma = gamma_s(galaxy, s)
    return math.sqrt(a_n_si**2 + a_n_si * 2.0 * C_SI * H0_SI * math.cos(gamma) / gamma)


def v_hmg_kms(galaxy: Galaxy, s: float) -> float:
    r_si = galaxy.r_sys_kpc * KPC_SI
    return math.sqrt(g_hmg_si(galaxy, s) * r_si) / KMPS


def v_mond_kms(galaxy: Galaxy) -> float:
    r_si = galaxy.r_sys_kpc * KPC_SI
    g_n_si = (v_newton_spherical_kms(galaxy) * KMPS) ** 2 / r_si
    x = g_n_si / A0_MOND_SI
    nu_hat = 1.0 / (1.0 - math.exp(-math.sqrt(x)))
    g_mond_si = nu_hat * g_n_si
    return math.sqrt(g_mond_si * r_si) / KMPS


def sigma_eff_kms(model_minus_obs_kms: float, galaxy: Galaxy) -> float:
    return galaxy.sigma_plus_kms if model_minus_obs_kms >= 0.0 else galaxy.sigma_minus_kms


def chi2_hmg(s: float) -> float:
    total = 0.0
    for galaxy in DATA:
        model = v_hmg_kms(galaxy, s)
        residual = model - galaxy.v_obs_kms
        sigma = sigma_eff_kms(residual, galaxy)
        total += (residual / sigma) ** 2
    return total


def chi2_fixed(model_name: str) -> float:
    total = 0.0
    for galaxy in DATA:
        if model_name == "newton":
            model = v_newton_kms(galaxy)
        elif model_name == "mond":
            model = v_mond_kms(galaxy)
        else:
            raise ValueError(model_name)
        residual = model - galaxy.v_obs_kms
        sigma = sigma_eff_kms(residual, galaxy)
        total += (residual / sigma) ** 2
    return total


def combined_sigma_z(v_model: float, sigma_model_plus: float, sigma_model_minus: float, galaxy: Galaxy) -> float:
    if v_model >= galaxy.v_obs_kms:
        sigma_comb = math.sqrt(galaxy.sigma_plus_kms**2 + sigma_model_minus**2)
    else:
        sigma_comb = math.sqrt(galaxy.sigma_minus_kms**2 + sigma_model_plus**2)
    return (v_model - galaxy.v_obs_kms) / sigma_comb


def golden_section_minimum(func, a: float, b: float, tol: float = 1e-6) -> tuple[float, float]:
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = func(c)
    fd = func(d)
    while abs(b - a) > tol:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = func(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = func(d)
    x = 0.5 * (a + b)
    return x, func(x)


def bisect_root(func, a: float, b: float, tol: float = 1e-8, max_iter: int = 200) -> float:
    fa = func(a)
    fb = func(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        return float("nan")
    for _ in range(max_iter):
        m = 0.5 * (a + b)
        fm = func(m)
        if abs(fm) < tol or abs(b - a) < tol:
            return m
        if fa * fm <= 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def lower_limit_s_95(s_hat: float, chi2_min: float) -> float:
    target = chi2_min + DELTA_CHI2_95_ONE_SIDED
    a = 1.0
    b = s_hat
    fa = chi2_hmg(a) - target
    fb = chi2_hmg(b) - target
    steps = 0
    while fa * fb > 0.0 and steps < 60:
        a = max(1e-6, a / 2.0)
        fa = chi2_hmg(a) - target
        steps += 1
        if a <= 1e-6:
            break
    return bisect_root(lambda x: chi2_hmg(x) - target, a, b)


def asymptotic_hmg_model_errors() -> dict[str, tuple[float, float]]:
    # Model-side asymmetric intervals used in the final manuscript table.
    return {
        "114905": (5.72, 4.03),
        "122966": (2.64, 2.17),
        "219533": (7.05, 4.79),
        "248945": (3.96, 2.88),
        "334315": (4.80, 3.48),
        "749290": (3.75, 2.96),
    }


def mond_model_errors() -> dict[str, tuple[float, float]]:
    return {
        "114905": (9.50, 8.15),
        "122966": (6.07, 5.68),
        "219533": (12.0, 10.2),
        "248945": (8.28, 7.30),
        "334315": (8.21, 7.04),
        "749290": (6.97, 6.44),
    }


def project_dir() -> str:
    return os.path.dirname(os.path.abspath(__file__))


def format_pm(value: float, plus: float, minus: float, decimals: int = 1) -> str:
    fmt = f"{{:.{decimals}f}}"
    return f"${fmt.format(value)}^{{+{fmt.format(plus)}}}_{{-{fmt.format(minus)}}}$"


def write_table2_csv(rows: list[dict]) -> str:
    path = os.path.join(project_dir(), "table2_repro.csv")
    fieldnames = [
        "AGC",
        "v_obs",
        "v_obs_plus",
        "v_obs_minus",
        "v_newton",
        "v_newton_plus",
        "v_newton_minus",
        "z_newton",
        "v_mond",
        "v_mond_plus",
        "v_mond_minus",
        "z_mond",
        "v_hmg",
        "v_hmg_plus",
        "v_hmg_minus",
        "z_hmg",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_table2_tex(rows: list[dict]) -> str:
    path = os.path.join(project_dir(), "table2_repro.tex")
    lines = [
        "% Reproducible Table 2 generated by fit_hmg_udg.py",
        "\\begin{tabular}{lccccccc}",
        "\\toprule",
        "Galaxy & $v_{\\rm obs}$ & $v_\\text{\\tiny N }$ & $z_\\text{\\tiny N }$ & $v_{\\rm MOND}$ & $z_{\\rm MOND}$ & $v_{\\rm HMG,\\infty}$ & $z_{\\rm HMG,\\infty}$ \\\\",
        " & (km s$^{-1}$) & (km s$^{-1}$) & ($\\sigma$) & (km s$^{-1}$) & ($\\sigma$) & (km s$^{-1}$) & ($\\sigma$) \\\\",
        "\\midrule",
    ]
    for row in rows:
        lines.append(
            f"{row['AGC']} & "
            f"{format_pm(row['v_obs'], row['v_obs_plus'], row['v_obs_minus'], 0)} & "
            f"{format_pm(row['v_newton'], row['v_newton_plus'], row['v_newton_minus'], 0)} & "
            f"{row['z_newton']:.2f} & "
            f"{format_pm(row['v_mond'], row['v_mond_plus'], row['v_mond_minus'], 1)} & "
            f"{row['z_mond']:.2f} & "
            f"{format_pm(row['v_hmg'], row['v_hmg_plus'], row['v_hmg_minus'], 1)} & "
            f"{row['z_hmg']:.2f} \\\\"
        )
    lines += ["\\bottomrule", "\\end{tabular}", ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    return path


def write_beardplot_svg(rows: list[dict]) -> str:
    path = os.path.join(project_dir(), "vcirc_beardplot_repro.svg")
    width, height = 1100, 620
    margin_left, margin_right, margin_top, margin_bottom = 90, 40, 70, 70
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom
    y_min, y_max = 10.0, 95.0
    labels = [r["AGC"] for r in rows]
    colors = {
        "Observed": "#F8766D",
        "Newtonian": "#7CAE00",
        "MOND": "#00BFC4",
        "HMG": "#C77CFF",
    }
    offsets = {
        "Observed": -0.27,
        "Newtonian": -0.09,
        "MOND": 0.09,
        "HMG": 0.27,
    }

    def x_pos(index: int, offset: float) -> float:
        step = plot_w / len(rows)
        return margin_left + (index + 0.5 + offset * 0.7) * step

    def y_pos(value: float) -> float:
        frac = (value - y_min) / (y_max - y_min)
        return margin_top + plot_h * (1.0 - frac)

    series_rows = []
    for row in rows:
        series_rows.extend([
            ("Observed", row["AGC"], row["v_obs"], row["v_obs_plus"], row["v_obs_minus"]),
            ("Newtonian", row["AGC"], row["v_newton"], row["v_newton_plus"], row["v_newton_minus"]),
            ("MOND", row["AGC"], row["v_mond"], row["v_mond_plus"], row["v_mond_minus"]),
            ("HMG", row["AGC"], row["v_hmg"], row["v_hmg_plus"], row["v_hmg_minus"]),
        ])

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;} .small{font-size:14px;} .axis{font-size:16px;} .title{font-size:22px;font-weight:bold;} .legend{font-size:14px;}</style>',
        f'<text x="{width/2}" y="34" text-anchor="middle" class="title">Beardplot: Observed vs Newtonian, MOND, and HMG circular speeds</text>',
    ]
    for y_tick in [25, 50, 75]:
        y = y_pos(y_tick)
        parts.append(f'<line x1="{margin_left}" y1="{y}" x2="{width-margin_right}" y2="{y}" stroke="#dddddd" stroke-width="1"/>')
        parts.append(f'<text x="{margin_left-12}" y="{y+5}" text-anchor="end" class="axis">{y_tick}</text>')
    for i, label in enumerate(labels):
        x = x_pos(i, 0.0)
        parts.append(f'<line x1="{x}" y1="{margin_top}" x2="{x}" y2="{height-margin_bottom}" stroke="#e6e6e6" stroke-width="1"/>')
        parts.append(f'<text x="{x}" y="{height-margin_bottom+28}" text-anchor="middle" class="axis">{label}</text>')
    parts.append(f'<line x1="{margin_left}" y1="{height-margin_bottom}" x2="{width-margin_right}" y2="{height-margin_bottom}" stroke="black" stroke-width="1.5"/>')
    parts.append(f'<line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{height-margin_bottom}" stroke="black" stroke-width="1.5"/>')
    parts.append(f'<text x="{width/2}" y="{height-18}" text-anchor="middle" class="axis">UDG (AGC)</text>')
    parts.append(f'<g transform="translate(24 {height/2}) rotate(-90)"><text text-anchor="middle" class="axis">Circular speed at Rsys (km s^-1)</text></g>')

    legend_x = width / 2 - 210
    legend_y = 58
    for i, (name, color) in enumerate(colors.items()):
        lx = legend_x + i * 110
        parts.append(f'<circle cx="{lx}" cy="{legend_y}" r="5" fill="{color}"/>')
        parts.append(f'<text x="{lx+12}" y="{legend_y+5}" class="legend">{name}</text>')

    for i, row in enumerate(rows):
        for name, value, plus, minus in [
            ("Observed", row["v_obs"], row["v_obs_plus"], row["v_obs_minus"]),
            ("Newtonian", row["v_newton"], row["v_newton_plus"], row["v_newton_minus"]),
            ("MOND", row["v_mond"], row["v_mond_plus"], row["v_mond_minus"]),
            ("HMG", row["v_hmg"], row["v_hmg_plus"], row["v_hmg_minus"]),
        ]:
            x = x_pos(i, offsets[name])
            y = y_pos(value)
            y_top = y_pos(value + plus)
            y_bot = y_pos(value - minus)
            color = colors[name]
            parts.append(f'<line x1="{x}" y1="{y_top}" x2="{x}" y2="{y_bot}" stroke="{color}" stroke-width="3"/>')
            parts.append(f'<line x1="{x-6}" y1="{y_top}" x2="{x+6}" y2="{y_top}" stroke="{color}" stroke-width="3"/>')
            parts.append(f'<line x1="{x-6}" y1="{y_bot}" x2="{x+6}" y2="{y_bot}" stroke="{color}" stroke-width="3"/>')
            parts.append(f'<circle cx="{x}" cy="{y}" r="5" fill="{color}"/>')
    parts.append("</svg>")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(parts))
    return path


def main() -> None:
    s_hat, chi2_min = golden_section_minimum(chi2_hmg, 1.0, 1e6)
    s95 = lower_limit_s_95(s_hat, chi2_min)

    chi2_newton = chi2_fixed("newton")
    chi2_mond = chi2_fixed("mond")

    print("# Global fit summary")
    print(f"s_hat_at_scan_edge = {s_hat:.6f}")
    print(f"chi2_hmg_asymptotic = {chi2_min:.6f}")
    print(f"s_95_one_sided_lower = {s95:.6f}")
    print(f"chi2_newton = {chi2_newton:.6f}")
    print(f"chi2_mond_mls = {chi2_mond:.6f}")
    print()

    print("# Per-galaxy velocities at the asymptotic branch")
    print("AGC\tv_obs\tv_newton\tv_mond\tv_hmg_inf")
    for galaxy in DATA:
        print(
            f"{galaxy.agc}\t{galaxy.v_obs_kms:.3f}\t{v_newton_kms(galaxy):.3f}\t"
            f"{v_mond_kms(galaxy):.3f}\t{v_hmg_kms(galaxy, 1e6):.3f}"
        )
    print()

    hmg_err = asymptotic_hmg_model_errors()
    mond_err = mond_model_errors()
    chi2_comb_newton = 0.0
    chi2_comb_hmg = 0.0
    chi2_comb_mond = 0.0
    print("# Combined tensions z = (model - obs) / sigma_comb")
    print("AGC\tz_newton\tz_mond\tz_hmg")
    for galaxy in DATA:
        z_newton = combined_sigma_z(
            galaxy.v_newton_tab_kms,
            galaxy.v_newton_plus_kms,
            galaxy.v_newton_minus_kms,
            galaxy,
        )
        z_mond = combined_sigma_z(
            v_mond_kms(galaxy),
            mond_err[galaxy.agc][0],
            mond_err[galaxy.agc][1],
            galaxy,
        )
        z_hmg = combined_sigma_z(
            v_hmg_kms(galaxy, 1e6),
            hmg_err[galaxy.agc][0],
            hmg_err[galaxy.agc][1],
            galaxy,
        )
        chi2_comb_newton += z_newton**2
        chi2_comb_mond += z_mond**2
        chi2_comb_hmg += z_hmg**2
        print(f"{galaxy.agc}\t{z_newton:.2f}\t{z_mond:.2f}\t{z_hmg:.2f}")
    print()
    print("# Combined squared tensions")
    print(f"chi2_comb_newton = {chi2_comb_newton:.6f}")
    print(f"chi2_comb_hmg = {chi2_comb_hmg:.6f}")
    print(f"chi2_comb_mond = {chi2_comb_mond:.6f}")

    table_rows = []
    for galaxy in DATA:
        row = {
            "AGC": galaxy.agc,
            "v_obs": galaxy.v_obs_kms,
            "v_obs_plus": galaxy.sigma_plus_kms,
            "v_obs_minus": galaxy.sigma_minus_kms,
            "v_newton": galaxy.v_newton_tab_kms,
            "v_newton_plus": galaxy.v_newton_plus_kms,
            "v_newton_minus": galaxy.v_newton_minus_kms,
            "z_newton": combined_sigma_z(
                galaxy.v_newton_tab_kms,
                galaxy.v_newton_plus_kms,
                galaxy.v_newton_minus_kms,
                galaxy,
            ),
            "v_mond": v_mond_kms(galaxy),
            "v_mond_plus": mond_err[galaxy.agc][0],
            "v_mond_minus": mond_err[galaxy.agc][1],
            "z_mond": combined_sigma_z(
                v_mond_kms(galaxy),
                mond_err[galaxy.agc][0],
                mond_err[galaxy.agc][1],
                galaxy,
            ),
            "v_hmg": v_hmg_kms(galaxy, 1e6),
            "v_hmg_plus": hmg_err[galaxy.agc][0],
            "v_hmg_minus": hmg_err[galaxy.agc][1],
            "z_hmg": combined_sigma_z(
                v_hmg_kms(galaxy, 1e6),
                hmg_err[galaxy.agc][0],
                hmg_err[galaxy.agc][1],
                galaxy,
            ),
        }
        table_rows.append(row)

    csv_path = write_table2_csv(table_rows)
    tex_path = write_table2_tex(table_rows)
    svg_path = write_beardplot_svg(table_rows)
    print()
    print("# Exported files")
    print(csv_path)
    print(tex_path)
    print(svg_path)


if __name__ == "__main__":
    main()
