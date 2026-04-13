"""
beng207_analysis.py
===================
I really hate the notion that we will have to pay for this commercial software!! 
This is totally against my lab's policy!!!

BENG207 Organoid KV Model — Data Analysis Pipeline
 
Reads the 10 CSV files produced by main.c (organoid_kv_model.c) and generates:
 
  PART A (single-source, uniaxial):
    Fig 1  — Single pulse: single KV vs modified KV, fast/slow branch decomposition
    Fig 2  — Pulse train: ratcheting accumulation
    Fig 3  — Peak strain vs pulse duration (log scale)
    Fig 4  — Parameter sensitivity (E1, E2, tau1, tau2)
    Fig 5  — Protocol comparison (short / long / burst)
    Fig 6  — Diameter change (compressed + extended axes)
    Fig 7  — Stress sweep: strain and diameter change vs p_ac
 
  PART B (dual PZT Config B, 741 kHz, biaxial):
    Fig 8  — Stress estimator comparison (old alpha=0.15 vs Gorkov dual-PZT)
    Fig 9  — Biaxial deformation: single pulse (9a) + pulse train (9b)
    Fig 10 — Phase sweep: sigma0(phi), eps_peak(phi), diameter vs phi
 
Usage
-----
    python beng207_analysis.py [--csv_dir PATH] [--out_dir PATH] [--show]
 
    --csv_dir   Directory containing the 10 CSV files  [default: ./PLOT_CSV]
    --out_dir   Directory for output PNG figures        [default: ./PLOT_FIGURES]
    --show      Show interactive plots (blocks until closed)
 
Requires: numpy, scipy, matplotlib, pandas
Install:  pip install numpy scipy matplotlib pandas
"""
 
import argparse
import os
import sys
import warnings
 
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
 
warnings.filterwarnings("ignore", category=RuntimeWarning)
 
# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
matplotlib.rcParams.update({
    "font.family":       "DejaVu Sans",
    "font.size":         10,
    "axes.labelsize":    11,
    "axes.titlesize":    12,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "lines.linewidth":   1.8,
    "legend.frameon":    False,
    "figure.dpi":        120,
    "savefig.dpi":       150,
    "savefig.bbox":      "tight",
})
 
# Color palette (colorblind-friendly)
C = {
    "single":  "#2166ac",   # blue      — single KV
    "mod":     "#d73027",   # red       — modified KV (total)
    "fast":    "#f46d43",   # orange    — fast branch
    "slow":    "#74add1",   # light blue— slow branch
    "uni":     "#4dac26",   # green     — uniaxial (Part A)
    "bi":      "#7b2d8b",   # purple    — biaxial (Part B)
    "old":     "#999999",   # grey      — old alpha=0.15 estimate
    "shade":   "#ffffcc",   # pale yellow — US-on shading
    "Dx":      "#d73027",   # red       — compressed axis
    "Dy":      "#2166ac",   # blue      — extended axis
    "phase":   "#1b7837",   # dark green— phase sweep
}
 
# Physical constants (must match main.c)
D0_um    = 200.0     # organoid diameter, µm
tau1     = 0.10      # s — fast branch
tau2     = 10.0      # s — slow branch
E1       = 500.0     # Pa
E2       = 30.0      # Pa
f0_B_kHz = 741.0     # kHz
lambda_mm = 1483.0 / (f0_B_kHz * 1e3) * 1e3  # mm  ~2.001
 

 
# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
 
def shade_us_on(ax, t, us_on_col, alpha=0.12):
    """Shade the regions where ultrasound is ON."""
    on = us_on_col.values.astype(bool)
    i = 0
    while i < len(on):
        if on[i]:
            j = i
            while j < len(on) and on[j]:
                j += 1
            ax.axvspan(t.iloc[i], t.iloc[j - 1], color=C["shade"],
                       alpha=alpha, linewidth=0, zorder=0)
            i = j
        else:
            i += 1
 
 
def save_fig(fig, out_dir, name):
    path = os.path.join(out_dir, name)
    fig.savefig(path)
    print(f"  Saved: {path}")
 
 
def modified_kv_pulse(t, A1, t1, A2, t2, T_on):
    """Two-branch KV response to a rectangular pulse — for scipy fitting."""
    eps = np.zeros_like(t, dtype=float)
    on  = t < T_on
    off = ~on
    eps[on]  = A1 * (1 - np.exp(-t[on] / t1)) + A2 * (1 - np.exp(-t[on] / t2))
    e1T = A1 * (1 - np.exp(-T_on / t1))
    e2T = A2 * (1 - np.exp(-T_on / t2))
    eps[off] = e1T * np.exp(-(t[off] - T_on) / t1) + \
               e2T * np.exp(-(t[off] - T_on) / t2)
    return eps
 
 
# ---------------------------------------------------------------------------
# PART A — Fig 1: Single pulse response
# ---------------------------------------------------------------------------
 
def plot_fig1(csv_dir, out_dir):
    df = pd.read_csv(os.path.join(csv_dir, "org_fig1_single_pulse.csv"))
    t  = df["t_s"]
 
    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    fig.suptitle("Fig 1 — Single Pulse Response (Part A, σ₀ = 17.1 Pa, T_on = 500 ms)",
                 fontweight="bold")
 
    # Top: total strain comparison
    ax = axes[0]
    shade_us_on(ax, t, df["sigma_norm"])
    ax.plot(t, df["eps_single_pct"],  color=C["single"], label="Single KV",    ls="--")
    ax.plot(t, df["eps_mod_pct"],     color=C["mod"],    label="Modified KV (total)")
    ax.axvline(0.5, color="k", lw=0.8, ls=":", alpha=0.5)
    ax.set_ylabel("Strain (%)")
    ax.set_title("Single vs. Modified KV")
    ax.legend()
    ax.annotate("US OFF", xy=(0.5, ax.get_ylim()[1] * 0.9),
                xytext=(0.7, ax.get_ylim()[1] * 0.9),
                fontsize=8, color="grey",
                arrowprops=dict(arrowstyle="->", color="grey", lw=0.8))
 
    # Bottom: branch decomposition
    ax = axes[1]
    shade_us_on(ax, t, df["sigma_norm"])
    ax.plot(t, df["eps_mod_pct"],      color=C["mod"],    label="Total (modified KV)")
    ax.plot(t, df["eps_mod_fast_pct"], color=C["fast"],   label=f"Fast branch (τ₁={tau1} s)",  ls="--")
    ax.plot(t, df["eps_mod_slow_pct"], color=C["slow"],   label=f"Slow branch (τ₂={tau2} s)",  ls="-.")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Strain (%)")
    ax.set_title("Branch Decomposition")
    ax.legend()
 
    # Annotate fast/slow recovery
    ax.annotate(f"Fast recovery\n~{5*tau1:.1f} s",
                xy=(0.5 + 5*tau1, 0.05), xytext=(1.5, 1.5),
                fontsize=8, color=C["fast"],
                arrowprops=dict(arrowstyle="->", color=C["fast"], lw=0.8))
    ax.annotate(f"Slow branch\nretains strain\nfor ~{5*tau2:.0f} s",
                xy=(4.5, df["eps_mod_slow_pct"].iloc[-1]),
                xytext=(2.5, 3.0), fontsize=8, color=C["slow"],
                arrowprops=dict(arrowstyle="->", color=C["slow"], lw=0.8))
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig1_single_pulse.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART A — Fig 2: Pulse train ratcheting
# ---------------------------------------------------------------------------
 
def plot_fig2(csv_dir, out_dir):
    df = pd.read_csv(os.path.join(csv_dir, "org_fig2_pulse_train.csv"))
    t  = df["t_s"]
 
    fig, ax = plt.subplots(figsize=(10, 4.5))
    fig.suptitle("Fig 2 — Pulse Train Accumulation (10 pulses, 500 ms ON / 2 s OFF)",
                 fontweight="bold")
 
    shade_us_on(ax, t, df["sigma_on"])
    ax.plot(t, df["eps_single_pct"], color=C["single"], label="Single KV",    ls="--", lw=1.2)
    ax.plot(t, df["eps_mod_pct"],    color=C["mod"],    label="Modified KV")
 
    # Mark peak of each pulse for modified KV
    peaks, _ = find_peaks(df["eps_mod_pct"].values, distance=50)
    ax.scatter(t.iloc[peaks], df["eps_mod_pct"].iloc[peaks],
               color=C["mod"], s=30, zorder=5, label="Modified KV peaks")
 
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Strain (%)")
    ax.legend()
 
    # Annotate ratcheting
    if len(peaks) >= 2:
        y1 = df["eps_mod_pct"].iloc[peaks[0]]
        y2 = df["eps_mod_pct"].iloc[peaks[-1]]
        ax.annotate("", xy=(t.iloc[peaks[-1]], y2), xytext=(t.iloc[peaks[0]], y1),
                    arrowprops=dict(arrowstyle="->", color="k", lw=1))
        ax.text((t.iloc[peaks[0]] + t.iloc[peaks[-1]]) / 2,
                (y1 + y2) / 2 + 0.3, "Ratcheting\n(slow branch)", fontsize=8,
                ha="center", color="k")
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig2_pulse_train.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART A — Fig 3: Peak strain vs pulse duration
# ---------------------------------------------------------------------------
 
def plot_fig3(csv_dir, out_dir):
    df = pd.read_csv(os.path.join(csv_dir, "org_fig3_peak_vs_duration.csv"))
 
    fig, ax = plt.subplots(figsize=(7, 4.5))
    fig.suptitle("Fig 3 — Peak Strain vs Pulse Duration", fontweight="bold")
 
    ax.semilogx(df["T_on_s"], df["peak_single_pct"],
                color=C["single"], ls="--", label="Single KV")
    ax.semilogx(df["T_on_s"], df["peak_mod_pct"],
                color=C["mod"], label="Modified KV (total)")
    ax.semilogx(df["T_on_s"], df["peak_mod_fast_pct"],
                color=C["fast"], ls="--", label="Fast branch only", lw=1.2)
    ax.semilogx(df["T_on_s"], df["peak_mod_slow_pct"],
                color=C["slow"], ls="-.", label="Slow branch only", lw=1.2)
 
    # Mark tau1 and tau2
    for tau, label, col in [(tau1, "τ₁ = 0.1 s", C["fast"]),
                             (tau2, "τ₂ = 10 s",  C["slow"])]:
        ax.axvline(tau, color=col, lw=1, ls=":", alpha=0.7)
        ax.text(tau * 1.1, ax.get_ylim()[0] + 0.5, label,
                color=col, fontsize=8, rotation=90, va="bottom")
 
    ax.set_xlabel("Pulse duration T_on (s)")
    ax.set_ylabel("Peak strain (%)")
    ax.legend(loc="upper left")
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig3_peak_vs_duration.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART A — Fig 4: Parameter sensitivity
# ---------------------------------------------------------------------------
 
def plot_fig4(csv_dir, out_dir):
    df   = pd.read_csv(os.path.join(csv_dir, "org_fig4_sensitivity.csv"))
    params = ["E1_Pa", "E2_Pa", "tau1_s", "tau2_s"]
    labels = ["E₁ (Pa)", "E₂ (Pa)", "τ₁ (s)", "τ₂ (s)"]
    colors = [C["fast"], C["slow"], C["fast"], C["slow"]]
 
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle("Fig 4 — Parameter Sensitivity (Modified KV)", fontweight="bold")
    axes = axes.flatten()
 
    for idx, (param, label, col) in enumerate(zip(params, labels, colors)):
        sub = df[df["param_name"] == param].copy()
        ax  = axes[idx]
        ax.semilogx(sub["param_value"], sub["peak_strain_pct"],
                    color=col,    label="Peak strain", lw=2)
        ax.semilogx(sub["param_value"], sub["residual_strain_pct"],
                    color=col,    label="Residual strain", ls="--", lw=1.5)
        ax.set_xlabel(label)
        ax.set_ylabel("Strain (%)")
        ax.set_title(f"Sweep: {label}")
        ax.legend(fontsize=8)
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig4_sensitivity.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART A — Fig 5: Protocol comparison
# ---------------------------------------------------------------------------
 
def plot_fig5(csv_dir, out_dir):
    df    = pd.read_csv(os.path.join(csv_dir, "org_fig5_protocols.csv"))
    protos = df["protocol"].unique()
    cols   = {"short_100ms": C["fast"], "long_2s": C["slow"], "burst_50ms": C["mod"]}
    names  = {"short_100ms": "Short (100 ms ON)", "long_2s": "Long (2 s ON)",
               "burst_50ms": "Burst (50 ms ON)"}
 
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=False)
    fig.suptitle("Fig 5 — Pulse Protocol Comparison", fontweight="bold")
 
    for ax, proto in zip(axes, protos):
        sub = df[df["protocol"] == proto].copy()
        t   = sub["t_s"]
        shade_us_on(ax, t, sub["sigma_on"])
        ax.plot(t, sub["eps_mod_pct"], color=cols.get(proto, "k"),
                label=names.get(proto, proto))
        ax.set_ylabel("Strain (%)")
        ax.set_title(names.get(proto, proto))
        ax.set_xlabel("Time (s)")
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig5_protocols.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART A — Fig 6: Diameter change (uniaxial)
# ---------------------------------------------------------------------------
 
def plot_fig6(csv_dir, out_dir):
    df = pd.read_csv(os.path.join(csv_dir, "org_fig6_diameter.csv"))
    t  = df["t_s"]
 
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    fig.suptitle(f"Fig 6 — Diameter Change Prediction (Part A, Uniaxial, D₀ = {D0_um:.0f} µm)",
                 fontweight="bold")
 
    ax = axes[0]
    shade_us_on(ax, t, df["sigma_on"])
    ax.plot(t, df["D_compressed_um"], color=C["Dx"], label="D_x (compressed, acoustic axis)")
    ax.plot(t, df["D_extended_um"],   color=C["Dy"], label="D_y (extended, transverse)")
    ax.axhline(D0_um, color="k", lw=0.8, ls=":", alpha=0.5, label=f"D₀ = {D0_um:.0f} µm")
    ax.set_ylabel("Diameter (µm)")
    ax.set_title("Absolute diameter")
    ax.legend(fontsize=8)
 
    ax = axes[1]
    shade_us_on(ax, t, df["sigma_on"])
    delta_comp = df["D_compressed_um"] - D0_um
    delta_ext  = df["D_extended_um"]   - D0_um
    ax.plot(t, delta_comp, color=C["Dx"], label="ΔD_x (compression)")
    ax.plot(t, delta_ext,  color=C["Dy"], label="ΔD_y (extension)")
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ΔDiameter (µm)")
    ax.set_title("Change from baseline")
    ax.legend(fontsize=8)
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig6_diameter.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART A — Fig 7: Stress sweep
# ---------------------------------------------------------------------------
 
def plot_fig7(csv_dir, out_dir):
    df = pd.read_csv(os.path.join(csv_dir, "org_fig7_stress_sweep.csv"))
 
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fig.suptitle("Fig 7 — Stress Sweep: Peak Deformation vs Acoustic Pressure (Part A)",
                 fontweight="bold")
 
    ax = axes[0]
    ax.loglog(df["p_ac_kPa"], df["peak_single_pct"],
              color=C["single"], ls="--", label="Single KV")
    ax.loglog(df["p_ac_kPa"], df["peak_mod_pct"],
              color=C["mod"], label="Modified KV")
    ax.axvline(500, color="grey", lw=0.8, ls=":", label="500 kPa (default)")
    ax.set_xlabel("Acoustic pressure (kPa)")
    ax.set_ylabel("Peak strain (%)")
    ax.set_title("Peak strain")
    ax.legend(fontsize=8)
 
    ax = axes[1]
    ax.loglog(df["p_ac_kPa"], df["D_change_um_single"],
              color=C["single"], ls="--", label="Single KV")
    ax.loglog(df["p_ac_kPa"], df["D_change_um_mod"],
              color=C["mod"], label="Modified KV")
    ax.axhline(0.3, color="grey", lw=0.8, ls=":",
               label="~0.3 µm (diffraction limit)")
    ax.axhline(1.0, color="grey", lw=0.8, ls="--",
               label="~1 µm (brightfield resolvable)")
    ax.axvline(500, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("Acoustic pressure (kPa)")
    ax.set_ylabel("Diameter change (µm)")
    ax.set_title("Diameter change")
    ax.legend(fontsize=8)
 
    # Annotate quadratic slope
    ax = axes[0]
    x0, x1 = 50, 500
    y0 = df.loc[(df["p_ac_kPa"] - x0).abs().idxmin(), "peak_mod_pct"]
    y1 = df.loc[(df["p_ac_kPa"] - x1).abs().idxmin(), "peak_mod_pct"]
    ax.annotate("slope = 2\n(strain ∝ p²)",
                xy=((x0*x1)**0.5, (y0*y1)**0.5),
                xytext=(20, 0.01), fontsize=8, color=C["mod"],
                arrowprops=dict(arrowstyle="->", color=C["mod"], lw=0.8))
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig7_stress_sweep.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART B — Fig 8: Stress estimator comparison
# ---------------------------------------------------------------------------
 
def plot_fig8(csv_dir, out_dir):
    df = pd.read_csv(os.path.join(csv_dir, "org_fig8_stress_comparison.csv"))
 
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    fig.suptitle("Fig 8 — Stress Estimator: Old (α=0.15) vs Gor'kov Dual-PZT (Part B)",
                 fontweight="bold")
 
    ax = axes[0]
    ax.loglog(df["p_ac_kPa"], df["sigma0_old_Pa"],
              color=C["old"], ls="--", label="Old: α=0.15 (single)")
    ax.loglog(df["p_ac_kPa"], df["sigma0_dual_Pa"],
              color=C["bi"],           label="Gor'kov Φ=0.059 (dual PZT)")
    ax.axvline(500, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("p_ac per PZT (kPa)")
    ax.set_ylabel("σ₀ (Pa)")
    ax.set_title("Radiation stress")
    ax.legend(fontsize=8)
 
    ax = axes[1]
    ax.semilogx(df["p_ac_kPa"], df["ratio"],
                color=C["bi"], lw=2)
    ax.axhline(df["ratio"].mean(), color="grey", lw=0.8, ls="--",
               label=f"Mean ratio = {df['ratio'].mean():.3f}")
    ax.axvline(500, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("p_ac per PZT (kPa)")
    ax.set_ylabel("σ₀_dual / σ₀_old")
    ax.set_title("Stress ratio (constant ~1.573)")
    ax.legend(fontsize=8)
 
    ax = axes[2]
    ax.loglog(df["p_ac_kPa"], df["dD_old_um"],
              color=C["old"], ls="--", label="Old estimate")
    ax.loglog(df["p_ac_kPa"], df["dD_dual_um"],
              color=C["bi"],           label="Gor'kov dual-PZT")
    ax.axhline(0.3, color="grey", lw=0.8, ls=":", label="0.3 µm (diffraction limit)")
    ax.axhline(1.0, color="grey", lw=0.8, ls="--", label="1 µm (resolvable)")
    ax.axvline(500, color="grey", lw=0.8, ls=":")
    ax.set_xlabel("p_ac per PZT (kPa)")
    ax.set_ylabel("ΔD (µm)")
    ax.set_title("Diameter change")
    ax.legend(fontsize=8)
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig8_stress_comparison.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART B — Fig 9: Biaxial deformation (9a single pulse, 9b pulse train)
# ---------------------------------------------------------------------------
 
def plot_fig9(csv_dir, out_dir):
 
    # ---- 9a: Single pulse ----
    df = pd.read_csv(os.path.join(csv_dir, "org_fig9a_biaxial_single_pulse.csv"))
    t  = df["t_s"]
 
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Fig 9a — Biaxial Deformation: Single Pulse (Part B vs Part A)",
                 fontweight="bold")
 
    # Strain comparison
    ax = axes[0, 0]
    shade_us_on(ax, t, df["sigma_on"])
    ax.plot(t, df["eps_uniaxial_pct"], color=C["uni"], ls="--",
            label=f"Uniaxial (σ₀=17.1 Pa)")
    ax.plot(t, df["eps_biaxial_pct"],  color=C["bi"],
            label=f"Biaxial (σ₀=26.8 Pa)")
    ax.set_ylabel("Strain (%)")
    ax.set_title("Strain: Uniaxial vs Biaxial")
    ax.legend(fontsize=8)
 
    # D_x comparison
    ax = axes[0, 1]
    shade_us_on(ax, t, df["sigma_on"])
    ax.plot(t, df["Dx_uni_um"], color=C["uni"], ls="--", label="D_x uniaxial")
    ax.plot(t, df["Dx_bi_um"],  color=C["bi"],           label="D_x biaxial")
    ax.axhline(D0_um, color="k", lw=0.8, ls=":", alpha=0.5)
    ax.set_ylabel("Diameter (µm)")
    ax.set_title("D_x (compressed, acoustic axis)")
    ax.legend(fontsize=8)
 
    # D_y comparison — key biaxial signature
    ax = axes[1, 0]
    shade_us_on(ax, t, df["sigma_on"])
    ax.plot(t, df["Dy_uni_um"], color=C["uni"], ls="--", label="D_y uniaxial (ε/2)")
    ax.plot(t, df["Dy_bi_um"],  color=C["bi"],           label="D_y biaxial (ε) — 2× larger")
    ax.axhline(D0_um, color="k", lw=0.8, ls=":", alpha=0.5)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Diameter (µm)")
    ax.set_title("D_y (extended, free transverse)")
    ax.legend(fontsize=8)
 
    # ΔD summary
    ax = axes[1, 1]
    shade_us_on(ax, t, df["sigma_on"])
    ax.plot(t, df["Dx_bi_um"] - D0_um,  color=C["Dx"], label="ΔD_x biaxial (compression)")
    ax.plot(t, df["Dy_bi_um"] - D0_um,  color=C["Dy"], label="ΔD_y biaxial (extension)")
    ax.plot(t, df["Dx_uni_um"] - D0_um, color=C["uni"], ls="--", label="ΔD_x uniaxial", lw=1)
    ax.plot(t, df["Dy_uni_um"] - D0_um, color=C["uni"], ls="-.", label="ΔD_y uniaxial", lw=1)
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ΔDiameter (µm)")
    ax.set_title("ΔD: biaxial vs uniaxial")
    ax.legend(fontsize=7)
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig9a_biaxial_single_pulse.png")
    plt.close(fig)
 
    # ---- 9b: Pulse train ----
    df = pd.read_csv(os.path.join(csv_dir, "org_fig9b_biaxial_pulse_train.csv"))
    t  = df["t_s"]
 
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    fig.suptitle("Fig 9b — Biaxial Deformation: Pulse Train (10 pulses, Part B vs Part A)",
                 fontweight="bold")
 
    ax = axes[0]
    shade_us_on(ax, t, df["sigma_on"])
    ax.plot(t, df["eps_uniaxial_pct"], color=C["uni"], ls="--", label="Uniaxial (Part A)")
    ax.plot(t, df["eps_biaxial_pct"],  color=C["bi"],           label="Biaxial (Part B)")
    ax.set_ylabel("Strain (%)")
    ax.set_title("Ratcheting: Biaxial shows stronger accumulation (higher σ₀)")
    ax.legend()
 
    ax = axes[1]
    shade_us_on(ax, t, df["sigma_on"])
    ax.plot(t, df["Dx_bi_um"] - D0_um, color=C["Dx"], label="ΔD_x biaxial (compression)")
    ax.plot(t, df["Dy_bi_um"] - D0_um, color=C["Dy"], label="ΔD_y biaxial (extension)")
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("ΔDiameter (µm)")
    ax.set_title("Biaxial diameter changes over pulse train")
    ax.legend()
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig9b_biaxial_pulse_train.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# PART B — Fig 10: Phase sweep
# ---------------------------------------------------------------------------
 
def plot_fig10(csv_dir, out_dir):
    df = pd.read_csv(os.path.join(csv_dir, "org_fig10_phase_sweep.csv"))
 
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    fig.suptitle(f"Fig 10 — Phase Sweep: σ₀(φ) and ε_peak(φ)  "
                 f"[λ = {lambda_mm:.3f} mm, λ/4 = {lambda_mm/4:.3f} mm]",
                 fontweight="bold")
 
    # σ₀(φ)
    ax = axes[0, 0]
    ax.plot(df["phi_deg"], df["sigma0_phi_Pa"], color=C["phase"], lw=2)
    phi_fit = np.linspace(0, 180, 300)
    sig_max = df["sigma0_phi_Pa"].max()
    ax.plot(phi_fit, sig_max * np.cos(np.radians(phi_fit) / 2)**2,
            color="k", ls="--", lw=1, label="Model: σ₀_max·cos²(φ/2)")
    ax.set_xlabel("Phase offset φ (°)")
    ax.set_ylabel("σ₀ (Pa)")
    ax.set_title("Radiation stress vs phase offset")
    ax.legend(fontsize=8)
    for phi_mark, label in [(0, "on-node"), (90, "−50%"), (180, "null")]:
        ax.axvline(phi_mark, color="grey", lw=0.8, ls=":")
        ax.text(phi_mark + 2, sig_max * 0.05, label, fontsize=7,
                color="grey", rotation=90, va="bottom")
 
    # ε_peak(φ)
    ax = axes[0, 1]
    ax.plot(df["phi_deg"], df["eps_peak_pct"], color=C["phase"], lw=2)
    eps_max = df["eps_peak_pct"].max()
    ax.plot(phi_fit, eps_max * np.cos(np.radians(phi_fit) / 2)**2,
            color="k", ls="--", lw=1, label="Model: ε_max·cos²(φ/2)")
    ax.set_xlabel("Phase offset φ (°)")
    ax.set_ylabel("Peak strain (%)")
    ax.set_title("Peak strain vs phase offset")
    ax.legend(fontsize=8)
 
    # Node shift
    ax = axes[1, 0]
    ax.plot(df["phi_deg"], df["node_shift_mm"], color=C["phase"], lw=2)
    ax.axhline(lambda_mm / 4, color="grey", lw=0.8, ls="--",
               label=f"λ/4 = {lambda_mm/4:.3f} mm")
    ax.set_xlabel("Phase offset φ (°)")
    ax.set_ylabel("Node shift (mm)")
    ax.set_title("Node position shift")
    ax.legend(fontsize=8)
 
    # Biaxial diameter at each phase
    ax = axes[1, 1]
    ax.plot(df["phi_deg"], df["Dx_bi_um"] - D0_um,
            color=C["Dx"], lw=2, label="ΔD_x (compression)")
    ax.plot(df["phi_deg"], df["Dy_bi_um"] - D0_um,
            color=C["Dy"], lw=2, label="ΔD_y (extension)")
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.axhline(-0.3, color="grey", lw=0.8, ls="--", label="±0.3 µm (diffraction limit)")
    ax.axhline(0.3,  color="grey", lw=0.8, ls="--")
    ax.set_xlabel("Phase offset φ (°)")
    ax.set_ylabel("ΔDiameter (µm)")
    ax.set_title("Biaxial diameter change vs phase\n(experimental calibration target)")
    ax.legend(fontsize=8)
 
    plt.tight_layout()
    save_fig(fig, out_dir, "fig10_phase_sweep.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# KV model fitter — run against Fig 1 data as validation
# ---------------------------------------------------------------------------
 
def fit_kv_model(csv_dir, out_dir):
    """
    Fit the two-branch KV model to Fig 1 (single pulse) data.
    This is the same function that will be applied to real experimental data.
    Here it is validated against synthetic (model-generated) data — the fit
    should recover the original parameters exactly.
    """
    df   = pd.read_csv(os.path.join(csv_dir, "org_fig1_single_pulse.csv"))
    t    = df["t_s"].values
    eps  = df["eps_mod_pct"].values / 100.0  # convert % to fraction
    T_on = 0.5  # s — known pulse duration
 
    # Fit — A1 = sigma0/E1, A2 = sigma0/E2
    try:
        popt, pcov = curve_fit(
            lambda t, A1, t1, A2, t2: modified_kv_pulse(t, A1, t1, A2, t2, T_on),
            t, eps,
            p0=[0.05, 0.1, 0.10, 10.0],
            bounds=([0, 0.01, 0, 0.5], [1.0, 2.0, 1.0, 100.0]),
            maxfev=10000
        )
        A1_fit, tau1_fit, A2_fit, tau2_fit = popt
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError as e:
        print(f"  [WARNING] KV fit did not converge: {e}")
        return
 
    # Recover material properties assuming known sigma0
    sigma0 = 17.051  # Pa — Part A value from model
    E1_fit   = sigma0 / A1_fit
    E2_fit   = sigma0 / A2_fit
    eta1_fit = E1_fit  * tau1_fit
    eta2_fit = E2_fit  * tau2_fit
 
    print("\n  --- KV Fit Results (validated against synthetic data) ---")
    print(f"  A1  = {A1_fit:.5f}  (σ₀/E₁, expected {sigma0/E1:.5f})")
    print(f"  τ₁  = {tau1_fit:.4f} s  (expected {tau1:.4f} s)")
    print(f"  A2  = {A2_fit:.5f}  (σ₀/E₂, expected {sigma0/E2:.5f})")
    print(f"  τ₂  = {tau2_fit:.4f} s  (expected {tau2:.4f} s)")
    print(f"  → E₁ = {E1_fit:.1f} Pa   (expected {E1:.1f} Pa)")
    print(f"  → E₂ = {E2_fit:.1f} Pa   (expected {E2:.1f} Pa)")
    print(f"  → η₁ = {eta1_fit:.1f} Pa·s")
    print(f"  → η₂ = {eta2_fit:.1f} Pa·s")
    print(f"  Parameter errors (1σ): A1±{perr[0]:.6f}, τ1±{perr[1]:.4f},"
          f" A2±{perr[2]:.6f}, τ2±{perr[3]:.4f}")
 
    # Plot fit vs data
    eps_fit = modified_kv_pulse(t, *popt, T_on=T_on)
 
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    fig.suptitle("KV Model Fit Validation (Synthetic Data — Should Match Exactly)",
                 fontweight="bold")
 
    ax = axes[0]
    ax.plot(t, eps * 100,     color=C["mod"], lw=2,   label="Data (from C model)")
    ax.plot(t, eps_fit * 100, color="k",       ls="--", label="Scipy fit", lw=1.5)
    ax.set_ylabel("Strain (%)")
    ax.set_title("Fit vs data")
    ax.legend()
 
    # Residuals
    ax = axes[1]
    residual_pct = (eps - eps_fit) * 100
    ax.plot(t, residual_pct, color=C["mod"], lw=1)
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Residual (%)")
    ax.set_title(f"Residuals — RMS = {np.sqrt(np.mean(residual_pct**2)):.4f}%")
 
    # Annotate recovered parameters
    txt = (f"τ₁ = {tau1_fit:.4f} s  (true: {tau1:.4f} s)\n"
           f"τ₂ = {tau2_fit:.4f} s  (true: {tau2:.4f} s)\n"
           f"E₁ = {E1_fit:.1f} Pa  (true: {E1:.1f} Pa)\n"
           f"E₂ = {E2_fit:.1f} Pa  (true: {E2:.1f} Pa)")
    axes[0].text(0.98, 0.95, txt, transform=axes[0].transAxes,
                 fontsize=8, va="top", ha="right",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey", alpha=0.8))
 
    plt.tight_layout()
    save_fig(fig, out_dir, "kv_fit_validation.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# Phase sweep fitter — validated against Fig 10 data
# ---------------------------------------------------------------------------
 
def fit_phase_sweep(csv_dir, out_dir):
    """
    Fit σ₀(φ) = σ₀_max · cos²(φ/2) to Fig 10 data.
    Applied to real data, this recovers σ₀_max and detects
    amplitude imbalance (asymmetric curve → p_L ≠ p_R).
    """
    df      = pd.read_csv(os.path.join(csv_dir, "org_fig10_phase_sweep.csv"))
    phi_rad = df["phi_rad"].values
    eps     = df["eps_peak_pct"].values
 
    def phase_model(phi, eps_max, phi_offset):
        """General form including phase offset — detects organoid off-node."""
        return eps_max * np.cos((phi + phi_offset) / 2)**2
 
    try:
        popt, pcov = curve_fit(
            phase_model, phi_rad, eps,
            p0=[eps.max(), 0.0],
            bounds=([0, -np.pi/2], [100, np.pi/2])
        )
        eps_max_fit, phi_offset_fit = popt
        perr = np.sqrt(np.diag(pcov))
    except RuntimeError as e:
        print(f"  [WARNING] Phase sweep fit did not converge: {e}")
        return
 
    print("\n  --- Phase Sweep Fit Results ---")
    print(f"  ε_max    = {eps_max_fit:.4f}%  (expected: {eps.max():.4f}%)")
    print(f"  φ_offset = {np.degrees(phi_offset_fit):.3f}°"
          f"  (0° = organoid perfectly on-node)")
    if abs(phi_offset_fit) > np.radians(5):
        print(f"  [NOTE] Non-zero offset suggests organoid is not at node center.")
    print(f"  Parameter errors: ε_max±{perr[0]:.4f}%, φ_offset±{np.degrees(perr[1]):.3f}°")
 
    # Plot
    phi_plot  = np.linspace(0, np.pi, 300)
    eps_model = phase_model(phi_plot, *popt)
 
    fig, ax = plt.subplots(figsize=(7, 4.5))
    fig.suptitle("Phase Sweep Fit Validation", fontweight="bold")
 
    ax.plot(np.degrees(phi_rad), eps,
            "o", color=C["phase"], ms=4, label="Data (from C model)")
    ax.plot(np.degrees(phi_plot), eps_model,
            color="k", ls="--", lw=1.5,
            label=f"Fit: ε_max={eps_max_fit:.2f}%, φ_offset={np.degrees(phi_offset_fit):.1f}°")
    ax.set_xlabel("Phase offset φ (°)")
    ax.set_ylabel("Peak strain (%)")
    ax.set_title("ε_peak(φ) = ε_max · cos²((φ + φ_offset)/2)\n"
                 "φ_offset ≠ 0 → organoid not at node center")
    ax.legend(fontsize=8)
 
    # Mark half-maximum
    phi_half = 90  # degrees
    ax.axvline(phi_half, color="grey", lw=0.8, ls=":")
    ax.axhline(eps_max_fit / 2, color="grey", lw=0.8, ls=":")
    ax.text(phi_half + 2, eps_max_fit / 2 + 0.1,
            f"φ=90°: ε={eps_max_fit/2:.2f}%\nnode shifted {lambda_mm/8:.3f} mm",
            fontsize=7, color="grey")
 
    plt.tight_layout()
    save_fig(fig, out_dir, "phase_sweep_fit_validation.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# Summary dashboard — one-page overview of key results
# ---------------------------------------------------------------------------
 
def plot_summary_dashboard(csv_dir, out_dir):
    """Single-page summary of the most important results for lab meetings."""
    fig = plt.figure(figsize=(14, 10))
    fig.suptitle("BENG207 Organoid KV Model — Results Summary Dashboard",
                 fontsize=14, fontweight="bold")
    gs = GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)
 
    # 1. Single pulse strain (Part A)
    ax1 = fig.add_subplot(gs[0, 0])
    df1 = pd.read_csv(os.path.join(csv_dir, "org_fig1_single_pulse.csv"))
    shade_us_on(ax1, df1["t_s"], df1["sigma_norm"])
    ax1.plot(df1["t_s"], df1["eps_mod_pct"],     color=C["mod"],  label="Modified KV")
    ax1.plot(df1["t_s"], df1["eps_single_pct"],  color=C["single"], ls="--", label="Single KV")
    ax1.set_title("Single pulse (Part A)", fontsize=9)
    ax1.set_xlabel("Time (s)", fontsize=8)
    ax1.set_ylabel("Strain (%)", fontsize=8)
    ax1.legend(fontsize=6)
 
    # 2. Stress sweep
    ax2 = fig.add_subplot(gs[0, 1])
    df7 = pd.read_csv(os.path.join(csv_dir, "org_fig7_stress_sweep.csv"))
    ax2.loglog(df7["p_ac_kPa"], df7["D_change_um_mod"], color=C["mod"])
    ax2.axhline(0.3, color="grey", lw=0.8, ls=":", label="0.3 µm limit")
    ax2.axhline(1.0, color="grey", lw=0.8, ls="--", label="1 µm limit")
    ax2.set_title("ΔD vs pressure (Part A)", fontsize=9)
    ax2.set_xlabel("p_ac (kPa)", fontsize=8)
    ax2.set_ylabel("ΔD (µm)", fontsize=8)
    ax2.legend(fontsize=6)
 
    # 3. Pulse train ratcheting
    ax3 = fig.add_subplot(gs[0, 2])
    df2 = pd.read_csv(os.path.join(csv_dir, "org_fig2_pulse_train.csv"))
    shade_us_on(ax3, df2["t_s"], df2["sigma_on"])
    ax3.plot(df2["t_s"], df2["eps_mod_pct"],    color=C["mod"],    label="Modified KV")
    ax3.plot(df2["t_s"], df2["eps_single_pct"], color=C["single"], ls="--", label="Single KV")
    ax3.set_title("Pulse train ratcheting (Part A)", fontsize=9)
    ax3.set_xlabel("Time (s)", fontsize=8)
    ax3.set_ylabel("Strain (%)", fontsize=8)
    ax3.legend(fontsize=6)
 
    # 4. Stress comparison (Part B)
    ax4 = fig.add_subplot(gs[1, 0])
    df8 = pd.read_csv(os.path.join(csv_dir, "org_fig8_stress_comparison.csv"))
    ax4.loglog(df8["p_ac_kPa"], df8["sigma0_old_Pa"],  color=C["old"], ls="--", label="Old α=0.15")
    ax4.loglog(df8["p_ac_kPa"], df8["sigma0_dual_Pa"], color=C["bi"],           label="Gor'kov dual")
    ax4.set_title("Stress estimator (Part B)", fontsize=9)
    ax4.set_xlabel("p_ac per PZT (kPa)", fontsize=8)
    ax4.set_ylabel("σ₀ (Pa)", fontsize=8)
    ax4.legend(fontsize=6)
 
    # 5. Biaxial vs uniaxial (Part B)
    ax5 = fig.add_subplot(gs[1, 1])
    df9 = pd.read_csv(os.path.join(csv_dir, "org_fig9a_biaxial_single_pulse.csv"))
    shade_us_on(ax5, df9["t_s"], df9["sigma_on"])
    ax5.plot(df9["t_s"], df9["eps_uniaxial_pct"], color=C["uni"], ls="--", label="Uniaxial")
    ax5.plot(df9["t_s"], df9["eps_biaxial_pct"],  color=C["bi"],           label="Biaxial")
    ax5.set_title("Biaxial vs uniaxial (Part B)", fontsize=9)
    ax5.set_xlabel("Time (s)", fontsize=8)
    ax5.set_ylabel("Strain (%)", fontsize=8)
    ax5.legend(fontsize=6)
 
    # 6. Phase sweep
    ax6 = fig.add_subplot(gs[1, 2])
    df10 = pd.read_csv(os.path.join(csv_dir, "org_fig10_phase_sweep.csv"))
    phi_fit = np.linspace(0, 180, 300)
    eps_max = df10["eps_peak_pct"].max()
    ax6.plot(df10["phi_deg"], df10["eps_peak_pct"], color=C["phase"], lw=2, label="Model")
    ax6.plot(phi_fit, eps_max * np.cos(np.radians(phi_fit) / 2)**2,
             color="k", ls="--", lw=1, label="cos²(φ/2) fit")
    ax6.set_title("Phase sweep ε(φ) (Part B)", fontsize=9)
    ax6.set_xlabel("Phase offset φ (°)", fontsize=8)
    ax6.set_ylabel("Peak strain (%)", fontsize=8)
    ax6.legend(fontsize=6)
 
    # 7. Peak strain vs duration
    ax7 = fig.add_subplot(gs[2, 0])
    df3 = pd.read_csv(os.path.join(csv_dir, "org_fig3_peak_vs_duration.csv"))
    ax7.semilogx(df3["T_on_s"], df3["peak_mod_pct"], color=C["mod"])
    ax7.axvline(tau1, color=C["fast"], lw=0.8, ls=":", label=f"τ₁={tau1}s")
    ax7.axvline(tau2, color=C["slow"], lw=0.8, ls=":", label=f"τ₂={tau2}s")
    ax7.set_title("Peak strain vs T_on (Part A)", fontsize=9)
    ax7.set_xlabel("T_on (s)", fontsize=8)
    ax7.set_ylabel("Peak strain (%)", fontsize=8)
    ax7.legend(fontsize=6)
 
    # 8. Phase sweep diameters
    ax8 = fig.add_subplot(gs[2, 1])
    ax8.plot(df10["phi_deg"], df10["Dx_bi_um"] - D0_um,
             color=C["Dx"], lw=2, label="ΔD_x (compression)")
    ax8.plot(df10["phi_deg"], df10["Dy_bi_um"] - D0_um,
             color=C["Dy"], lw=2, label="ΔD_y (extension)")
    ax8.axhline(0,    color="k",    lw=0.8, ls=":")
    ax8.axhline(-0.3, color="grey", lw=0.8, ls="--")
    ax8.axhline(0.3,  color="grey", lw=0.8, ls="--", label="±0.3 µm limit")
    ax8.set_title("Biaxial ΔD vs phase (Part B)", fontsize=9)
    ax8.set_xlabel("Phase offset φ (°)", fontsize=8)
    ax8.set_ylabel("ΔD (µm)", fontsize=8)
    ax8.legend(fontsize=6)
 
    # 9. Key numbers table
    ax9 = fig.add_subplot(gs[2, 2])
    ax9.axis("off")
    table_data = [
        ["Parameter",          "Part A",       "Part B"],
        ["σ₀ at 500 kPa",      "17.1 Pa",      "26.8 Pa"],
        ["ε_peak (500 ms)",    "6.2%",         "9.7%"],
        ["ΔD_x (500 ms)",      "12.4 µm",      "19.4 µm"],
        ["ΔD_y (500 ms)",      "6.2 µm",       "19.4 µm"],
        ["τ₁ (fast branch)",   "0.1 s",        "0.1 s"],
        ["τ₂ (slow branch)",   "10.0 s",       "10.0 s"],
        ["Node shift (φ=π)",   "N/A",          "0.50 mm (λ/4)"],
        ["Geometry",           "Uniaxial",     "Biaxial (disk)"],
        ["Couplant",           "US gel",       "None (immersed)"],
    ]
    tbl = ax9.table(cellText=table_data[1:], colLabels=table_data[0],
                    loc="center", cellLoc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.5)
    tbl.scale(1.1, 1.4)
    # Header styling
    for j in range(3):
        tbl[(0, j)].set_facecolor("#2c7bb6")
        tbl[(0, j)].set_text_props(color="white", fontweight="bold")
    # Part B column highlight
    for i in range(1, len(table_data)):
        tbl[(i, 2)].set_facecolor("#f0e6f6")
    ax9.set_title("Key Model Predictions", fontsize=9, fontweight="bold")
 
    save_fig(fig, out_dir, "summary_dashboard.png")
    plt.close(fig)
 
 
# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
 
def main():
    parser = argparse.ArgumentParser(description="BENG207 Analysis Pipeline")
    parser.add_argument("--csv_dir", default="./PLOT_CSV",
                        help="Directory containing the 10 CSV files from main.c")
    parser.add_argument("--out_dir", default="./PLOT_FIGURES",
                        help="Output directory for PNG figures")
    parser.add_argument("--show", action="store_true",
                        help="Show interactive plots (blocks until closed)")
    args = parser.parse_args()
 
    # Verify CSV directory
    if not os.path.isdir(args.csv_dir):
        print(f"ERROR: CSV directory not found: {args.csv_dir}")
        print(f"  Run main.c first to generate the CSV files, then re-run this script.")
        sys.exit(1)
 
    # Create output directory
    os.makedirs(args.out_dir, exist_ok=True)
 
    if args.show:
        matplotlib.use("TkAgg")
    else:
        matplotlib.use("Agg")
 
    print(f"\nBENG207 Analysis Pipeline")
    print(f"  CSV input:  {os.path.abspath(args.csv_dir)}")
    print(f"  Fig output: {os.path.abspath(args.out_dir)}")
    print()
 
    steps = [
        ("Part A — Fig 1: Single pulse",              plot_fig1),
        ("Part A — Fig 2: Pulse train",               plot_fig2),
        ("Part A — Fig 3: Peak vs duration",          plot_fig3),
        ("Part A — Fig 4: Parameter sensitivity",     plot_fig4),
        ("Part A — Fig 5: Protocol comparison",       plot_fig5),
        ("Part A — Fig 6: Diameter change",           plot_fig6),
        ("Part A — Fig 7: Stress sweep",              plot_fig7),
        ("Part B — Fig 8: Stress comparison",         plot_fig8),
        ("Part B — Fig 9: Biaxial deformation",       plot_fig9),
        ("Part B — Fig 10: Phase sweep",              plot_fig10),
        ("Validation — KV model fit",                 fit_kv_model),
        ("Validation — Phase sweep fit",              fit_phase_sweep),
        ("Summary dashboard",                         plot_summary_dashboard),
    ]
 
    for label, fn in steps:
        print(f"[{label}]")
        try:
            fn(args.csv_dir, args.out_dir)
        except FileNotFoundError as e:
            print(f"  [SKIP] CSV not found: {e}")
        except Exception as e:
            print(f"  [ERROR] {e}")
        print()
 
    print(f"\nDone. {len(steps)} steps complete.")
    print(f"Figures written to: {os.path.abspath(args.out_dir)}")
    print()
    print("To run against our experimental data:")
    print("  1. Replace CSV files in PLOT_CSV/ with our measured data")
    print("     (keep the same column names as the model output)")
    print("  2. Re-run: python beng207_analysis.py --csv_dir PLOT_CSV --out_dir PLOT_FIGURES")
    print("  3. The KV fit (kv_fit_validation.png) and phase sweep fit")
    print("     (phase_sweep_fit_validation.png) will recover E1, E2, tau1, tau2")
    print("     and sigma0_max from our experimental creep-recovery curves.")
 
 
if __name__ == "__main__":
    main()