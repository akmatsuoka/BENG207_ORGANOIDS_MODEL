/*
 * organoid_kv_model.c
 *
 * Organoid Viscoelastic Deformation Model Under Acoustic Loading
 *
 * Models:
 *   1) Single Kelvin-Voigt (KV):  E*eps + eta*eps_dot = sigma(t)
 *   2) Modified KV (2-timescale): two KV branches in series
 *
 * Organoid: 200 µm diameter iPSC-derived
 * Acoustic: standing wave at 741 kHz (dual PZT Config B)
 *
 * ---------------------------------------------------------------
 * PART A (Figs 1-7): Original single-source uniaxial model
 *                    UNCHANGED from previous version
 * PART B (Figs 8-10): Dual-PZT Config B additions
 *   Fig 8  — Stress comparison: old alpha=0.15 vs Gorkov dual-PZT
 *             as function of p_ac (10 kPa – 3 MPa)
 *   Fig 9  — Biaxial KV deformation: single pulse + pulse train
 *             Overlay with uniaxial (Part A) for comparison
 *             Includes x-compression and y/z-extension (disk geometry)
 *   Fig 10 — Phase sweep: 20 organoids, phi = 0 -> pi
 *             sigma0(phi) and eps_peak(phi) as node shifts by lambda/4
 * ---------------------------------------------------------------
 *
 * Compile:  gcc -std=c99 -O2 -o organoid_kv organoid_kv_model.c -lm
 * Run:      ./organoid_kv
 *
 * Outputs:  7 CSVs (Part A) + 3 CSVs (Part B) = 10 total
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Output directory — change to your local path */
#define OUTPUT_DIR "/home/shh/Projects/BENG207_3/PLOT_CSV/"

/* ============================================================
 * PART A PHYSICAL PARAMETERS  (original — do not modify!!!!!!!!!!!!!!!)
 * ============================================================ */

/* Organoid geometry */
static const double R_org      = 100.0e-6;   /* radius = 100 µm        */
static const double D_org      = 200.0e-6;   /* diameter = 200 µm      */

/* Fluid (water) */
static const double rho_f      = 1000.0;     /* kg/m^3                 */
static const double c_f        = 1483.0;     /* m/s                    */
static const double mu_f       = 1.0e-3;     /* Pa·s dynamic viscosity */

/* Acoustic — Part A uses 20 MHz for D/lambda reporting */
static const double f0         = 20.0e6;     /* 20 MHz                 */

/* Acoustic pressure from previous model */
static const double p_ac       = 0.5e6;      /* 0.5 MPa — moderate     */

/* ---- Baseline KV parameters (shared Part A and Part B) ---- */

/* Single KV */
static const double E_single   = 200.0;      /* Pa (soft organoid)     */
static const double eta_single = 50.0;       /* Pa·s                   */
/* tau_single = eta/E = 0.25 s */

/* Modified KV: fast branch (cortex / single-cell deformation) */
static const double E1_mod     = 500.0;      /* Pa                     */
static const double eta1_mod   = 50.0;       /* Pa·s, tau1 = 0.1 s    */

/* Modified KV: slow branch (ECM / intercellular rearrangement) */
static const double E2_mod     = 30.0;       /* Pa                     */
static const double eta2_mod   = 300.0;      /* Pa·s, tau2 = 10.0 s   */

/* ============================================================
 * PART B PHYSICAL PARAMETERS  (dual-PZT Config B, 741 kHz)
 * ============================================================ */

/* Operating frequency for Config B */
static const double f0_B       = 741.0e3;    /* 741 kHz                */

/* Gorkov acoustic contrast factor for organoid in water
 * f1 = 1 - (rho_f*c_f^2)/(rho_p*c_p^2)  monopole (compressibility)
 * f2 = 2*(rho_p - rho_f)/(2*rho_p + rho_f) dipole (density)
 * Phi = (1/3)*f1 - (1/2)*f2
 * Organoid: rho_p=1050 kg/m^3, c_p=1550 m/s
 * Pre-computed values from BENG207_2:
 */
static const double f1_gorkov  = 0.1282;     /* monopole factor        */
static const double f2_gorkov  = 0.0323;     /* dipole factor          */
static const double Phi_gorkov = 0.059;      /* acoustic contrast       */

/* Organoid acoustic properties */
static const double rho_p      = 1050.0;     /* kg/m^3                 */
static const double c_p        = 1550.0;     /* m/s                    */

/* Channel geometry */
static const double L_channel  = 20.0e-3;    /* 20 mm channel length   */
static const double N_organoids = 20;        /* ~20 nodes at 1mm spacing*/

/* Default pressure amplitude for each PZT (symmetric) */
static const double p_L        = 0.5e6;      /* left PZT, 500 kPa      */
static const double p_R        = 0.5e6;      /* right PZT, 500 kPa     */

/* ============================================================
 * PART A: Original stress estimator  (old alpha=0.15, unchanged)
 * Used ONLY by Part A figures to preserve exact original output
 * ============================================================ */

static double estimate_sigma0(double p_acoustic)
{
    double alpha_contrast = 0.15;  /* conservative hand-waved factor   */
    double E_ac = p_acoustic * p_acoustic / (2.0 * rho_f * c_f * c_f);
    return 2.0 * E_ac * alpha_contrast;
}

/* ============================================================
 * PART B: Gorkov-linked dual-PZT stress estimator
 *
 * For a perfect counter-propagating standing wave (dual PZT):
 *   p(x,t) = (p_L + p_R) * cos(kx) * cos(wt)
 *
 * Acoustic energy density at antinode:
 *   E_ac = p0^2 / (2 * rho_f * c_f^2)
 *   where p0 = p_L + p_R  (superposition of two traveling waves)
 *
 * Gorkov radiation stress (small-particle regime, D/lambda ~ 0.1):
 *   sigma0 = Phi * p0^2 / (rho_f * c_f^2)
 *
 * This replaces the old alpha=0.15 with the physically derived Phi=0.059.
 * Note: Phi_gorkov < alpha_contrast because Phi accounts for the actual
 * compressibility and density contrast of soft biological tissue, which
 * has much smaller contrast vs. water than polystyrene beads (Phi~0.22).
 *
 * p_left, p_right: individual PZT pressure amplitudes
 * Returns sigma0 in Pa
 * ============================================================ */

static double estimate_sigma0_dual(double p_left, double p_right)
{
    /* Superposition amplitude at antinode */
    double p0 = p_left + p_right;

    /* Acoustic energy density (time-averaged) */
    double E_ac = (p0 * p0) / (2.0 * rho_f * c_f * c_f);

    /* Gorkov radiation stress */
    return Phi_gorkov * p0 * p0 / (rho_f * c_f * c_f);

    (void)E_ac; /* suppress unused warning */
}

/* ============================================================
 * UTILITY
 * ============================================================ */

static void linspace(double a, double b, int n, double *out)
{
    for (int i = 0; i < n; i++)
        out[i] = a + (b - a) * i / (n - 1);
}

static void logspace(double a, double b, int n, double *out)
{
    for (int i = 0; i < n; i++)
        out[i] = pow(10.0, a + (b - a) * i / (n - 1));
}

/* ============================================================
 * KV MODEL RESPONSES (closed-form) — shared Part A and Part B
 * ============================================================ */

/* Single KV: strain during pulse (0 < t < T) */
static double kv_single_on(double t, double sigma0,
                            double E, double eta)
{
    double tau = eta / E;
    return (sigma0 / E) * (1.0 - exp(-t / tau));
}

/* Single KV: strain after pulse (t >= T) */
static double kv_single_off(double t, double T, double sigma0,
                             double E, double eta)
{
    double tau   = eta / E;
    double eps_T = (sigma0 / E) * (1.0 - exp(-T / tau));
    return eps_T * exp(-(t - T) / tau);
}

/* Modified KV (2-branch): strain during pulse */
static double kv_mod_on(double t, double sigma0,
                         double E1, double eta1,
                         double E2, double eta2)
{
    double tau1 = eta1 / E1;
    double tau2 = eta2 / E2;
    return (sigma0/E1) * (1.0 - exp(-t/tau1))
         + (sigma0/E2) * (1.0 - exp(-t/tau2));
}

/* Modified KV (2-branch): strain after pulse */
static double kv_mod_off(double t, double T, double sigma0,
                          double E1, double eta1,
                          double E2, double eta2)
{
    double tau1   = eta1 / E1;
    double tau2   = eta2 / E2;
    double eps1_T = (sigma0/E1) * (1.0 - exp(-T/tau1));
    double eps2_T = (sigma0/E2) * (1.0 - exp(-T/tau2));
    return eps1_T * exp(-(t-T)/tau1)
         + eps2_T * exp(-(t-T)/tau2);
}

/* ============================================================
 * PULSE TRAIN: superposition of shifted pulse responses
 * ============================================================ */

/* Single KV pulse train */
static double kv_single_train(double t, double sigma0,
                               double E, double eta,
                               double T_on, double T_off, int N_pulses)
{
    double eps      = 0.0;
    double T_period = T_on + T_off;

    for (int n = 0; n < N_pulses; n++) {
        double t_start = n * T_period;
        double t_end   = t_start + T_on;
        if (t < t_start) break;
        double t_local = t - t_start;
        if (t < t_end)
            eps += kv_single_on(t_local, sigma0, E, eta);
        else
            eps += kv_single_off(t - t_start, T_on, sigma0, E, eta);
    }
    return eps;
}

/* Modified KV pulse train */
static double kv_mod_train(double t, double sigma0,
                            double E1, double eta1,
                            double E2, double eta2,
                            double T_on, double T_off, int N_pulses)
{
    double eps      = 0.0;
    double T_period = T_on + T_off;

    for (int n = 0; n < N_pulses; n++) {
        double t_start = n * T_period;
        double t_end   = t_start + T_on;
        if (t < t_start) break;
        double t_local = t - t_start;
        if (t < t_end)
            eps += kv_mod_on(t_local, sigma0, E1, eta1, E2, eta2);
        else
            eps += kv_mod_off(t - t_start, T_on, sigma0,
                              E1, eta1, E2, eta2);
    }
    return eps;
}

/* ============================================================
 * ============================================================
 *  PART A — ORIGINAL OUTPUTS (Figs 1–7)  UNCHANGED
 * ============================================================
 * ============================================================ */

/* ============================================================
 * OUTPUT 1: Single pulse response — single KV vs modified KV
 * ============================================================ */

static void output_single_pulse(double sigma0)
{
    const char *fname = OUTPUT_DIR "org_fig1_single_pulse.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    double T_on    = 0.5;
    double T_total = 5.0;
    int    Nt      = 2000;
    double t_arr[2000];
    linspace(0, T_total, Nt, t_arr);

    double E_s   = E_single, eta_s = eta_single;
    double tau_s = eta_s / E_s;
    double tau1  = eta1_mod / E1_mod;
    double tau2  = eta2_mod / E2_mod;

    fprintf(fp, "t_s,sigma_norm,eps_single_pct,eps_mod_pct,"
                "eps_mod_fast_pct,eps_mod_slow_pct\n");

    for (int i = 0; i < Nt; i++) {
        double t   = t_arr[i];
        double sig = (t < T_on) ? 1.0 : 0.0;

        double eps_s, eps_m, eps_m1, eps_m2;

        if (t < T_on) {
            eps_s  = kv_single_on(t, sigma0, E_s, eta_s);
            eps_m1 = (sigma0/E1_mod) * (1.0 - exp(-t/tau1));
            eps_m2 = (sigma0/E2_mod) * (1.0 - exp(-t/tau2));
        } else {
            eps_s  = kv_single_off(t, T_on, sigma0, E_s, eta_s);
            double e1T = (sigma0/E1_mod) * (1.0 - exp(-T_on/tau1));
            double e2T = (sigma0/E2_mod) * (1.0 - exp(-T_on/tau2));
            eps_m1 = e1T * exp(-(t-T_on)/tau1);
            eps_m2 = e2T * exp(-(t-T_on)/tau2);
        }
        eps_m = eps_m1 + eps_m2;

        fprintf(fp, "%.6f,%.1f,%.6f,%.6f,%.6f,%.6f\n",
                t, sig, eps_s*100, eps_m*100, eps_m1*100, eps_m2*100);
    }

    fclose(fp);
    printf("  Written: org_fig1_single_pulse.csv"
           " (tau_single=%.3fs, tau1=%.3fs, tau2=%.1fs)\n",
           tau_s, tau1, tau2);
}

/* ============================================================
 * OUTPUT 2: Pulse train — residual strain accumulation
 * ============================================================ */

static void output_pulse_train(double sigma0)
{
    const char *fname = OUTPUT_DIR "org_fig2_pulse_train.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    int    N_pulses = 10;
    double T_on     = 0.5;
    double T_off    = 2.0;
    double T_total  = N_pulses * (T_on + T_off) + 5.0;

    int    Nt = 5000;
    double t_arr[5000];
    linspace(0, T_total, Nt, t_arr);

    fprintf(fp, "t_s,eps_single_pct,eps_mod_pct,sigma_on\n");

    for (int i = 0; i < Nt; i++) {
        double t = t_arr[i];

        double eps_s = kv_single_train(t, sigma0, E_single, eta_single,
                                        T_on, T_off, N_pulses);
        double eps_m = kv_mod_train(t, sigma0,
                                     E1_mod, eta1_mod,
                                     E2_mod, eta2_mod,
                                     T_on, T_off, N_pulses);

        double T_period    = T_on + T_off;
        int    pulse_idx   = (int)(t / T_period);
        double t_in_period = t - pulse_idx * T_period;
        int    us_on       = (pulse_idx < N_pulses && t_in_period < T_on)
                             ? 1 : 0;

        fprintf(fp, "%.6f,%.6f,%.6f,%d\n",
                t, eps_s*100, eps_m*100, us_on);
    }

    fclose(fp);
    printf("  Written: org_fig2_pulse_train.csv"
           " (%d pulses, T_on=%.1fs, T_off=%.1fs)\n",
           N_pulses, T_on, T_off);
}

/* ============================================================
 * OUTPUT 3: Effect of pulse duration — peak strain vs T_on
 * ============================================================ */

static void output_peak_vs_duration(double sigma0)
{
    const char *fname = OUTPUT_DIR "org_fig3_peak_vs_duration.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    int    Nd = 200;
    double T_on_arr[200];
    logspace(-2, 1, Nd, T_on_arr);  /* 10 ms to 10 s */

    double tau1 = eta1_mod / E1_mod;
    double tau2 = eta2_mod / E2_mod;

    fprintf(fp, "T_on_s,T_on_over_tau1,peak_single_pct,peak_mod_pct,"
                "peak_mod_fast_pct,peak_mod_slow_pct\n");

    for (int i = 0; i < Nd; i++) {
        double T = T_on_arr[i];

        double eps_s  = kv_single_on(T, sigma0, E_single, eta_single);
        double eps_m1 = (sigma0/E1_mod) * (1.0 - exp(-T/tau1));
        double eps_m2 = (sigma0/E2_mod) * (1.0 - exp(-T/tau2));
        double eps_m  = eps_m1 + eps_m2;

        fprintf(fp, "%.6f,%.4f,%.6f,%.6f,%.6f,%.6f\n",
                T, T/tau1,
                eps_s*100, eps_m*100, eps_m1*100, eps_m2*100);
    }

    fclose(fp);
    printf("  Written: org_fig3_peak_vs_duration.csv\n");
}

/* ============================================================
 * OUTPUT 4: Parameter sensitivity — sweep E1, E2, tau1, tau2
 * ============================================================ */

static void output_parameter_sensitivity(double sigma0)
{
    const char *fname = OUTPUT_DIR "org_fig4_sensitivity.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    double T_on    = 0.5;
    double T_off   = 2.0;
    int    N_pulses = 5;
    double t_peak  = (N_pulses-1)*(T_on+T_off) + T_on;
    double t_resid = t_peak + T_off;
    int    Np      = 50;

    fprintf(fp, "param_name,param_value,peak_strain_pct,residual_strain_pct\n");

    /* Sweep E1: 50 to 2000 Pa */
    double E1_arr[50];
    logspace(log10(50), log10(2000), Np, E1_arr);
    for (int i = 0; i < Np; i++) {
        double eps_peak = kv_mod_train(t_peak, sigma0,
                                        E1_arr[i], eta1_mod,
                                        E2_mod, eta2_mod,
                                        T_on, T_off, N_pulses);
        double eps_res  = kv_mod_train(t_resid, sigma0,
                                        E1_arr[i], eta1_mod,
                                        E2_mod, eta2_mod,
                                        T_on, T_off, N_pulses);
        fprintf(fp, "E1_Pa,%.2f,%.6f,%.6f\n",
                E1_arr[i], eps_peak*100, eps_res*100);
    }

    /* Sweep E2: 5 to 500 Pa */
    double E2_arr[50];
    logspace(log10(5), log10(500), Np, E2_arr);
    for (int i = 0; i < Np; i++) {
        double eps_peak = kv_mod_train(t_peak, sigma0,
                                        E1_mod, eta1_mod,
                                        E2_arr[i], eta2_mod,
                                        T_on, T_off, N_pulses);
        double eps_res  = kv_mod_train(t_resid, sigma0,
                                        E1_mod, eta1_mod,
                                        E2_arr[i], eta2_mod,
                                        T_on, T_off, N_pulses);
        fprintf(fp, "E2_Pa,%.2f,%.6f,%.6f\n",
                E2_arr[i], eps_peak*100, eps_res*100);
    }

    /* Sweep tau1: 0.01 to 2 s (vary eta1, keep E1 fixed) */
    double tau1_arr[50];
    logspace(-2, log10(2), Np, tau1_arr);
    for (int i = 0; i < Np; i++) {
        double eta1_v   = E1_mod * tau1_arr[i];
        double eps_peak = kv_mod_train(t_peak, sigma0,
                                        E1_mod, eta1_v,
                                        E2_mod, eta2_mod,
                                        T_on, T_off, N_pulses);
        double eps_res  = kv_mod_train(t_resid, sigma0,
                                        E1_mod, eta1_v,
                                        E2_mod, eta2_mod,
                                        T_on, T_off, N_pulses);
        fprintf(fp, "tau1_s,%.4f,%.6f,%.6f\n",
                tau1_arr[i], eps_peak*100, eps_res*100);
    }

    /* Sweep tau2: 1 to 100 s (vary eta2, keep E2 fixed) */
    double tau2_arr[50];
    logspace(0, 2, Np, tau2_arr);
    for (int i = 0; i < Np; i++) {
        double eta2_v   = E2_mod * tau2_arr[i];
        double eps_peak = kv_mod_train(t_peak, sigma0,
                                        E1_mod, eta1_mod,
                                        E2_mod, eta2_v,
                                        T_on, T_off, N_pulses);
        double eps_res  = kv_mod_train(t_resid, sigma0,
                                        E1_mod, eta1_mod,
                                        E2_mod, eta2_v,
                                        T_on, T_off, N_pulses);
        fprintf(fp, "tau2_s,%.4f,%.6f,%.6f\n",
                tau2_arr[i], eps_peak*100, eps_res*100);
    }

    fclose(fp);
    printf("  Written: org_fig4_sensitivity.csv"
           " (4 parameter sweeps x %d points)\n", Np);
}

/* ============================================================
 * OUTPUT 5: Varying pulse protocols — identifiability test
 * ============================================================ */

static void output_protocol_comparison(double sigma0)
{
    const char *fname = OUTPUT_DIR "org_fig5_protocols.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    typedef struct {
        double T_on; double T_off; int N; const char *name;
    } Proto;

    Proto protos[] = {
        {0.1,  2.0,  5,  "short_100ms"},
        {2.0,  10.0, 3,  "long_2s"},
        {0.05, 0.2,  20, "burst_50ms"}
    };
    int Np = 3;

    fprintf(fp, "t_s,protocol,eps_mod_pct,sigma_on\n");

    for (int ip = 0; ip < Np; ip++) {
        Proto *pr     = &protos[ip];
        double T_total = pr->N * (pr->T_on + pr->T_off) + 5.0;
        int    Nt     = 3000;
        double dt     = T_total / (Nt - 1);

        for (int i = 0; i < Nt; i++) {
            double t    = i * dt;
            double eps_m = kv_mod_train(t, sigma0,
                                         E1_mod, eta1_mod,
                                         E2_mod, eta2_mod,
                                         pr->T_on, pr->T_off, pr->N);

            double T_period = pr->T_on + pr->T_off;
            int    pidx     = (int)(t / T_period);
            double t_in     = t - pidx * T_period;
            int    us_on    = (pidx < pr->N && t_in < pr->T_on) ? 1 : 0;

            fprintf(fp, "%.6f,%s,%.6f,%d\n",
                    t, pr->name, eps_m*100, us_on);
        }
    }

    fclose(fp);
    printf("  Written: org_fig5_protocols.csv (3 protocols)\n");
}

/* ============================================================
 * OUTPUT 6: Diameter change prediction (uniaxial, Part A)
 * ============================================================ */

static void output_diameter_change(double sigma0)
{
    const char *fname = OUTPUT_DIR "org_fig6_diameter.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    double T_on    = 0.5;
    double T_off   = 5.0;
    int    N_pulses = 5;
    double T_total = N_pulses * (T_on + T_off) + 10.0;
    int    Nt      = 4000;
    double dt      = T_total / (Nt - 1);
    double D0      = D_org * 1e6;  /* µm */

    fprintf(fp, "t_s,D_compressed_um,D_extended_um,"
                "eps_axial_pct,sigma_on\n");

    for (int i = 0; i < Nt; i++) {
        double t   = i * dt;
        double eps = kv_mod_train(t, sigma0,
                                   E1_mod, eta1_mod,
                                   E2_mod, eta2_mod,
                                   T_on, T_off, N_pulses);

        /* Uniaxial: one free transverse axis extends by eps/2 */
        double D_comp = D0 * (1.0 - eps);
        double D_ext  = D0 * (1.0 + eps / 2.0);

        double T_period = T_on + T_off;
        int    pidx     = (int)(t / T_period);
        double t_in     = t - pidx * T_period;
        int    us_on    = (pidx < N_pulses && t_in < T_on) ? 1 : 0;

        fprintf(fp, "%.6f,%.4f,%.4f,%.6f,%d\n",
                t, D_comp, D_ext, eps*100, us_on);
    }

    fclose(fp);
    printf("  Written: org_fig6_diameter.csv (D0 = %.0f µm)\n", D0);
}

/* ============================================================
 * OUTPUT 7: Stress sweep — peak deformation vs acoustic pressure
 * ============================================================ */

static void output_stress_sweep(void)
{
    const char *fname = OUTPUT_DIR "org_fig7_stress_sweep.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    int    Np = 100;
    double p_ac_arr[100];
    logspace(4, 6.5, Np, p_ac_arr);  /* 10 kPa to ~3 MPa */

    double T_on = 0.5;

    fprintf(fp, "p_ac_kPa,sigma0_Pa,peak_single_pct,peak_mod_pct,"
                "D_change_um_single,D_change_um_mod\n");

    for (int i = 0; i < Np; i++) {
        double sig   = estimate_sigma0(p_ac_arr[i]);
        double eps_s = kv_single_on(T_on, sig, E_single, eta_single);
        double eps_m = kv_mod_on(T_on, sig,
                                  E1_mod, eta1_mod,
                                  E2_mod, eta2_mod);

        double dD_s = D_org * eps_s * 1e6;  /* µm */
        double dD_m = D_org * eps_m * 1e6;

        fprintf(fp, "%.4f,%.6f,%.6f,%.6f,%.4f,%.4f\n",
                p_ac_arr[i]/1e3, sig,
                eps_s*100, eps_m*100,
                dD_s, dD_m);
    }

    fclose(fp);
    printf("  Written: org_fig7_stress_sweep.csv\n");
}

/* ============================================================
 * ============================================================
 *  PART B — DUAL-PZT CONFIG B  (Figs 8–10)
 * ============================================================
 * ============================================================ */

/* ------------------------------------------------------------------
 * BIAXIAL KV HELPERS
 *
 * Physics: dual PZT squeezes organoid symmetrically from ±x.
 * The organoid is free to expand in y (height) and z (depth).
 * For an incompressible organoid (nu = 0.5):
 *
 *   eps_x  = -eps          (compression along acoustic axis)
 *   eps_y  = +eps          (extension, both free axes expand equally)
 *   eps_z  = +eps
 *
 * Volume conservation: (1-eps)(1+eps)^2 ~ 1  (first order in eps)
 *
 * The constitutive equation (KV) is unchanged — sigma0 still drives
 * the same KV creep. What changes is:
 *   1) sigma0 is computed via Gorkov (estimate_sigma0_dual)
 *   2) the diameter conversion uses biaxial geometry
 *
 * The strain eps returned by kv_mod_on/off is the AXIAL strain
 * (along x). The transverse strain is +eps (not +eps/2 as in uniaxial).
 * ------------------------------------------------------------------*/

/* Biaxial diameter conversion — returns D_x (compressed) and D_y (extended)
 * D_x = D0*(1 - eps)      (acoustic axis, squeezed from both sides)
 * D_y = D0*(1 + eps)      (transverse, free — disk expands outward)
 */
static void biaxial_diameters(double eps, double D0_um,
                               double *D_x_um, double *D_y_um)
{
    *D_x_um = D0_um * (1.0 - eps);
    *D_y_um = D0_um * (1.0 + eps);
}

/* ============================================================
 * OUTPUT 8: Stress comparison — old alpha=0.15 vs Gorkov dual-PZT
 *
 * Sweeps p_ac from 10 kPa to 3 MPa.
 * Columns:
 *   p_ac_kPa        — acoustic pressure amplitude (each PZT)
 *   sigma0_old_Pa   — Part A estimate (alpha=0.15, single source)
 *   sigma0_dual_Pa  — Part B estimate (Gorkov Phi=0.059, dual PZT)
 *   ratio           — sigma0_dual / sigma0_old
 *   peak_old_pct    — peak strain, old sigma0, modified KV, 500ms pulse
 *   peak_dual_pct   — peak strain, dual sigma0, modified KV, 500ms pulse
 *   dD_old_um       — diameter change (µm), old
 *   dD_dual_um      — diameter change (µm), dual
 *
 * Key physics note:
 *   Old alpha=0.15: sigma0 = 2 * 0.15 * p^2/(2*rho*c^2) = 0.15*p^2/(rho*c^2)
 *   New Gorkov dual: sigma0 = Phi * (2p)^2 / (rho*c^2) = 4*Phi*p^2/(rho*c^2)
 *   Ratio = 4*Phi / alpha = 4*0.059/0.15 = 1.573
 *   Dual PZT gives ~57% MORE stress than the old single-source estimate
 *   because: (a) two waves superpose at the node, (b) Phi from Gorkov
 *   theory vs hand-waved alpha nearly cancel, but the 4x superposition wins.
 * ============================================================ */

static void output_B_stress_comparison(void)
{
    const char *fname = OUTPUT_DIR "org_fig8_stress_comparison.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    int    Np = 120;
    double p_arr[120];
    logspace(4, 6.5, Np, p_arr);  /* 10 kPa to ~3 MPa per PZT */

    double T_on = 0.5;  /* 500 ms pulse */

    fprintf(fp, "p_ac_kPa,sigma0_old_Pa,sigma0_dual_Pa,ratio,"
                "peak_old_pct,peak_dual_pct,dD_old_um,dD_dual_um\n");

    for (int i = 0; i < Np; i++) {
        double p        = p_arr[i];
        double sig_old  = estimate_sigma0(p);          /* alpha=0.15, single */
        double sig_dual = estimate_sigma0_dual(p, p);  /* Gorkov, dual PZT   */
        double ratio    = sig_dual / sig_old;

        /* Peak strain at end of 500 ms pulse — modified KV */
        double eps_old  = kv_mod_on(T_on, sig_old,
                                     E1_mod, eta1_mod,
                                     E2_mod, eta2_mod);
        double eps_dual = kv_mod_on(T_on, sig_dual,
                                     E1_mod, eta1_mod,
                                     E2_mod, eta2_mod);

        double dD_old  = D_org * eps_old  * 1e6;  /* µm */
        double dD_dual = D_org * eps_dual * 1e6;

        fprintf(fp, "%.4f,%.6f,%.6f,%.4f,%.6f,%.6f,%.4f,%.4f\n",
                p/1e3,
                sig_old, sig_dual, ratio,
                eps_old*100, eps_dual*100,
                dD_old, dD_dual);
    }

    fclose(fp);

    /* Print key values at p_ac = 500 kPa */
    double sig_old_500  = estimate_sigma0(0.5e6);
    double sig_dual_500 = estimate_sigma0_dual(0.5e6, 0.5e6);
    printf("  Written: org_fig8_stress_comparison.csv\n");
    printf("    At 500 kPa per PZT:\n");
    printf("      sigma0 old  (alpha=0.15) = %.3f Pa\n", sig_old_500);
    printf("      sigma0 dual (Gorkov Phi) = %.3f Pa\n", sig_dual_500);
    printf("      Ratio = %.3f  (dual PZT gives ~%.0f%% more stress)\n",
           sig_dual_500/sig_old_500,
           (sig_dual_500/sig_old_500 - 1.0)*100);
}

/* ============================================================
 * OUTPUT 9: Biaxial KV deformation — single pulse + pulse train
 *
 * Uses sigma0_dual (Gorkov, dual PZT).
 * Outputs:
 *   (A) Single pulse (500 ms): uniaxial vs biaxial strain + diameters
 *   (B) Pulse train (10 pulses, 500 ms ON / 2 s OFF): ratcheting
 *
 * Biaxial geometry (dual PZT, organoid as incompressible disk):
 *   D_x = D0*(1 - eps)   [compressed, acoustic axis]
 *   D_y = D0*(1 + eps)   [extended, free transverse axis]
 *   D_z = D0*(1 + eps)   [extended, depth axis — same as D_y by symmetry]
 *
 * Compare with uniaxial (Part A, Fig 6):
 *   D_x = D0*(1 - eps)
 *   D_y = D0*(1 + eps/2) [only one free axis in uniaxial model]
 *
 * The acoustic axis compression is IDENTICAL between uniaxial and
 * biaxial for the same sigma0. The difference shows up in:
 *   1) sigma0 is larger in dual PZT (Fig 8)
 *   2) the free transverse expansion is larger (eps vs eps/2)
 *      This matters for organoid-organoid spacing and imaging.
 * ============================================================ */

static void output_B_biaxial_deformation(void)
{
    double sig_dual = estimate_sigma0_dual(p_L, p_R);
    double sig_old  = estimate_sigma0(p_ac);
    double D0       = D_org * 1e6;  /* µm */

    /* ---- Part 9A: Single pulse ---- */
    {
        const char *fname = OUTPUT_DIR "org_fig9a_biaxial_single_pulse.csv";
        FILE *fp = fopen(fname, "w");
        if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

        double T_on    = 0.5;
        double T_total = 8.0;
        int    Nt      = 3000;
        double dt      = T_total / (Nt - 1);

        fprintf(fp, "t_s,"
                    "eps_uniaxial_pct,eps_biaxial_pct,"
                    "Dx_uni_um,Dy_uni_um,"
                    "Dx_bi_um,Dy_bi_um,"
                    "sigma_on\n");

        for (int i = 0; i < Nt; i++) {
            double t = i * dt;

            /* Uniaxial strain (old sigma0, Part A) */
            double eps_uni = (t < T_on)
                ? kv_mod_on(t, sig_old, E1_mod, eta1_mod, E2_mod, eta2_mod)
                : kv_mod_off(t, T_on, sig_old,
                              E1_mod, eta1_mod, E2_mod, eta2_mod);

            /* Biaxial strain (dual sigma0, Part B) */
            double eps_bi = (t < T_on)
                ? kv_mod_on(t, sig_dual, E1_mod, eta1_mod, E2_mod, eta2_mod)
                : kv_mod_off(t, T_on, sig_dual,
                              E1_mod, eta1_mod, E2_mod, eta2_mod);

            /* Uniaxial diameters */
            double Dx_uni = D0 * (1.0 - eps_uni);
            double Dy_uni = D0 * (1.0 + eps_uni / 2.0);

            /* Biaxial diameters */
            double Dx_bi, Dy_bi;
            biaxial_diameters(eps_bi, D0, &Dx_bi, &Dy_bi);

            int us_on = (t < T_on) ? 1 : 0;

            fprintf(fp, "%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f,%d\n",
                    t,
                    eps_uni*100, eps_bi*100,
                    Dx_uni, Dy_uni,
                    Dx_bi, Dy_bi,
                    us_on);
        }

        fclose(fp);
        printf("  Written: org_fig9a_biaxial_single_pulse.csv\n");
    }

    /* ---- Part 9B: Pulse train ---- */
    {
        const char *fname = OUTPUT_DIR "org_fig9b_biaxial_pulse_train.csv";
        FILE *fp = fopen(fname, "w");
        if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

        int    N_pulses = 10;
        double T_on     = 0.5;
        double T_off    = 2.0;
        double T_total  = N_pulses * (T_on + T_off) + 5.0;
        int    Nt       = 6000;
        double dt       = T_total / (Nt - 1);

        fprintf(fp, "t_s,"
                    "eps_uniaxial_pct,eps_biaxial_pct,"
                    "Dx_uni_um,Dy_uni_um,"
                    "Dx_bi_um,Dy_bi_um,"
                    "sigma_on\n");

        for (int i = 0; i < Nt; i++) {
            double t = i * dt;

            double eps_uni = kv_mod_train(t, sig_old,
                                           E1_mod, eta1_mod,
                                           E2_mod, eta2_mod,
                                           T_on, T_off, N_pulses);
            double eps_bi  = kv_mod_train(t, sig_dual,
                                           E1_mod, eta1_mod,
                                           E2_mod, eta2_mod,
                                           T_on, T_off, N_pulses);

            double Dx_uni = D0 * (1.0 - eps_uni);
            double Dy_uni = D0 * (1.0 + eps_uni / 2.0);
            double Dx_bi, Dy_bi;
            biaxial_diameters(eps_bi, D0, &Dx_bi, &Dy_bi);

            double T_period    = T_on + T_off;
            int    pulse_idx   = (int)(t / T_period);
            double t_in_period = t - pulse_idx * T_period;
            int    us_on       = (pulse_idx < N_pulses &&
                                  t_in_period < T_on) ? 1 : 0;

            fprintf(fp, "%.6f,%.6f,%.6f,%.4f,%.4f,%.4f,%.4f,%d\n",
                    t,
                    eps_uni*100, eps_bi*100,
                    Dx_uni, Dy_uni,
                    Dx_bi, Dy_bi,
                    us_on);
        }

        fclose(fp);
        printf("  Written: org_fig9b_biaxial_pulse_train.csv\n");
    }

    /* Print summary */
    printf("    Biaxial vs uniaxial at 500 kPa, t=0.5s:\n");
    printf("      sigma0 uniaxial = %.3f Pa  ->  eps = %.4f%%\n",
           sig_old,
           kv_mod_on(0.5, sig_old, E1_mod, eta1_mod,
                     E2_mod, eta2_mod)*100);
    printf("      sigma0 biaxial  = %.3f Pa  ->  eps = %.4f%%\n",
           sig_dual,
           kv_mod_on(0.5, sig_dual, E1_mod, eta1_mod,
                     E2_mod, eta2_mod)*100);
}

/* ============================================================
 * OUTPUT 10: Phase sweep — 20 organoids, phi = 0 -> pi
 *
 * In dual-PZT Config B, the phase offset phi between the two
 * transducers shifts all pressure nodes simultaneously by:
 *   delta_x = phi / (2*k) = phi * lambda / (4*pi)
 *
 * At 741 kHz: lambda = c_f / f0_B = 2.0 mm
 *   phi = 0   -> nodes at x = n * lambda/2 = n * 1.0 mm
 *   phi = pi  -> nodes shift by lambda/4 = 0.5 mm
 *
 * Each organoid starts at a pressure node (x_n = n * lambda/2).
 * As phi increases, the node moves away from the organoid.
 * The organoid experiences a REDUCED stress because it is no
 * longer at the pressure antinode of the standing wave.
 *
 * Stress as function of phase offset for organoid n at x_n:
 *
 *   p(x, phi) = (p_L + p_R) * cos(k*x - phi/2) * cos(wt)
 *                  [standing wave with phase-shifted node]
 *
 *   At organoid location x_n = n * lambda/2  (original node):
 *   p_eff(phi) = (p_L + p_R) * |cos(phi/2)|
 *
 *   sigma0(phi) = Phi_gorkov * p_eff^2 / (rho_f * c_f^2)
 *               = sigma0_max * cos^2(phi/2)
 *
 * This gives maximum stress at phi=0 (organoid ON-node) and
 * zero stress at phi=pi (organoid at pressure node of antiphase
 * standing wave — equivalent to acoustic antinode of displacement).
 *
 * NOTE: All 20 organoids shift together because the phase offset
 * moves ALL nodes equally. The organoid index n only matters if
 * the organoid positions are fixed (e.g., stuck to substrate)
 * while the node moves — then organoid n feels a different phase
 * than organoid m.
 *
 * Columns:
 *   phi_deg         — phase offset 0 to 180 degrees
 *   phi_rad         — same in radians
 *   node_shift_mm   — node position shift (lambda/4 * phi/pi)
 *   p_eff_kPa       — effective pressure at original node positions
 *   sigma0_phi_Pa   — radiation stress at phase phi
 *   eps_peak_pct    — peak strain (500 ms pulse) at this phi
 *   Dx_bi_um        — compressed diameter (µm)
 *   Dy_bi_um        — extended diameter (µm)
 *   organoid_n      — organoid index (1 to N_organoids)
 *                     (all organoids identical for free-floating case)
 * ============================================================ */

static void output_B_phase_sweep(void)
{
    const char *fname = OUTPUT_DIR "org_fig10_phase_sweep.csv";
    FILE *fp = fopen(fname, "w");
    if (!fp) { fprintf(stderr, "ERROR: cannot open %s -- does OUTPUT_DIR exist?\n", fname); return; }

    double lambda  = c_f / f0_B;      /* wavelength at 741 kHz = ~2.0 mm */
    double k       = 2.0 * M_PI / lambda;
    double p0_each = p_L;              /* symmetric: p_L = p_R             */
    double p0_total = p_L + p_R;       /* superposition at antinode         */

    /* Maximum stress (phi=0, organoid perfectly on-node) */
    double sigma0_max = estimate_sigma0_dual(p_L, p_R);

    double T_on    = 0.5;  /* 500 ms reference pulse */
    double D0      = D_org * 1e6;  /* µm */

    int N_phi = 181;   /* 0 to 180 degrees in 1-degree steps */

    fprintf(fp, "phi_deg,phi_rad,node_shift_mm,p_eff_kPa,"
                "sigma0_phi_Pa,eps_peak_pct,Dx_bi_um,Dy_bi_um\n");

    for (int i = 0; i < N_phi; i++) {
        double phi_deg     = (double)i;
        double phi_rad     = phi_deg * M_PI / 180.0;

        /* Node shifts by phi/(2k) = phi*lambda/(4*pi) */
        double node_shift  = phi_rad * lambda / (4.0 * M_PI);  /* meters */

        /* Effective pressure at original node positions */
        double p_eff       = p0_total * fabs(cos(phi_rad / 2.0));

        /* Gorkov stress at this phase */
        double sigma0_phi  = Phi_gorkov * p_eff * p_eff
                             / (rho_f * c_f * c_f);

        /* Peak strain at end of 500 ms pulse */
        double eps_peak    = kv_mod_on(T_on, sigma0_phi,
                                        E1_mod, eta1_mod,
                                        E2_mod, eta2_mod);

        /* Biaxial diameter change */
        double Dx_bi, Dy_bi;
        biaxial_diameters(eps_peak, D0, &Dx_bi, &Dy_bi);

        fprintf(fp, "%.1f,%.6f,%.4f,%.4f,%.6f,%.6f,%.4f,%.4f\n",
                phi_deg, phi_rad,
                node_shift * 1e3,           /* mm */
                p_eff / 1e3,                /* kPa */
                sigma0_phi,
                eps_peak * 100,
                Dx_bi, Dy_bi);
    }

    fclose(fp);

    /* Print key landmarks */
    printf("  Written: org_fig10_phase_sweep.csv\n");
    printf("    lambda at 741 kHz = %.3f mm\n", lambda*1e3);
    printf("    lambda/2 spacing  = %.3f mm  (node-to-node)\n",
           lambda*1e3/2.0);
    printf("    phi=0  : sigma0 = %.3f Pa, eps = %.4f%%\n",
           sigma0_max,
           kv_mod_on(T_on, sigma0_max,
                     E1_mod, eta1_mod, E2_mod, eta2_mod)*100);
    printf("    phi=90 : sigma0 = %.3f Pa, node shifted %.3f mm\n",
           sigma0_max * 0.5,
           lambda*1e3/8.0);
    printf("    phi=180: sigma0 = 0.000 Pa (organoid at pressure node)\n");
    printf("    N_organoids = %.0f  (%.0f mm channel / 1 mm spacing)\n",
           N_organoids, L_channel*1e3);

    (void)k;
    (void)p0_each;
}

/* ============================================================
 * MAIN
 * ============================================================ */

int main(void)
{
    /* ---- Part A sigma0 (original, unchanged) ---- */
    double sigma0_A = estimate_sigma0(p_ac);

    /* ---- Part B sigma0 (Gorkov, dual PZT) ---- */
    double sigma0_B = estimate_sigma0_dual(p_L, p_R);

    double lambda_A      = c_f / f0;
    double lambda_B      = c_f / f0_B;
    double D_over_lam_A  = D_org / lambda_A;
    double D_over_lam_B  = D_org / lambda_B;

    printf("=====================================================\n");
    printf(" Organoid Viscoelastic Deformation Model\n");
    printf("=====================================================\n");
    printf(" Organoid diameter   = %.0f µm\n", D_org*1e6);
    printf("\n");
    printf(" PART A (original, single-source, 20 MHz)\n");
    printf("   f0               = %.0f MHz\n", f0/1e6);
    printf("   lambda_water     = %.1f µm\n", lambda_A*1e6);
    printf("   D / lambda       = %.2f  (not small-particle)\n",
           D_over_lam_A);
    printf("   p_ac             = %.0f kPa\n", p_ac/1e3);
    printf("   sigma0 (alpha)   = %.3f Pa  [alpha=0.15]\n", sigma0_A);
    printf("\n");
    printf(" PART B (dual PZT Config B, 741 kHz)\n");
    printf("   f0_B             = %.0f kHz\n", f0_B/1e3);
    printf("   lambda_water     = %.3f mm\n", lambda_B*1e3);
    printf("   D / lambda       = %.4f  (small-particle, Gorkov valid)\n",
           D_over_lam_B);
    printf("   p_L = p_R        = %.0f kPa each\n", p_L/1e3);
    printf("   p_total at node  = %.0f kPa  (superposition)\n",
           (p_L+p_R)/1e3);
    printf("   Phi_gorkov       = %.4f  (f1=%.4f, f2=%.4f)\n",
           Phi_gorkov, f1_gorkov, f2_gorkov);
    printf("   sigma0 (Gorkov)  = %.3f Pa\n", sigma0_B);
    printf("   Ratio B/A        = %.3f\n", sigma0_B/sigma0_A);
    printf("   N_organoids      = %.0f  (20mm / 1mm spacing)\n",
           N_organoids);
    printf("\n");
    printf(" KV parameters (shared)\n");
    printf("   Single KV:  E=%.0f Pa, eta=%.0f Pa·s, tau=%.3f s\n",
           E_single, eta_single, eta_single/E_single);
    printf("   Fast branch: E1=%.0f Pa, tau1=%.3f s\n",
           E1_mod, eta1_mod/E1_mod);
    printf("   Slow branch: E2=%.0f Pa, tau2=%.1f s\n",
           E2_mod, eta2_mod/E2_mod);
    printf("=====================================================\n\n");

    /* ---- Part A outputs (original, unchanged) ---- */
    printf("PART A — Generating original outputs (Figs 1-7)...\n\n");
    output_single_pulse(sigma0_A);
    output_pulse_train(sigma0_A);
    output_peak_vs_duration(sigma0_A);
    output_parameter_sensitivity(sigma0_A);
    output_protocol_comparison(sigma0_A);
    output_diameter_change(sigma0_A);
    output_stress_sweep();

    /* ---- Part B outputs (dual PZT) ---- */
    printf("\nPART B — Generating dual-PZT outputs (Figs 8-10)...\n\n");
    output_B_stress_comparison();
    output_B_biaxial_deformation();
    output_B_phase_sweep();

    printf("\nDone. 10 CSV files generated (7 Part A + 3 Part B).\n");
    printf("Part A: org_fig1 through org_fig7\n");
    printf("Part B: org_fig8_stress_comparison\n");
    printf("        org_fig9a_biaxial_single_pulse\n");
    printf("        org_fig9b_biaxial_pulse_train\n");
    printf("        org_fig10_phase_sweep\n");

    return 0;
}
