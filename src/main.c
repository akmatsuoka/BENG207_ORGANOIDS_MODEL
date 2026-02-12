/*
 * organoid_kv_model.c
 *
 * Organoid Viscoelastic Deformation Model Under Acoustic Loading
 *
 * Models:
 *   1) Single Kelvin-Voigt (KV):  E*eps + eta*eps_dot = sigma(t)
 *   2) Modified KV (2-timescale): two parallel KV branches
 *
 * Organoid: 200 µm diameter iPSC-derived
 * Acoustic: 20 MHz standing wave, same stack as previous models
 *
 * Outputs CSV files for your MATLAB plotting.
 *
 * Compile:  gcc -std=c99 -O2 -o organoid_kv organoid_kv_model.c -lm
 * Run:      ./organoid_kv
 *
 * Once you run this on your computer and plot them on your platform, we can discuss results and then design our next steps.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Output directory — change YOUR_USERNAME to your actual username */
#define OUTPUT_DIR "/home/shh/Projects/BENG207_3/PLOT_CSV/"
/* ============================================================
 * Physical parameters
 * ============================================================ */

/* Organoid geometry */
static const double R_org   = 100.0e-6;    /* radius = 100 µm */
static const double D_org   = 200.0e-6;    /* diameter = 200 µm */

/* Fluid (water) */
static const double rho_f   = 1000.0;      /* kg/m^3 */
static const double c_f     = 1483.0;      /* m/s */
static const double mu_f    = 1.0e-3;      /* Pa·s dynamic viscosity */

/* Acoustic */
static const double f0      = 20.0e6;      /* 20 MHz */

/* Acoustic pressure from previous model (achievable in channel) */
static const double p_ac    = 0.5e6;       /* 0.5 MPa — moderate */

/* ---- Baseline KV parameters for organoids ---- */
/* These are representative; we should fit from data and so on */

/* Single KV */
static const double E_single  = 200.0;     /* Pa (soft organoid) */
static const double eta_single = 50.0;     /* Pa·s */
/* tau_single = eta/E = 0.25 s */

/* Modified KV: fast branch */
static const double E1_mod   = 500.0;      /* Pa */
static const double eta1_mod = 50.0;       /* Pa·s */
/* tau1 = 0.1 s (fast: cortex/cell deformation) */

/* Modified KV: slow branch */
static const double E2_mod   = 30.0;       /* Pa */
static const double eta2_mod = 300.0;      /* Pa·s */
/* tau2 = 10 s (slow: intercellular rearrangement, ECM) */

/* ============================================================
 * Acoustic stress estimation
 *
 * For organoid several wavelengths across (D/lambda ~ 2.7 at 20 MHz),
 * we use effective stress from radiation pressure:
 *   sigma_0 ~ p_ac^2 / (rho_f * c_f^2) * geometric_factor
 *
 * Or from force balance:
 *   F_ac ~ pi * R^2 * sigma_0
 *   sigma_0 calibrated from experiment
 *
 * Here we compute an estimate and also treat sigma_0 as sweepable. Does it make sense to ya?
 * ============================================================ */

static double estimate_sigma0(double p_acoustic)
{
    /* Acoustic radiation stress scale */
    /* sigma ~ p^2 / (2 * rho * c^2) * acoustic_contrast_factor */
    /* For organoid (density ~1050, c~1550), contrast factor ~0.1-0.3 */
    double alpha_contrast = 0.15;  /* conservative */
    double E_ac = p_acoustic * p_acoustic / (2.0 * rho_f * c_f * c_f);
    return 2.0 * E_ac * alpha_contrast;  /* effective uniaxial stress */
}

/* ============================================================
 * Utility
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
 * KV MODEL RESPONSES (closed-form)
 * ============================================================ */

/* Single KV: strain during pulse (0 < t < T) */
static double kv_single_on(double t, double sigma0, double E, double eta)
{
    double tau = eta / E;
    return (sigma0 / E) * (1.0 - exp(-t / tau));
}

/* Single KV: strain after pulse (t >= T) */
static double kv_single_off(double t, double T, double sigma0,
                             double E, double eta)
{
    double tau = eta / E;
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
    double tau1 = eta1 / E1;
    double tau2 = eta2 / E2;
    double eps1_T = (sigma0/E1) * (1.0 - exp(-T/tau1));
    double eps2_T = (sigma0/E2) * (1.0 - exp(-T/tau2));
    return eps1_T * exp(-(t-T)/tau1)
         + eps2_T * exp(-(t-T)/tau2);
}

/* ============================================================
 * Pulse train: compute strain for a sequence of ON/OFF pulses
 * Uses superposition of shifted pulse responses
 * ============================================================ */

/* Single KV pulse train */
static double kv_single_train(double t, double sigma0, double E, double eta,
                               double T_on, double T_off, int N_pulses)
{
    double eps = 0.0;
    double T_period = T_on + T_off;

    for (int n = 0; n < N_pulses; n++) {
        double t_start = n * T_period;
        double t_end   = t_start + T_on;

        if (t < t_start) break;  /* future pulse, no contribution */

        double t_local = t - t_start;

        if (t < t_end) {
            /* We are during this pulse */
            eps += kv_single_on(t_local, sigma0, E, eta);
        } else {
            /* Pulse has ended, in recovery */
            eps += kv_single_off(t - t_start, T_on, sigma0, E, eta);
        }
    }
    return eps;
}

/* Modified KV pulse train */
static double kv_mod_train(double t, double sigma0,
                            double E1, double eta1,
                            double E2, double eta2,
                            double T_on, double T_off, int N_pulses)
{
    double eps = 0.0;
    double T_period = T_on + T_off;

    for (int n = 0; n < N_pulses; n++) {
        double t_start = n * T_period;
        double t_end   = t_start + T_on;

        if (t < t_start) break;

        double t_local = t - t_start;

        if (t < t_end) {
            eps += kv_mod_on(t_local, sigma0, E1, eta1, E2, eta2);
        } else {
            eps += kv_mod_off(t - t_start, T_on, sigma0,
                              E1, eta1, E2, eta2);
        }
    }
    return eps;
}

/* ============================================================
 * OUTPUT 1: Single pulse response — single KV vs modified KV
 * ============================================================ */

static void output_single_pulse(double sigma0)
{
    const char *fname = "org_fig1_single_pulse.csv";
    FILE *fp = fopen(fname, "w");

    double T_on = 0.5;     /* 500 ms pulse */
    double T_total = 5.0;  /* observe for 5 s */
    int Nt = 2000;
    double t_arr[2000];
    linspace(0, T_total, Nt, t_arr);

    double E_s = E_single, eta_s = eta_single;
    double tau_s = eta_s / E_s;
    double tau1 = eta1_mod / E1_mod;
    double tau2 = eta2_mod / E2_mod;

    fprintf(fp, "t_s,sigma_norm,eps_single_pct,eps_mod_pct,"
                "eps_mod_fast_pct,eps_mod_slow_pct\n");

    for (int i = 0; i < Nt; i++) {
        double t = t_arr[i];
        double sig = (t < T_on) ? 1.0 : 0.0;

        double eps_s, eps_m, eps_m1, eps_m2;

        if (t < T_on) {
            eps_s  = kv_single_on(t, sigma0, E_s, eta_s);
            eps_m1 = (sigma0/E1_mod) * (1.0 - exp(-t/tau1));
            eps_m2 = (sigma0/E2_mod) * (1.0 - exp(-t/tau2));
        } else {
            eps_s  = kv_single_off(t, T_on, sigma0, E_s, eta_s);
            double e1T = (sigma0/E1_mod)*(1.0-exp(-T_on/tau1));
            double e2T = (sigma0/E2_mod)*(1.0-exp(-T_on/tau2));
            eps_m1 = e1T * exp(-(t-T_on)/tau1);
            eps_m2 = e2T * exp(-(t-T_on)/tau2);
        }
        eps_m = eps_m1 + eps_m2;

        fprintf(fp, "%.6f,%.1f,%.6f,%.6f,%.6f,%.6f\n",
                t, sig, eps_s*100, eps_m*100, eps_m1*100, eps_m2*100);
    }

    fclose(fp);
    printf("  Written: %s (tau_single=%.3fs, tau1=%.3fs, tau2=%.1fs)\n",
           fname, tau_s, tau1, tau2);
}

/* ============================================================
 * OUTPUT 2: Pulse train — residual strain accumulation
 * ============================================================ */

static void output_pulse_train(double sigma0)
{
    const char *fname = "org_fig2_pulse_train.csv";
    FILE *fp = fopen(fname, "w");

    int N_pulses = 10;
    double T_on  = 0.5;    /* 500 ms ON */
    double T_off = 2.0;    /* 2 s OFF */
    double T_total = N_pulses * (T_on + T_off) + 5.0;  /* extra recovery */

    int Nt = 5000;
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

        /* Is ultrasound on? */
        double T_period = T_on + T_off;
        int pulse_idx = (int)(t / T_period);
        double t_in_period = t - pulse_idx * T_period;
        int us_on = (pulse_idx < N_pulses && t_in_period < T_on) ? 1 : 0;

        fprintf(fp, "%.6f,%.6f,%.6f,%d\n",
                t, eps_s*100, eps_m*100, us_on);
    }

    fclose(fp);
    printf("  Written: %s (%d pulses, T_on=%.1fs, T_off=%.1fs)\n",
           fname, N_pulses, T_on, T_off);
}

/* ============================================================
 * OUTPUT 3: Effect of pulse duration — peak strain vs T_on
 * ============================================================ */

static void output_peak_vs_duration(double sigma0)
{
    const char *fname = "org_fig3_peak_vs_duration.csv";
    FILE *fp = fopen(fname, "w");

    int Nd = 200;
    double T_on_arr[200];
    logspace(-2, 1, Nd, T_on_arr);  /* 10 ms to 10 s */

    double tau_s = eta_single / E_single;
    double tau1  = eta1_mod / E1_mod;
    double tau2  = eta2_mod / E2_mod;

    fprintf(fp, "T_on_s,T_on_over_tau1,peak_single_pct,peak_mod_pct,"
                "peak_mod_fast_pct,peak_mod_slow_pct\n");

    for (int i = 0; i < Nd; i++) {
        double T = T_on_arr[i];

        /* Peak strain at end of pulse */
        double eps_s = kv_single_on(T, sigma0, E_single, eta_single);

        double eps_m1 = (sigma0/E1_mod) * (1.0 - exp(-T/tau1));
        double eps_m2 = (sigma0/E2_mod) * (1.0 - exp(-T/tau2));
        double eps_m  = eps_m1 + eps_m2;

        fprintf(fp, "%.6f,%.4f,%.6f,%.6f,%.6f,%.6f\n",
                T, T/tau1,
                eps_s*100, eps_m*100, eps_m1*100, eps_m2*100);
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * OUTPUT 4: Parameter sensitivity — sweep E1, E2, tau1, tau2
 *           Show peak strain and residual strain
 * ============================================================ */

static void output_parameter_sensitivity(double sigma0)
{
    const char *fname = "org_fig4_sensitivity.csv";
    FILE *fp = fopen(fname, "w");

    double T_on = 0.5;
    double T_off = 2.0;
    int N_pulses = 5;
    /* Evaluate residual after last pulse OFF */
    double t_peak = (N_pulses-1)*(T_on+T_off) + T_on;  /* end of last pulse */
    double t_resid = t_peak + T_off;                     /* after recovery */

    /* Sweep E1 */
    int Np = 50;

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
        double eta1_v = E1_mod * tau1_arr[i];
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
        double eta2_v = E2_mod * tau2_arr[i];
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
    printf("  Written: %s (4 parameter sweeps x %d points)\n", fname, Np);
}

/* ============================================================
 * OUTPUT 5: Varying pulse protocols — identifiability test
 *           Short pulse vs long pulse responses
 * ============================================================ */

static void output_protocol_comparison(double sigma0)
{
    const char *fname = "org_fig5_protocols.csv";
    FILE *fp = fopen(fname, "w");

    /* Protocol A: short pulse (100 ms ON, 2 s OFF) */
    /* Protocol B: long pulse (2 s ON, 10 s OFF) */
    /* Protocol C: rapid bursts (50 ms ON, 200 ms OFF) , we can always change these thou*/

    typedef struct { double T_on; double T_off; int N; const char *name; } Proto;
    Proto protos[] = {
        {0.1,  2.0,  5, "short_100ms"},
        {2.0,  10.0, 3, "long_2s"},
        {0.05, 0.2,  20, "burst_50ms"}
    };
    int Np = 3;

    fprintf(fp, "t_s,protocol");
    fprintf(fp, ",eps_mod_pct,sigma_on\n");

    for (int ip = 0; ip < Np; ip++) {
        Proto *pr = &protos[ip];
        double T_total = pr->N * (pr->T_on + pr->T_off) + 5.0;
        int Nt = 3000;
        double dt = T_total / (Nt - 1);

        for (int i = 0; i < Nt; i++) {
            double t = i * dt;

            double eps_m = kv_mod_train(t, sigma0,
                                         E1_mod, eta1_mod,
                                         E2_mod, eta2_mod,
                                         pr->T_on, pr->T_off, pr->N);

            double T_period = pr->T_on + pr->T_off;
            int pidx = (int)(t / T_period);
            double t_in = t - pidx * T_period;
            int us_on = (pidx < pr->N && t_in < pr->T_on) ? 1 : 0;

            fprintf(fp, "%.6f,%s,%.6f,%d\n",
                    t, pr->name, eps_m*100, us_on);
        }
    }

    fclose(fp);
    printf("  Written: %s (3 protocols)\n", fname);
}

/* ============================================================
 * OUTPUT 6: Diameter change prediction
 *           Convert strain to observable diameter change
 * ============================================================ */

static void output_diameter_change(double sigma0)
{
    const char *fname = "org_fig6_diameter.csv";
    FILE *fp = fopen(fname, "w");

    double T_on = 0.5;
    double T_off = 5.0;
    int N_pulses = 5;
    double T_total = N_pulses * (T_on + T_off) + 10.0;
    int Nt = 4000;
    double dt = T_total / (Nt - 1);

    double D0 = D_org * 1e6;  /* µm */

    fprintf(fp, "t_s,D_compressed_um,D_extended_um,"
                "eps_axial_pct,sigma_on\n");

    for (int i = 0; i < Nt; i++) {
        double t = i * dt;

        double eps = kv_mod_train(t, sigma0,
                                   E1_mod, eta1_mod,
                                   E2_mod, eta2_mod,
                                   T_on, T_off, N_pulses);

        /* Compressed axis (along acoustic force) */
        double D_comp = D0 * (1.0 - eps);

        /* Extended axis (assuming incompressibility: a*b^2 ~ const) */
        /* If eps_axial = -eps, then eps_transverse ~ +eps/2 */
        double D_ext = D0 * (1.0 + eps / 2.0);

        double T_period = T_on + T_off;
        int pidx = (int)(t / T_period);
        double t_in = t - pidx * T_period;
        int us_on = (pidx < N_pulses && t_in < T_on) ? 1 : 0;

        fprintf(fp, "%.6f,%.4f,%.4f,%.6f,%d\n",
                t, D_comp, D_ext, eps*100, us_on);
    }

    fclose(fp);
    printf("  Written: %s (D0 = %.0f µm)\n", fname, D0);
}

/* ============================================================
 * OUTPUT 7: Stress sweep — peak deformation vs acoustic pressure
 * ============================================================ */

static void output_stress_sweep(void)
{
    const char *fname = "org_fig7_stress_sweep.csv";
    FILE *fp = fopen(fname, "w");

    int Np = 100;
    double p_ac_arr[100];
    logspace(4, 6.5, Np, p_ac_arr);  /* 10 kPa to ~3 MPa */

    double T_on = 0.5;

    fprintf(fp, "p_ac_kPa,sigma0_Pa,peak_single_pct,peak_mod_pct,"
                "D_change_um_single,D_change_um_mod\n");

    for (int i = 0; i < Np; i++) {
        double sig = estimate_sigma0(p_ac_arr[i]);

        double eps_s = kv_single_on(T_on, sig, E_single, eta_single);
        double eps_m = kv_mod_on(T_on, sig, E1_mod, eta1_mod,
                                  E2_mod, eta2_mod);

        double dD_s = D_org * eps_s * 1e6;  /* µm */
        double dD_m = D_org * eps_m * 1e6;

        fprintf(fp, "%.4f,%.6f,%.6f,%.6f,%.4f,%.4f\n",
                p_ac_arr[i]/1e3, sig,
                eps_s*100, eps_m*100,
                dD_s, dD_m);
    }

    fclose(fp);
    printf("  Written: %s\n", fname);
}

/* ============================================================
 * MAIN
 * ============================================================ */

int main(void)
{
    double sigma0 = estimate_sigma0(p_ac);
    double lambda_f = c_f / f0;
    double D_over_lambda = D_org / lambda_f;

    printf("=====================================================\n");
    printf(" Organoid Viscoelastic Deformation Model\n");
    printf("=====================================================\n");
    printf(" Organoid diameter  = %.0f µm\n", D_org*1e6);
    printf(" Frequency          = %.0f MHz\n", f0/1e6);
    printf(" lambda_water       = %.1f µm\n", lambda_f*1e6);
    printf(" D / lambda         = %.2f  (>> 1: not small-particle)\n",
           D_over_lambda);
    printf(" Acoustic pressure  = %.0f kPa\n", p_ac/1e3);
    printf(" Estimated sigma_0  = %.2f Pa\n", sigma0);
    printf("-----------------------------------------------------\n");
    printf(" Single KV:  E=%.0f Pa, eta=%.0f Pa·s, tau=%.3f s\n",
           E_single, eta_single, eta_single/E_single);
    printf(" Mod KV fast: E1=%.0f Pa, eta1=%.0f Pa·s, tau1=%.3f s\n",
           E1_mod, eta1_mod, eta1_mod/E1_mod);
    printf(" Mod KV slow: E2=%.0f Pa, eta2=%.0f Pa·s, tau2=%.1f s\n",
           E2_mod, eta2_mod, eta2_mod/E2_mod);
    printf(" Max single-pulse strain (0.5s): %.2f%%\n",
           kv_mod_on(0.5, sigma0, E1_mod, eta1_mod,
                     E2_mod, eta2_mod) * 100);
    printf(" Diameter change at peak:  %.2f µm\n",
           D_org * kv_mod_on(0.5, sigma0, E1_mod, eta1_mod,
                             E2_mod, eta2_mod) * 1e6);
    printf("=====================================================\n\n");

    printf("Generating output files...\n\n");

    output_single_pulse(sigma0);
    output_pulse_train(sigma0);
    output_peak_vs_duration(sigma0);
    output_parameter_sensitivity(sigma0);
    output_protocol_comparison(sigma0);
    output_diameter_change(sigma0);
    output_stress_sweep();

    printf("\nDone. 7 CSV files generated.\n");

    return 0;
}
