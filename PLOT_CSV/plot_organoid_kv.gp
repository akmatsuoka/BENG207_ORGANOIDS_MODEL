#!/usr/bin/gnuplot
# ============================================================
# plot_organoid_kv.gp
#
# Plots all 7 organoid viscoelastic deformation CSV outputs.
#
# Install:  sudo dnf install gnuplot
# Run:      gnuplot plot_organoid_kv.gp
# ============================================================

set datafile separator ','
set key autotitle columnheader
set grid
set style line 1 lw 2 lc rgb '#2B5797'
set style line 2 lw 2 lc rgb '#C0392B'
set style line 3 lw 2 lc rgb '#27AE60'
set style line 4 lw 2 lc rgb '#F39C12'
set style line 5 lw 1.5 lc rgb '#8E44AD' dt 2
set style line 6 lw 1.5 lc rgb '#16A085' dt 2

# ============================================================
# FIGURE 1: Single pulse — KV vs modified KV
# ============================================================

set terminal qt 0 title 'Fig 1: Single Pulse' size 800,500
set title 'Single Pulse Response — KV vs Modified KV' font ',14'
set xlabel 'Time (s)'
set ylabel 'Strain (%)'

plot 'organoid_fig1_single_pulse.csv' \
     using 1:($3*100) with lines ls 1 title 'Modified KV (total)', \
     '' using 1:($4*100) with lines ls 5 title 'Fast branch ({/Symbol t}_1)', \
     '' using 1:($5*100) with lines ls 6 title 'Slow branch ({/Symbol t}_2)', \
     '' using 1:($2*100) with lines ls 2 title 'Single KV'

pause -1 "Press Enter for next figure..."

# ============================================================
# FIGURE 2: Pulse train — residual accumulation
# ============================================================

set terminal qt 1 title 'Fig 2: Pulse Train' size 1000,500
set title '10-Pulse Train — Residual Strain Accumulation' font ',14'
set xlabel 'Time (s)'
set ylabel 'Strain (%)'

plot 'organoid_fig2_pulse_train.csv' \
     using 1:($3*100) with lines ls 1 title 'Modified KV', \
     '' using 1:($2*100) with lines ls 2 title 'Single KV'

pause -1 "Press Enter for next figure..."

# ============================================================
# FIGURE 3: Peak strain vs pulse duration
# ============================================================

set terminal qt 2 title 'Fig 3: Peak vs Duration' size 800,500
set title 'Peak Strain vs Pulse Duration' font ',14'
set xlabel 'Pulse Duration T_{on} (s)'
set ylabel 'Peak Strain (%)'
set logscale x

plot 'organoid_fig3_peak_vs_duration.csv' \
     using 1:($3*100) with lines ls 1 title 'Modified KV', \
     '' using 1:($2*100) with lines ls 2 title 'Single KV'

unset logscale x
pause -1 "Press Enter for next figure..."

# ============================================================
# FIGURE 4: Parameter sensitivity
# ============================================================

set terminal qt 3 title 'Fig 4: Sensitivity' size 1000,700
set multiplot layout 2,2 title 'Parameter Sensitivity' font ',14'

# Read the sensitivity CSV — it has columns:
# param_name, param_value, peak_strain, residual_strain

set xlabel 'E_1 (Pa)'; set ylabel 'Peak Strain (%)'
set title 'E_1 sensitivity'
plot '< awk -F, ''NR>1 && $1=="E1"'' organoid_fig4_sensitivity.csv' \
     using 2:($3*100) with linespoints ls 1 pt 7 ps 0.8 notitle

set xlabel 'E_2 (Pa)'
set title 'E_2 sensitivity'
plot '< awk -F, ''NR>1 && $1=="E2"'' organoid_fig4_sensitivity.csv' \
     using 2:($3*100) with linespoints ls 2 pt 7 ps 0.8 notitle

set xlabel '{/Symbol t}_1 (s)'
set title '{/Symbol t}_1 sensitivity'
plot '< awk -F, ''NR>1 && $1=="tau1"'' organoid_fig4_sensitivity.csv' \
     using 2:($3*100) with linespoints ls 3 pt 7 ps 0.8 notitle

set xlabel '{/Symbol t}_2 (s)'
set title '{/Symbol t}_2 sensitivity'
plot '< awk -F, ''NR>1 && $1=="tau2"'' organoid_fig4_sensitivity.csv' \
     using 2:($3*100) with linespoints ls 4 pt 7 ps 0.8 notitle

unset multiplot
pause -1 "Press Enter for next figure..."

# ============================================================
# FIGURE 5: Protocol comparison (short / long / burst)
# ============================================================

set terminal qt 4 title 'Fig 5: Protocols' size 1000,700
set multiplot layout 3,1 title 'Protocol Comparison' font ',14'

set ylabel 'Strain (%)'

set title 'Short pulse (100 ms on, 2 s off)'
set xlabel 'Time (s)'
plot 'organoid_fig5_protocols.csv' \
     using 1:($2*100) with lines ls 1 notitle

set title 'Long pulse (2 s on, 30 s off)'
plot 'organoid_fig5_protocols.csv' \
     using 1:($3*100) with lines ls 2 notitle

set title 'Burst (50 ms on, 200 ms off, x20)'
plot 'organoid_fig5_protocols.csv' \
     using 1:($4*100) with lines ls 3 notitle

unset multiplot
pause -1 "Press Enter for next figure..."

# ============================================================
# FIGURE 6: Diameter change prediction
# ============================================================

set terminal qt 5 title 'Fig 6: Diameter Change' size 800,500
set title 'Predicted Diameter Change' font ',14'
set xlabel 'Time (s)'
set ylabel '{/Symbol D}D ({/Symbol m}m)'

# Detection limit line
set arrow from graph 0,first 1 to graph 1,first 1 nohead ls 4 dt 3
set label 'Detection limit (~1 {/Symbol m}m)' at graph 0.6,first 1.3

plot 'organoid_fig6_diameter_change.csv' \
     using 1:($2*1e6) with lines ls 1 title 'Modified KV'

unset arrow; unset label
pause -1 "Press Enter for next figure..."

# ============================================================
# FIGURE 7: Stress sweep — peak strain & diameter vs pressure
# ============================================================

set terminal qt 6 title 'Fig 7: Stress Sweep' size 900,500
set multiplot layout 1,2 title 'Acoustic Pressure Requirements' font ',14'

set title 'Peak Strain vs Pressure'
set xlabel 'Acoustic Pressure (kPa)'
set ylabel 'Peak Strain (%)'
plot 'organoid_fig7_stress_sweep.csv' \
     using ($1/1e3):($2*100) with lines ls 1 notitle

set title 'Peak {/Symbol D}D vs Pressure'
set ylabel '{/Symbol D}D ({/Symbol m}m)'
set arrow from graph 0,first 1 to graph 1,first 1 nohead ls 4 dt 3
set label '1 {/Symbol m}m limit' at graph 0.5,first 1.5

plot 'organoid_fig7_stress_sweep.csv' \
     using ($1/1e3):($3*1e6) with lines ls 1 notitle

unset arrow; unset label
unset multiplot
pause -1 "Press Enter to close all..."
