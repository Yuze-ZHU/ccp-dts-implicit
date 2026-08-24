# Reproduction data

This directory contains the compact numerical data used for the convergence, validation, and two-dimensional comparisons in the associated manuscript.

## Contents

- time_step_refinement: cycle-start spatial averages for T/dt = 25, 50, 100, 200, and 400.
- inner_iterations: cycle-start spatial averages for 25, 50, 100, and 200 inner iterations.
- pseudo_cfl: pseudo-time residual histories for pseudo-CFL = 100, 500, 1000, 5000, 10000, 15000, and 20000.
- phase_profiles: phase-resolved electron-density profiles from the present solver and digitized 1993 reference profiles.
- two_dimensional: final fields and cycle-start averages for the lower-symmetry and lower-wall configurations, plus the compact phase-resolved electron-density field.

The DAT files retain their Tecplot-compatible headers and dimensional units. The CSV files contain comma-separated reference coordinates and values.

## Reference data

Files under phase_profiles/reference_1993 were digitized from:

D. P. Lymberopoulos and D. J. Economou, Fluid simulations of glow discharges: Effect of metastable atoms in argon, Journal of Applied Physics 73 (1993) 3668-3679. https://doi.org/10.1063/1.352926

These digitized values are included solely to reproduce the manuscript comparison. Copyright in the original publication remains with its publisher.

## Exclusions

Restart states, cluster logs, Tecplot layout files, duplicate meshes, and intermediate iteration dumps are intentionally excluded because they are not required to reproduce the reported comparisons.
