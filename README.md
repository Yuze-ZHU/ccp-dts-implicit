# CCP Dual-Time Implicit Solver

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22074319.svg)](https://doi.org/10.5281/zenodo.22074319)

This repository contains the time-domain reference implementation used for the one-dimensional argon capacitively coupled plasma benchmark. It contains only the BDF2 dual-time-stepping solver. Harmonic-balance and modified-Jacobi components are not included.

## Model

The conserved variables are electron number density, ion number density, and electron energy density:

$$
\mathbf{Q}=(n_e,n_i,\varepsilon_e)^T.
$$

The drift-diffusion equations are coupled to Poisson's equation. The electron-impact model contains direct ionization and ground-state excitation energy loss:

$$
R_{\mathrm{ion}}=k_{\mathrm{ion}}N_gn_e, 
R_{\mathrm{exc}}=k_{\mathrm{exc}}N_gn_e.
$$

The source terms are

$$
S_{n_e}=R_{\mathrm{ion}}, 
S_{n_i}=R_{\mathrm{ion}}, 
S_{\varepsilon_e}=-H_{\mathrm{ion}}R_{\mathrm{ion}}
-H_{\mathrm{exc}}R_{\mathrm{exc}}.
$$

The physical-time derivative is discretized with BDF2. Backward Euler is used for the first physical step. At each physical step, the nonlinear BDF residual is converged in pseudo time with a block implicit iteration. Poisson's equation is solved at every pseudo-time iteration and uses the pseudo-time semi-implicit correction employed in the manuscript.

## Benchmark

The supplied case reproduces the Lymberopoulos and Economou argon CCP configuration:

- discharge gap: 2.54 cm
- pressure: 1 Torr
- gas temperature: 300 K
- RF frequency: 13.56 MHz
- RF voltage amplitude: 100 V
- mesh: 90 nonuniform cells
- wall model: Lymberopoulos boundary condition
- wall electron temperature: 0.5 eV
- secondary-electron coefficient: 0.01
- surface-loss velocity: 1.19e5 m/s
- physical resolution: T/dt = 200
- pseudo-CFL number: 10000
- inner iterations per physical step: 300
- simulated duration: 1500 RF periods

The mesh is represented as a 90 by 1 finite-volume grid. The two transverse patches are symmetry boundaries, so the calculation is equivalent to the one-dimensional benchmark while retaining the original solver data structure.

## Two-dimensional paper case

cases/paper2d_halfheight.case provides the half-height 90 by 45 two-dimensional configuration used for the wall-effect study. The left RF electrode, right grounded electrode, and bottom grounded boundary are Lymberopoulos walls with a wall electron temperature of 0.5 eV. Only the upper boundary is a symmetry plane.

The supplied settings use BDF2 with T/dt = 100, a pseudo-CFL number of 12000, and 100 inner iterations per physical step. Run it with

```bash
./bin/ccp-dts-implicit -c cases/paper2d_halfheight.case
```

This two-dimensional calculation is substantially more expensive than the 90 by 1 benchmark and has a separate output directory, results/paper2d_halfheight.

## Dependencies

A C++17 compiler, Eigen 3.3.9 or newer, and a Trilinos installation containing Sacado and TeuchosCore are required.

The default paths match the parent project:

```text
DEPS_DIR=../../ccpDeps/install
EIGEN_DIR=../../eigen-3.3.9
```

Alternative installations can be selected on the command line.

## Build

```bash
make -j
```

With nondefault dependency paths:

```bash
make -j DEPS_DIR=/path/to/trilinos EIGEN_DIR=/path/to/eigen
```

The executable is written to `bin/ccp-dts-implicit`.

## Quick test

```bash
./bin/ccp-dts-implicit -c cases/smoke.case
```

The smoke case uses the same mesh, chemistry, transport model, and boundary conditions as the benchmark, but advances only 0.01 RF period with five inner iterations. It verifies compilation and the complete solver path; it is not a converged physical result.

## Full reproduction

```bash
./bin/ccp-dts-implicit -c cases/lymberopoulos1993.case
```

The full calculation is computationally expensive because it performs 200 physical steps per RF period and 300 pseudo-time iterations per physical step. Output is written to `results/lymberopoulos1993`.

After completion, verify the cycle-start spatial averages with

```bash
python3 tools/verify_benchmark.py
```

The reference values at the beginning of RF cycle 1500 are

```text
NeAve = 3.17696652e15 m^-3
NiAve = 3.30382591e15 m^-3
EeAve = 2.93689536e-3 J m^-3
```

The verification script uses a relative tolerance of 0.5 percent by default. A different tolerance can be supplied with `--tolerance`.

## Output

- `ResidualHistory.dat`: residual after each selected physical step
- `InnerResidualHistory.dat`: optional pseudo-time residual zones
- `NeNiEeAveDT.dat`: cycle-start volume averages
- `UnsteadyFlowFieldDT.dat`: phase-resolved field in the final RF period
- `finalSolution.dat`: final dimensional field
- `restart.dat`: final nondimensional cell state
- `mesh_geom.dat`: finite-volume mesh geometry

All field and history files use Tecplot-compatible ASCII formatting.

## Stop control

During a calculation, write the value `1` to `stop.dat` in the selected output directory. The solver checks this file before each physical time step.

## Reproduction data

The data directory contains the compact numerical histories and fields used for the manuscript convergence and validation comparisons. See data/README.md for file descriptions and reference-data attribution.

## License

The solver is distributed under the BSD 3-Clause License. See LICENSE.

## Citation

If you use this code, please cite the associated *Computer Physics Communications* paper and the archived software release:

Y. Zhu, *CCP Dual-Time Implicit Solver*, version 1.0.0, Zenodo (2026). https://doi.org/10.5281/zenodo.22074320

The version-specific DOI above identifies the exact `v1.0.0` release used for the manuscript. The DOI for all versions of the software is https://doi.org/10.5281/zenodo.22074319.
