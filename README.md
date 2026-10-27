# spGARCH

R/C++ package for spatial and spatiotemporal ARCH/GARCH models. It builds on Philipp Otto's spGARCH 0.2.2 (CRAN) and adds sparse C++ simulation, spatial GARCH estimation, spatiotemporal models and parallel batch simulation.

## Contents

### Simulation (`R/Simulations.R`, `src/Simulations.cpp`)
| Function | Models | Notes |
|---|---|---|
| `sim.spARCH()` | `spARCH`, `log-spARCH`, `e-spARCH`, `complex-spARCH` | Sparse LU solves in C++ |
| `sim.spGARCH()` | `spGARCH`, `log-spGARCH`, `e-spGARCH`, `h-spGARCH`, `complex-spGARCH` | Analytical solve for e-spGARCH; O(nnz) input checks |
| `sim.spARCH_batch()`, `sim.spGARCH_batch()` | same as above | `n_sim` replications, parallel via `foreach`/`doRNG` with reproducible seeds |
| `sim.stARCH()`, `sim.stGARCH()` | spatiotemporal ARCH/GARCH | `control$burnin` defaults to 50 |

### Estimation (`R/Estimation.R`, `src/Estimation.cpp`)
| Function | Method | Notes |
|---|---|---|
| `qml.spARCH()` | QML for `spARCH`, `log-spARCH` | Numerically stable log-determinant |
| `qml.SARspARCH()` | QML for SAR model with spARCH errors | |
| `qml.spGARCH()` | QML for `spGARCH`, `e-spGARCH`, `log-spGARCH` | Covariates not supported yet |
| `nlse.spGARCH()` | Non-linear least squares | L-BFGS-B; data-driven start values if `start = NULL` |
| `gmm.stARCH()` | GMM for spatiotemporal ARCH | Online Helmert transform; `alpha` is not identified by design |

The fitted objects have the S3 methods `print`, `summary`, `plot`, `fitted`, `residuals`, `logLik` and `extractAIC` (class `spARCH`).

Data: `LondonCrime`, `LondonCrimeRate`, `LondonMap`, `prostate_cancer`.

### Model conventions
- At least two spatial locations are required (`n >= 2`).
- `alpha` is a scalar in `sim.spARCH` and `sim.spGARCH`.
- log-spARCH: ln h = alpha + rho W ln Y^2.

## Layout
```
R/, src/, man/, data/, inst/   package sources
tests/testthat/                test-simulations.R (interface), test-backend.R (C++ vs model equations)
examples/                      simulation_estimation.R (end-to-end workflow), London crime Rmds
dev/                           implementation_guide.txt (models, call chains, known issues), literature.md
```
Run the scripts in `examples/` from the package root.

## Development
```r
devtools::load_all()   # compile + load
devtools::test()       # run tests
devtools::document()   # regenerate man/ and NAMESPACE (roxygen)
```

Known open issues (details in `dev/implementation_guide.txt`, section 8):
- `qml.spGARCH` does not support covariates
- no `summary` method with proper labels for spGARCH fits
- `gmm.stARCH`: start values are hard-coded and the `p` argument is ignored (lag 1 only)
- model naming is inconsistent (`.LL_EspARCH` is log-spARCH)
- stationarity is enforced differently across model types

## References
Otto, Schmid & Garthoff (2018), *Spatial Statistics* 26, 125–145. Otto & Schmid (2019), arXiv:1908.08320. Otto (2019), *The R Journal*, arXiv:1812.01871.
