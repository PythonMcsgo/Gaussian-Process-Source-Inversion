# Gaussian-Process-Source-Inversion
2D Poisson Equation Source Inversion by Gaussian Process

- `01_main_experiment/`：Main Experiement (Fixed Observations 7x7 Design), including the Calibration Plots
- `02_mcmc/`：6 Sources Adaptive RandomWalk MH MCMC;
- `03_active_learning/`：Active Learning (6 sources) including maximum-variance, and gradient-based approach.

The eaiast way to rull all experiment:
```text
Rscript RUN_ALL.R
```
MCMC takes long time to run, it is recommended to run each part independently.
```text
Rscript 01_main_experiment/run.R
Rscript 02_mcmc/run.R
Rscript 03_active_learning/run.R
```

## Note
- For `01_main_experiment/`, you should run `01_main_experiment/` first, then run 01b and 01c.
- `00_MCMC_core.R` cannot run independently.
