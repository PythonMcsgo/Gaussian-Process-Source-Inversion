# Gaussian-Process-Source-Inversion
2D Poisson Equation Source Inversion by Gaussian Process

- `01_main_experiment/`：Main Experiement (Fixed Observations 7x7 Design), including the Calibration Plots
- `02_mcmc/`：6 Sources Adaptive RandomWalk MH MCMC;
- `03_active_learning/`：Active Learning (6 sources) including maximum-variance, and gradient-based approach.


For each experimetns, you can click the run.R file to realize the whole part.
```text
Rscript 01_main_experiment/run.R
Rscript 02_mcmc/run.R
Rscript 03_active_learning/run.R
```

## Note:
results. folder contains all the results and plots.
