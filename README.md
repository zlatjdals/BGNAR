# BGNAR replication code

This archive contains the simulation, wind-data, and reporting code used for
the final BGNAR analysis.  Generated checkpoints, logs, tables, figures, and
private machine paths are intentionally excluded.

## Default experiment

The default configuration is stored in `config.env` and matches the final
reported analysis:

- DGPs: `M11`, `M12`, `M13`, and `M2`;
- graphs: Erdős-Rényi (`ER`), four-block stochastic block model (`SBM`), and
  small-world network (`SWN`);
- training lengths: 20, 40, 80, and 160;
- 100 replications per cell, 20 nodes, and five-step forecasts;
- common intercept 0.1 and innovations with standard deviation
  `sqrt(0.1)`;
- BGNAR fixed envelope `p=5`, `s=(3,3,3,3,3)` with hierarchical own-lag
  effects and `kappa_alpha ~ IG(3,2)`;
- reduced SOC plus DIO with adaptive, scale-calibrated tightness
  (`a_tau=2.5`, `b_tau=1`);
- centered GNAR with BIC selection over temporal lags up to 5 and network
  orders up to 3;
- BVAR fitted at the fixed lag order `p=5` (no BIC selection).

The coefficient designs are

| DGP | own-lag coefficients | network coefficients |
|---|---|---|
| M11 | `alpha=(0.70)` | `beta_1=(0.20)` |
| M12 | `alpha=(0.45)` | `beta_1=(0.45)` |
| M13 | `alpha=(0.20)` | `beta_1=(0.70)` |
| M2 | `alpha=(0.30,0.20)` | `beta_1=(0.20,0.10)`, `beta_2=(0.10,0)` |

For the graphs, ER uses edge probability 0.2.  SBM uses four equal blocks,
within-block probability 0.70, and between-block probability 0.075.  SWN uses
a one-dimensional ring with two neighbors on each side and rewiring
probability 0.1.  Disconnected draws are rejected.

## File order

The two main model files retain their established names and order:

1. `01_bgnar_model.R` - BGNAR Gibbs sampler, dummy observations, prediction,
   and posterior summaries.
2. `02_compare_methods.R` - wrappers for BGNAR, fixed-order BVAR,
   BIC-selected GNAR, and random walk.

The remaining files are:

- `00_utils.R` - graph, design-matrix, forecast, and evaluation utilities;
- `00_check_environment.R` - package and R-version check;
- `03_cluster_functions.R` - DGPs, graph generators, configuration, and
  method dispatch;
- `04_run_one.R` - one simulation cell;
- `05_collect_results.R` - simulation checkpoint aggregation;
- `06_validate_results.R` - simulation grid and failure validation;
- `07_run_parallel.sh` - parallel simulation launcher;
- `08_collect_only.sh` - simulation collection without refitting;
- `09_wind_functions.R` - CAPEL-centered 50-station wind design;
- `10_run_wind_one.R` - one rolling forecast origin;
- `11_collect_wind_results.R` - wind aggregation and validation;
- `12_run_wind_parallel.sh` - parallel wind launcher;
- `13_collect_wind_only.sh` - wind collection without refitting;
- `14_make_simulation_tables.R` - simulation reporting tables;
- `15_make_wind_outputs.R` - wind reporting table and basic figure;
- `16_make_manuscript_figures.R` - final simulation and wind figures;
- `config.env` - all user-editable defaults.

## Requirements

R packages: `igraph`, `GNAR`, and `BVAR`.

```bash
Rscript 00_check_environment.R
```

To install missing packages into the archive-local `.Rlib` directory:

```bash
INSTALL_MISSING=1 Rscript 00_check_environment.R
```

All scripts determine their directory at runtime.  No user-specific or server
path is hard-coded.

## Simulation

```bash
chmod +x 07_run_parallel.sh 08_collect_only.sh
nohup ./07_run_parallel.sh > simulation_master.out 2> simulation_master.err &
```

Before a full run, inspect the job grid without fitting:

```bash
DRY_RUN=1 DGP_LIST=M11 GRAPH_LIST=ER TRAIN_GRID=20 REPLICATIONS=2 \
  ./07_run_parallel.sh
```

One small fitting check can be run as follows:

```bash
source ./config.env
METHODS=bgnar_soc_dio_adaptive,gnar \
MCMC_ITER=200 MCMC_BURN=100 MCMC_THIN=1 \
Rscript 04_run_one.R M11 ER 20 1
```

Existing checkpoints are skipped.  Use `OVERWRITE=1` only when intentionally
replacing them.  To aggregate existing checkpoints without refitting:

```bash
./08_collect_only.sh
Rscript 14_make_simulation_tables.R
```

## Wind application

The wind analysis uses `GNAR::vswind`, chooses the 50-node connected
subnetwork by breadth-first search from CAPEL, and retains raw observations
for BGNAR and BVAR.  Only GNAR uses training-window nodewise centering.  The
default rolling design has a 240-observation window, 15 origins beginning at
541 and spaced by 12, and a five-step forecast horizon.

```bash
chmod +x 12_run_wind_parallel.sh 13_collect_wind_only.sh
nohup ./12_run_wind_parallel.sh > wind_master.out 2> wind_master.err &
```

To rebuild summaries without refitting:

```bash
./13_collect_wind_only.sh
Rscript 15_make_wind_outputs.R
```

After simulation and wind summaries both exist under `results/`, generate the
final manuscript figures with:

```bash
Rscript 16_make_manuscript_figures.R
```

Optional positional arguments are simulation result root, figure output
directory, and wind result root, in that order.

## Outputs

Each parallel job writes one atomic RDS checkpoint.  Collectors write CSV
summaries under `results/summary` and `results/wind/summary`.  Logs, full fits,
and generated results are excluded from version control by `.gitignore`.

Additional BGNAR dummy variants and random walk are available through the
`METHODS` setting; see the documented list in `config.env`.
