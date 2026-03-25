# Joint Inference for Two-Phase Nonlinear Mixed Effects Models with Change Points, Censoring, and Measurement Errors

This repository contains the R code and data processing pipeline for the joint modeling of longitudinal viral load and CD4 cell counts using Non-Linear Mixed-Effects (NLME) models. The framework is designed to estimate latent change points in biomarker trajectories following treatment interruption.

## Repository Structure
The project is organized into five functional directories:
- `R/`: Contains core helper functions and the main implementation of the Joint NLME algorithm.
- `Data/`: Includes raw and processed datasets:
  - `raw_data.csv`: The original unedited dataset.
  - `cleaned_data.csv`: The master cleaned dataset.
  - `cleaned_data_decay`, `_rebound`, `_cd4`, `_trans`: Specialized subsets for different modeling phases.
- `Figures/`: Graphical outputs of the raw data.
  - `Viral load all patients.png`: Population level spaghetti plot.
  - `Viral load 3 patients.png`: Annotated comparison of selected trajectories.
  - `Ind_plot/`: Individual subject plots showing viral load and CD4 dynamics with potential change point ranges.
- `Data Analysis/`: Scripts for applying the model to the clinical data.
  - `analysis.R`: Main script to fit the Joint NLME model.
  - `plot_fitted.R`: Generates visualizations of the estimated trajectories.
  - `Results/`: Contains `Results.RDS`, `Table.tex` (for paper inclusion), and `FittedPlots/` showing estimated change points.
- `Simulation/`: Scripts and results for the simulation study.
  - `SimResults I/`: Performance metrics for the lower censoring rate setting.
  - `SimResults II/`: Performance metrics for the higher censoring rate setting.

## Getting Started
### Prerequisites
Ensure you have **R (>= 4.0.0)** installed. The following packages are required:

```{r}
install.packages(c("tidyverse", "patchwork", "scales", "here", "berryFunctions"))
```

### Installation
1. Clone this repository to your local machine.
2. Open the `.Rproj` file or set your working directory to the root folder. The project utilizes the `here` package for robust file path management.

## Usage Workflow
### 1. Data Preparation
Run the `DataProcess.R` script located in the root directory. This script performs:
- Outlier removal and data cleaning.
- Variable transformation (log-scaling).
- Identification of treatment interruption windows.

### 2. Exploratory Visualization
Run `Plot.R` to generate the visual summaries found in the `Figures/` folder, including individual trajectories for all patients.

### 3. Model Fitting
Navigate to `Data Analysis/` and run `analysis.R`. This will:
- Load the cleaned data.
- Execute the Joint NLME estimation.
- Save the outputs to `Data Analysis/Results/`.

### 4. Simulation Study
To reproduce the simulation results, run the scripts within the `Simulation/` folder. Results are partitioned by censoring rate settings to validate model.
