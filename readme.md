DEBtool_M is Matlab code that can be used for application of Dynamic Energy Budget (DEB) theory
https://bio.vu.nl/thb/deb/deblab/

add-my-pet.github.io/debtool_m/docs/

## DEBtool_M (Custom Fork & Zoo Collection)

Please note that the entirety of the code in this repository—with the strict exception of the zoo_collection directory—is a direct fork of DEBtool_M.

If you are looking to understand the core mechanics, modify the base mathematical models, or need guidance on standard DEBtool development, you must refer to the official documentation of the original repository.

Original Repository & Documentation: https://github.com/add-my-pet/DEBtool_M

## `zoo_collection` Documentation

The `zoo_collection` directory is the unique, custom addition to this fork.
It acts as an isolated workspace to house all datasets, parameter estimations,
and analysis pipelines specific to a curated collection of species, keeping
the core DEBtool codebase clean and untouched.

### Overview

This module is designed to manage the data and run the AmP (Add-my-Pet) 
estimation procedures specifically tailored for target animals.
It centralizes custom routines so that researchers can easily execute models,
analyze physiological traits, and plot results without navigating 
through the base `DEBtool_M` directories.

### Directory Structure & Contents

Here is a precise breakdown of what you will find inside the `zoo_collection` module:

*   **`data/` (Species Data Files):**
    * Contains the standard `mydata_*.m` files that hold the zero-variate (e.g., age at birth, maximum weight) and uni-variate data (e.g., length-weight curves) for each specific species in our collection.
    * Includes supplemental reference data and environmental variation logs (such as CSV files tracking food availability variations).
*   **`pars/` (Parameter Initializations):**
    * Houses the `pars_init_*.m` scripts. These files define the initial parameter guesses and set the pseudo-data required by the Nelder-Mead optimization algorithms to estimate standard DEB parameters.
*   **`scripts/` (Execution & Extended Analysis):**
    * **Run Scripts:** Contains the `run_*.m` execution files used to configure `estim_options` and launch `estim_pars` for the species in the collection.
    * **Post-Processing (R/Markdown):** Includes R-based analysis pipelines (e.g., `dataAnalysis.Rmd` and compiled HTML reports) designed to perform extended statistical analyses on the MATLAB outputs.
*   **`results/` (Model Outputs):**
    * The default output destination for all model runs. This folder stores the generated `.mat` files containing the final estimated parameters, implied traits, and the auto-generated plots comparing our model predictions against the empirical data.