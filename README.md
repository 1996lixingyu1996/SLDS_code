# M-survival Learner

This repository provides the R code for the **M-survival Learner** method.

Full paper link: https://arxiv.org/abs/2603.13464

## Requirements

Please ensure that the `xgboost` package version **1.7.9.1** is installed before running the code.  
Other required packages can be installed using standard R installation procedures.

## Real Data Analysis

The real data analysis can be executed directly without installing the R package by running the script: decision_tree.R  located in the “code” folder.

On a standard personal computer, the real data analysis requires approximately **5 minutes** to complete.

## Simulation Studies

The `code` folder contains scripts used for the simulation studies.

Each simulation scenario requires approximately **1 hour** to run on a standard personal computer.  
Due to the computational cost, users may adjust the number of simulation replicates if faster execution is desired.

## File Description

- `code/`  
  Contains scripts for simulation studies and data analysis.

- `code/HIV_example`  
  Script for the real data analysis, which produces the subgroup identification results reported in the manuscript.

- `code/decision_tree.R`  
  Main script for running the real data analysis without installing the package.

## Support and Contact

If you encounter any issues or would like to report bugs, please contact:

**Xingyu Li**  
Email: lxingyu1996@gmail.com or xli36@mdanderson.org

## Funding

This work was supported by **Amgen Inc.**

