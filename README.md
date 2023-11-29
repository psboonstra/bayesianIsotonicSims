# `bayesianIsotonic` Numerical Studies
Code to accompany Boonstra, Owen, and Kang (2023)

This repository accompanies "Shrinkage priors for isotonic probability vectors
and binary data modeling" by Philip S. Boonstra, Daniel R. Owen, and Jian Kang
(2023). It allows for reproducing the numerical studies presented in Section 3
of that manuscript. 

You should be familiar with the R package `isotonicBayes`, which implements the
methodology. you can install via
`remotes::install_github("psboonstra/isotonicBayes", build_vignettes = TRUE)`.
You can read the introduction to the R functions in by running
`vignette("introduction", "isotonicBayes")` and you can see the code for
conducting Tutorial 1: Dose Escalation Example (Section 2.3 of the paper) by
running `vignette("tutorial1", "isotonicBayes")`.

There are five main parts to this

1. `fixed_data_evaluation.Rmd` / `fixed_data_evaluation.html` /
`fixed_data_evaluation.R` provide commands for conducting the fixed-data
evaluation in Section 3.1. Read through the documented code on
`fixed_data_evaluation.html` (which is the knitted version of
`fixed_data_evaluation.Rmd`), and when you are ready to run it your self, use
`fixed_data_evaluation.R` (which is the result of
`knitr::purl("fixed_data_evaluation.Rmd")`). This takes quite a while ($>1$ day)
to run start-to-finish on a personal computer. You can speed this up by
distributing the tasks.

2. `varying_data_run_sims.Rmd` / `varying_data_run_sims.html` /
`varying_data_run_sims.R` provide commands for conducting the varying-data
evaluation in Section 3.2, analogous to the previous item. Each instance of the
script can generate and analyze an arbitrary number of simulated datasets for
one of the eight data generating mechanisms presented in the manuscript (four
probability curves crossed with two sample sizes). However, for some of the data
generating mechanisms it takes many minutes or even hours to analyze a single
simulated dataset. Thus, although you can run this code locally on your own
computer, if you wish to conduct the full simulation study from the manuscript,
you will need to use an high-performance compute cluster and distribute multiple
instances of this script in parallel. The file `slurm_template.txt` gives the
skeleton of the SLURM batch script we used; you will need to fill in the
particulars before you can use it.


3. Having successfully run `varying_data_run_sims.R` and saved all of the outputs into a
folder called `out` (which will be done automatically above), read through and
run `varying_data_process_sims.Rmd` / `varying_data_process_sims.html` /
`varying_data_process_sims.R` to turn the raw results into the tables and figures in
the manuscript. You can do this step on your local computer, but you obviously
first need to download `out`.

4. The folder `sim_functions` contains R scripts pertaining to running the
varying-data evaluation.  It does not need to be directly called by the user. 

