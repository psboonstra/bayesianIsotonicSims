## ----setup, include = TRUE, echo = FALSE---------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)
klippy::klippy(position = c("top", "right"), tooltip_message = "Copy")


## ---- message = F, warning = F-------------------------------
library(tidyverse);
library(Iso);
library(cir);
library(binom);
library(isotonicBayes);
library(brms);
library(rstan);


## ------------------------------------------------------------
#Helpful to print warnings when they occur for debugging
options(warn = 1);
# Other options
rstan_options(auto_write = FALSE);
options(mc.cores = parallel::detectCores());


## ------------------------------------------------------------
my_computer = FALSE  
which_batch = 1 # Choose from 1, 2, 3, 4


## ------------------------------------------------------------
jobs_per_scenario = 100;
n_sim = 2;
if(my_computer) {
  #Choose from between 1-400 if varying_data_generate_params.R is used as-is
  array_id = 2;
} else {
  array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));
}


## ------------------------------------------------------------
permute_array_ids = TRUE


## ------------------------------------------------------------
array_id_offset = (which_batch - 1) * (jobs_per_scenario * 4);

if(which_batch %in% c(1,3)) {#first/third batches
  
  # Request 59 minutes running time for this batch (using 2 nodes)

  # This is the main simulation study;
  do_hs = TRUE;
  hs_scales = c(
    horseshoe1 = 0.001,
    horseshoe2 = 10)
  hs_prior_types = c(
    horseshoe1 = "horseshoe",
    horseshoe2 = "horseshoe")
  hs_slab_precisions = c(
    horseshoe1 = 1,
    horseshoe2 = 1)
  do_ga = TRUE;
  ga_shapes = c(gamma1 = quote(0.5 / length_unique_x),
                gamma5 = 1)
  ga_lower_bounds = c(gamma1 = .Machine$double.eps,
                      gamma5 = 0);
  
  include_comparators = TRUE;
  
} else if(which_batch %in% c(2,4)) {#second/fourth batches
  
  # Request 2 days running time for this batch (using 1 node)
  
  # This batch studies the effect of pushing the gamma lower bound closer to zero
  do_hs = FALSE
  do_ga = TRUE;
  ga_shapes = c(gamma2 = quote(0.5 / length_unique_x),
                gamma3 = quote(0.5 / length_unique_x),
                gamma4 = quote(0.5 / length_unique_x))
  ga_lower_bounds = c(gamma2 = 10^(log10(.Machine$double.eps) - 1), 
                      gamma3 = 10^(log10(.Machine$double.eps) - 2), 
                      gamma4 = 0);
  include_comparators = FALSE;
  
} else {
  stop("'which_batch' must be in 1,2,3,4")
}

if(which_batch %in% c(1,2)) {#first/second batches
  true_prob_curve_list = list(
    # Curve 1: simple linear increase
    `1` = function(x) {pmin(1, pmax(0, x));},
    # Curve 2: Two sharp increases
    `2` = function(x) {
      0.35 * plogis(x, location = 0.25, scale = 0.015) +
        0.65 * plogis(x, location = 0.75, scale = 0.015);
    }
  )
} else {#third/fourth batches
  true_prob_curve_list = list(
    # Curve 3: One sharp increase
    `3` = function(x) {
      .3 * plogis(x, location = 0.25, scale = 0.015);
    },
    # Curve 4: Smooth, positive first and second derivatives
    `4` = function(x) {pmin(1, pmax(0, x^4));}
  )
} 


## ------------------------------------------------------------
source("sim_functions/varying_data_generate_params.R");
source("sim_functions/varying_data_functions.R");


## ------------------------------------------------------------
if(permute_array_ids) {
  permute_array_ids = seq(1, jobs_per_scenario * length(arglist), by = jobs_per_scenario)
  set.seed(2);
  permute_array_ids = 
    c(permute_array_ids, 
      setdiff(sample(jobs_per_scenario * length(arglist)), permute_array_ids));
  
} else {
  permute_array_ids = 1:(jobs_per_scenario * length(arglist));
}


## ------------------------------------------------------------
# It's wasteful, but we only need the single id for this job
curr_args = arglist[[ceiling(permute_array_ids[array_id]/jobs_per_scenario)]];
# Ensure that the datasets are identical within array_ids and different 
# between array_ids
curr_args[["random_seed"]] = array_id;

# This is the actual call to the simulator function
assign(paste0("job",array_id_offset + array_id),
       do.call("simulator",args = curr_args));

# Save the entire workspace in case you want to dig in
do.call("save",list(paste0("job",array_id_offset + array_id),
                    file = paste0("out/job",array_id_offset + array_id,".RData"),
                    precheck = FALSE));
# Also just save the performance as a csv 
write_csv(get(paste0("job",array_id_offset + array_id))$summarized_performance,
          path = paste0("out/job",array_id_offset + array_id,"_performance.csv"),
          append = FALSE);

write_csv(get(paste0("job",array_id_offset + array_id))$summarized_bayesian_performance,
          path = paste0("out/job",array_id_offset + array_id,"_bayesian_performance.csv"),
          append = FALSE);

if(curr_args[["return_summarized_models"]]) {
  write_csv(get(paste0("job",array_id_offset + array_id))$summarized_models,
            path = paste0("out/job",array_id_offset + array_id,"_models.csv"),
            append = FALSE);
}


