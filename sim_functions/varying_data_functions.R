
rgamma_trunc = function(n, shape, rate, lb = 0) {
  qgamma(runif(n) * pgamma(lb, shape, rate, lower = FALSE), shape, rate, lower = FALSE);
}


validate_bayesian_isotonic = function(x_validation_category, 
                                      true_prob_validation,
                                      draws_xi) {
  
  matrix_fitted_probs = 
    draws_xi[,x_validation_category,drop = FALSE];
  if(any(is.nan(matrix_fitted_probs))) {
    matrix_fitted_probs = 
      matrix_fitted_probs[which(rowSums(is.nan(matrix_fitted_probs)) == 0), , drop = FALSE]
  }
  if(nrow(matrix_fitted_probs)) {
  tiny_positive = .Machine$double.eps;
  matrix_fitted_probs = pmax(pmin(matrix_fitted_probs, 1 - tiny_positive), tiny_positive);
  
  matrix_true_probs = 
    matrix(rep(true_prob_validation, each = nrow(matrix_fitted_probs)), 
           nrow = nrow(matrix_fitted_probs));
  
  matrix_true_probs = pmax(pmin(matrix_true_probs, 1 - tiny_positive), tiny_positive);
  
  matrix_bias_probs = 
    matrix_true_probs - matrix_fitted_probs;
  
  curr_rmse = 
    sqrt(mean((matrix_bias_probs)^2));
  
  curr_bias = 
    mean(matrix_bias_probs);
  
  curr_loglik = 
    mean((matrix_true_probs * log(matrix_fitted_probs)) + 
           ((1 - matrix_true_probs) * log(1 - matrix_fitted_probs)) - 
           tiny_positive * log(tiny_positive));
  
  curr_kl_div = 
    mean((matrix_true_probs * log(matrix_true_probs)) + 
           ((1 - matrix_true_probs) * log(1 - matrix_true_probs)) -
           tiny_positive * log(tiny_positive)) - 
    curr_loglik;
  
  curr_coverage = 
    mean((apply(matrix_bias_probs, 2, quantile, p = 0.25) < 0) &
           ((apply(matrix_bias_probs, 2, quantile, p = 0.75) > 0)))
  
  summary_fitted_isotonic = 
    tibble(rmse = curr_rmse, 
           bias = curr_bias,
           loglik = curr_loglik,
           kl_div = curr_kl_div,
           coverage = curr_coverage);
  
  returned_predictions =
    tibble(x_cat = x_validation_category, 
           true_prob = true_prob_validation,
           fitted_prob = apply(matrix_fitted_probs, 2, median),
           lower50 = apply(matrix_fitted_probs, 2, quantile, p = 0.25),
           upper50 = apply(matrix_fitted_probs, 2, quantile, p = 0.75))
  } else {
    summary_fitted_isotonic = 
      tibble(rmse = NA, 
             bias = NA,
             loglik = NA,
             kl_div = NA,
             coverage = NA);
    
    returned_predictions =
      tibble(x_cat = x_validation_category, 
             true_prob = true_prob_validation,
             fitted_prob = NA,
             lower50 = NA,
             upper50 = NA)
  }
  
  list(summary_fitted_isotonic = summary_fitted_isotonic,
       returned_predictions = returned_predictions)
}


validate_isotonic = function(x_validation, 
                             true_prob_validation,
                             fitted_isotonic) {
  
  tiny_positive = .Machine$double.eps;
  
  fitted_isotonic <- 
    tibble(x_validation = x_validation, 
           true_prob_validation = true_prob_validation) %>%
    group_by(x_validation) %>%
    summarize(true_prob_validation = mean(true_prob_validation),
              number_validate = n()) %>%
    full_join(fitted_isotonic, 
              by = c("x_validation" = "x")) %>%
    arrange(x_validation) %>%
    fill(y, lower50conf, upper50conf, 
         .direction = "downup") %>%
    mutate(fitted_prob_validation =
             pmax(pmin(y, 1 - tiny_positive), tiny_positive),
           true_prob_validation =
             pmax(pmin(true_prob_validation, 1 - tiny_positive), tiny_positive), 
           bias_prob = 
             true_prob_validation - fitted_prob_validation,
           loglik = 
             (true_prob_validation * log(fitted_prob_validation)) + 
             ((1 - true_prob_validation) * log(1 - fitted_prob_validation)) - 
             tiny_positive * log(tiny_positive), 
           kl_div = 
             (true_prob_validation * log(true_prob_validation)) + 
             ((1 - true_prob_validation) * log(1 - true_prob_validation)) -
             tiny_positive * log(tiny_positive) - 
             loglik, 
           coverage = 
             (lower50conf <= true_prob_validation) * 
             (upper50conf >= true_prob_validation));
  
  
  summary_fitted_isotonic <- 
    fitted_isotonic %>%
    summarize(rmse = 
                sqrt(weighted.mean((bias_prob)^2, 
                                   w = number_validate)),
              bias = 
                weighted.mean(bias_prob, 
                              w = number_validate),
              loglik = 
                weighted.mean(loglik, 
                              w = number_validate),
              kl_div = 
                weighted.mean(kl_div, 
                              w = number_validate), 
              coverage = 
                weighted.mean(coverage, 
                              w = number_validate));
  
  list(fitted_isotonic = fitted_isotonic,
       summary_fitted_isotonic = summary_fitted_isotonic);
}




simulator = function(array_id = 1,
                     scenario_id = 1,
                     n_sim = 10,
                     n_training = 250,
                     n_validation = 2000,
                     true_prob_curve = function(x) {0.25 + 0.50 * (x > 0)}, 
                     true_prob_curve_id = 1,
                     predictor_dist = function(n) {sample(x = seq(0.05, 0.95, length = 10), size = n, replace = TRUE)},
                     predictor_dist_id = 1,
                     do_hs = TRUE,
                     hs_scales = c(horseshoe1 = 0.0025),
                     hs_slab_precisions = c(horseshoe1 = 1),
                     hs_prior_types = c("horseshoe"),
                     do_ga = TRUE,
                     ga_lower_bounds = c(gamma1 = .Machine$double.eps,
                                         gamma5 = 0),
                     ga_shapes = c(gamma1 = quote(0.5 / length_unique_x),
                                   gamma5 = 1),
                     mc_warmup = 2.5e3, 
                     mc_samps = 5e3, 
                     mc_chains = 1, 
                     mc_thin = 1, 
                     mc_stepsize = 0.1, 
                     mc_adapt_delta = 0.99,
                     mc_max_treedepth = 15,
                     ntries_per_iter = 1,
                     random_seed = sample.int(.Machine$integer.max,1),
                     data_seeds = NULL,#If non-null, should have length equal to n_sim
                     stan_seeds = NULL,#If non-null, should have length equal to n_sim
                     include_comparators = TRUE,
                     return_summarized_models = TRUE,
                     nonbayes_added_weight = c(0, 1/64), 
                     dynamic_run = T) 
{ 
  
  begin_all = Sys.time();
  set.seed(random_seed);
  if(!length(data_seeds)) {
    data_seeds = sample.int(.Machine$integer.max,n_sim);
  } else {
    if(length(data_seeds) != n_sim) {
      stop("'data_seeds' was provided but must have length equal to 'n_sim'");
    }  
  }
  if(!length(stan_seeds)) {
    stan_seeds = sample.int(.Machine$integer.max,n_sim);
  } else {
    if(length(stan_seeds) != n_sim) {
      stop("'stan_seeds' was provided but must have length equal to 'n_sim'");
    }  
  }
  
  # If both do_hs and do_ga are FALSE, the 
  # simulator will be a dry run with no actual methods fit
  
  if(!do_hs && !do_ga) {
    do_these_priors = NULL;
    fit_methods = include_comparators;
  } else {
    fit_methods = T;
    if(do_hs && is.null(names(hs_scales))) {
      names(hs_scales) = 
        paste0("horseshoe", 1:length(hs_scales));
    }
    if(do_hs) {
      hs_names = names(hs_scales)
    } else {
      hs_names = NULL
    }
    
    if(do_ga && length(ga_lower_bounds) != length(ga_shapes))  {
      stop("length of 'ga_lower_bounds' should be equal to length 
             of 'ga_shapes'");
    }
    if(do_ga && !all(names(ga_lower_bounds) == names(ga_shapes)))  {
      stop("names of 'ga_lower_bounds' should be identically equal to names 
             of 'ga_shapes'");
    }
    
    if(do_ga && is.null(names(ga_lower_bounds))) {
      names(ga_lower_bounds) = 
        names(ga_shapes) = 
        paste0("gamma",1:length(ga_lower_bounds));
    }
    
    if(do_ga) {
      ga_names = names(ga_lower_bounds)
    } else {
      ga_names = NULL
    }
    
    do_these_priors = c(hs_names, ga_names);
    
  }
  
  if(include_comparators) {
    nonbayes_names = paste0("nonbayes", seq_along(nonbayes_added_weight));
    names(nonbayes_added_weight) = nonbayes_names;
    do_these_priors = c(nonbayes_names, "brms_mo", do_these_priors);
    brm_fit = brm(y | trials(n) ~ mo(x_cat), 
                  data =   tibble(n = c(2, 2, 2, 2, 2), 
                                  y = c(0, 1, 1, 1, 2), 
                                  x_cat = 1:5),
                  family = "binomial", 
                  chains = mc_chains,
                  iter = mc_samps + mc_warmup,
                  warmup = mc_warmup,
                  control = list(adapt_delta = mc_adapt_delta, 
                                 max_treedepth = mc_max_treedepth));
  } else {
    nonbayes_added_weight = NA;
  }
  
  tiny_positive = .Machine$double.eps;
  do_these_priors = factor(do_these_priors) %>% fct_inorder
  
  summarized_performance = 
    crossing(priors = do_these_priors, 
             crossing(sim_id = 1:n_sim,
                      array_id = array_id,
                      scenario_id = scenario_id, 
                      true_prob_curve_id = true_prob_curve_id, 
                      predictor_dist_id = predictor_dist_id,
                      n_training = n_training,
                      # Were sparse categories observed?
                      obs_sparse_group = NA_real_,
                      # What is probability of sparse categories?
                      prob_sparse_group = NA_real_,
                      #root mean squared error of pointwise probabilities:
                      rmse = NA_real_,
                      #bias using pointwise probabilities:
                      bias = NA_real_,
                      #oos pointwise log-likelihood using model:
                      loglik = NA_real_,
                      #oos pointwise kl divergence:
                      kl_div = NA_real_,
                      #oos log-likelihood using empiric mean of y:
                      loglik_intercept = NA_real_,
                      #oos 50% coverage 
                      coverage_50 = NA_real_));
  
  summarized_bayesian_performance = 
    crossing(priors = do_these_priors, 
             crossing(sim_id = 1:n_sim,
                      array_id = array_id,
                      scenario_id = scenario_id, 
                      true_prob_curve_id = true_prob_curve_id, 
                      predictor_dist_id = predictor_dist_id,
                      n_training = n_training,
                      # Were sparse categories observed?
                      obs_sparse_group = NA_real_,
                      # What is probability of sparse categories?
                      prob_sparse_group = NA_real_,
                      #Number of breaks in the predictor used:
                      n_breaks = NA_real_,
                      #Either shape parameter (gamma prior) or scale parameter (horseshoe prior):
                      tuning_param_val = NA_real_,
                      #Bayesian RMSE
                      rmse = NA_real_,
                      #Bayesian bias
                      bias = NA_real_,
                      #Bayesian oos log-likehood
                      loglik = NA_real_,
                      #Bayesian KL divergence:
                      kl_div = NA_real_,
                      #Bayesian coverage
                      coverage_50 = NA_real_, 
                      #how many divergent transitions were there?
                      number_divergences = NA_real_,
                      #what was the value of the gelman-rubin diagnostic?
                      rhat = NA_real_, 
                      #where there NaNs? (another symptom of poor mixing)
                      any_nan = NA_real_,
                      #time required to fit each method
                      run_time_secs = NA_real_,
                      #ratio of run time of fastest versus slowest chains
                      chain_relative_diff = NA_real_)) %>%
    filter(!str_detect(priors, "nonbayes")) %>%
    mutate(priors = fct_drop(priors))
  
  
  if(return_summarized_models) {
    
    prespecified_x_validation = 
      predictor_dist(1e6) %>% unique() %>% sort();
    true_prob_prespecified_validation = 
      true_prob_curve(prespecified_x_validation);
    summarized_models = 
      crossing(priors = do_these_priors, 
               crossing(sim_id = 1:n_sim,
                        x = prespecified_x_validation, 
                        array_id = array_id, 
                        scenario_id = scenario_id, 
                        true_prob_curve_id = true_prob_curve_id, 
                        predictor_dist_id = predictor_dist_id,
                        n_training = n_training,
                        true_prob = NA_real_, 
                        fitted_prob = NA_real_,
                        lower50 = NA_real_,
                        upper50 = NA_real_));
  } else {
    summarized_models = NULL;
  }
  
  begin_sim = Sys.time();
  i=1;
  
  if(fit_methods) {
    cat("\n######################################\n");
    cat("# the true probability curve is defined by the function: \n");
    cat("# ")
    print(true_prob_curve);
    cat("# the true distribution of the predictor is defined by the function: \n");
    cat("# ")
    print(predictor_dist);
    cat("# the number of observations used for training is", n_training, "\n");
    cat("# the value of 'random_seed' is", random_seed,"\n");
    cat("# the value of 'n_sim' is", n_sim,"\n");
    cat("######################################\n\n");
  }
  
  for(i in 1:n_sim) {
    
    if(fit_methods) {
      
      # Draw data ----
      set.seed(data_seeds[i]);
      x = predictor_dist(n_training);
      length_unique_x = length(unique(x));
      true_prob = true_prob_curve(x);
      if(any(true_prob > 1 | true_prob < 0)) {
        stop(paste0("The function 'true_prob_curve' generated invalid probabilities, i.e. outside [0,1], given the following inputs: ", paste0(formatC(x[which(true_prob > 1 | true_prob < 0)],format = "f", digits = 4), collapse = ", ")));
      }
      y = rbinom(n_training, 1, prob = true_prob);
      mean_y = pmax(pmin(mean(y), 1 - tiny_positive), tiny_positive);
      
      x_validation = predictor_dist(n_validation);
      true_prob_validation = true_prob_curve(x_validation);
      if(any(true_prob_validation > 1 | true_prob_validation < 0)) {
        stop(paste0("The function 'true_prob_curve' generated invalid probabilities, i.e. outside [0,1], given the following inputs: ", paste0(formatC(true_prob_validation[which(true_prob_validation > 1 | true_prob_validation < 0)],format = "f", digits = 4), collapse = ", ")));
      }
      true_prob_validation = pmax(pmin(true_prob_validation, 1 - tiny_positive), tiny_positive)
      
      curr_loglik_intercept = 
        mean(true_prob_validation * log(mean_y) + 
               (1 - true_prob_validation)  * log(1 - mean_y) -
               tiny_positive * log(tiny_positive));
      
      # Non-Bayesian Isotonic Regression first ----
      if(include_comparators) {
        for(curr_prior in nonbayes_names) {
          
          cat("\n######################################\n");
          cat("# iteration:",i, "/", n_sim,"\n");
          cat("# the current data seed is",data_seeds[i],"\n");
          cat("#", curr_prior,"\n");
          cat("######################################\n\n");
          
          curr_row_performance = 
            with(summarized_performance, 
                 which(sim_id == i & 
                         priors == curr_prior));
          
          if(return_summarized_models) {
            curr_rows_model = 
              with(summarized_models, 
                   which(sim_id == i & 
                           priors == curr_prior));
          }
          
          
          data_grouped_iso = 
            tibble(x = x,
                   y = y) %>%
            group_by(x) %>%
            summarize(wt = length(y),
                      unweighted_y = mean(y), 
                      weighted_y = (sum(y) + (0.5 * nonbayes_added_weight[curr_prior])) /
                        (length(y) + nonbayes_added_weight[curr_prior])) %>%
            arrange(x)
          
          if(nrow(data_grouped_iso) > 1) {
            curr_fit <- 
              quickIsotone(doseResponse(y = pull(data_grouped_iso, unweighted_y), 
                                        x = pull(data_grouped_iso, x), 
                                        wt = pull(data_grouped_iso, wt)), 
                           estfun = oldPAVA, 
                           intfun = wilsonCI,
                           conf = 0.5);
          } else {
            # If just one category, reduces to Wilson's confidence interval, i.e.
            # inverted score
            foo = binom.confint(x = pull(data_grouped_iso, unweighted_y) * 
                                  pull(data_grouped_iso, wt),
                                n = pull(data_grouped_iso, wt), 
                                conf.level = 0.5, 
                                method = "wilson");
            curr_fit <-
              tibble(x = pull(data_grouped_iso, x), 
                     y = pull(data_grouped_iso, unweighted_y),
                     lower50conf = foo$lower[1],
                     upper50conf = foo$upper[1])
          }
          
          model_performance = 
            validate_isotonic(x_validation = x_validation, 
                              true_prob_validation = true_prob_validation, 
                              fitted_isotonic = curr_fit);
          
          summarized_performance[curr_row_performance, c("rmse","bias","loglik","kl_div","coverage_50")] = 
            model_performance$summary_fitted_isotonic[,c("rmse","bias","loglik","kl_div","coverage")];
          summarized_performance[curr_row_performance, "loglik_intercept"] = 
            curr_loglik_intercept;
          
          #Model summary ----
          if(return_summarized_models) {
            model_summary = 
              validate_isotonic(x_validation = prespecified_x_validation, 
                                true_prob_validation = true_prob_prespecified_validation, 
                                fitted_isotonic = curr_fit);
            
            summarized_models[curr_rows_model,c("x","true_prob","fitted_prob","lower50","upper50")] = 
              model_summary$fitted_isotonic[,c("x_validation","true_prob_validation","y","lower50conf","upper50conf")] %>%
              arrange(x_validation);
            rm(curr_rows_model, model_summary);
          }
          
          rm(data_grouped_iso, model_performance, curr_fit, curr_row_performance, curr_prior);
        }
      }
      
      if(length(c(hs_names, ga_names)) > 0) {
        
        # PHIL: DELETE THESE LINES BELOW ----
        # curr_prior = hs_names[1];
        # PHIL: DELETE THESE LINES ABOVE ----
        
        data_grouped = 
          tibble(x = x, 
                 y = y) %>%
          group_by(x) %>%
          summarize(n = length(x),
                    y = sum(y)) %>%
          arrange(x) %>%
          mutate(x_cat = 1:n())
        
        validation_data_ungrouped <- 
          left_join(
            tibble(x = x_validation) %>%
              mutate(row_num = 1:n()),
            data_grouped %>%
              select(x, x_cat)) %>%
          arrange(x) %>%
          fill(x_cat, 
               .direction = "downup") %>%
          arrange(row_num)
        
        prespecified_validation_data_ungrouped <- 
          left_join(
            tibble(x = prespecified_x_validation) %>%
              mutate(row_num = 1:n()),
            data_grouped %>%
              select(x, x_cat), 
            by = "x") %>%
          arrange(x) %>%
          fill(x_cat, 
               .direction = "downup") %>%
          arrange(row_num) 
        
        for(curr_prior in c(hs_names, ga_names)) {
          
          cat("\n######################################\n");
          cat("# iteration:",i, "/", n_sim,"\n");
          cat("# the current data seed is",data_seeds[i],"\n");
          cat("# the current stan seed is",stan_seeds[i],"\n");
          cat("#", curr_prior, "prior","\n");
          cat("######################################\n\n");
          
          curr_row_performance = 
            with(summarized_performance, 
                 which(sim_id == i & 
                         priors == curr_prior));
          
          curr_row_bayesian_performance = 
            with(summarized_bayesian_performance, 
                 which(sim_id == i & 
                         priors == curr_prior));
          
          if(return_summarized_models) {
            curr_rows_model = 
              with(summarized_models, 
                   which(sim_id == i & 
                           priors == curr_prior));
          }
          
          
          if(curr_prior %in% hs_names) {
            # horseshoe
            stan_args = 
              list(
                local_dof_stan = 1, 
                global_dof_stan = 1,
                alpha_scale_stan = eval(hs_scales[[curr_prior]]),
                slab_precision_stan = hs_slab_precisions[[curr_prior]]
              );
            prior_type = hs_prior_types[[curr_prior]];
            
            summarized_bayesian_performance[curr_row_bayesian_performance, "tuning_param_val"] = 
              stan_args$alpha_scale_stan;
            
          } else {
            # gamma
            stan_args = 
              list(
                alpha_shape_stan = eval(ga_shapes[[curr_prior]]),
                tiny_positive_stan = ga_lower_bounds[[curr_prior]]);
            prior_type = "gamma"
            
            summarized_bayesian_performance[curr_row_bayesian_performance, "tuning_param_val"] = 
              stan_args$alpha_shape_stan;
            
          }
          
          summarized_performance[curr_row_performance, "obs_sparse_group"] = 
            summarized_bayesian_performance[curr_row_bayesian_performance, "obs_sparse_group"] =
            any(pull(data_grouped, y) == 0 | 
                  pull(data_grouped, y) == pull(data_grouped, n))
          
          summarized_performance[curr_row_performance, "prob_sparse_group"] =  
            summarized_bayesian_performance[curr_row_bayesian_performance, "prob_sparse_group"] =
            tibble(x = x, 
                   true_prob = true_prob) %>%
            group_by(x) %>%
            summarize(prob_homogenous = 
                        prod(true_prob) + prod(1 - true_prob)) %>%
            summarize(prob_any_homogenous = 1 - prod(1 - prob_homogenous)) %>%
            pull(prob_any_homogenous)
          
          #Fit Stan model ----
          curr_fit = bayesian_isotonic(data_grouped = data_grouped,
                                       prior_type = prior_type, 
                                       stan_args = stan_args, 
                                       # 'conf_level' must remain hard coded for this simulator
                                       # because the function 'validate_bayesian_isotonic' assumes 
                                       # this level
                                       conf_level = 0.50, 
                                       conf_level_direction = "both",
                                       sample_from_prior_only = F,
                                       mc_warmup = mc_warmup, 
                                       mc_samps = mc_samps, 
                                       mc_chains = mc_chains, 
                                       mc_thin = mc_thin, 
                                       mc_stepsize = mc_stepsize, 
                                       mc_adapt_delta = mc_adapt_delta,
                                       mc_max_treedepth = mc_max_treedepth,
                                       verbose = T,
                                       stan_seed = stan_seeds[i]);
          
          #Non-Bayesian Model performance ----
          
          fitted_isotonic =
            bind_cols(x = pull(data_grouped, x), 
                      y = apply(curr_fit$all_draws$xi, 2, quantile, p = 0.50, na.rm = TRUE),  
                      lower50conf = apply(curr_fit$all_draws$xi, 2, quantile, p = 0.25, na.rm = TRUE), 
                      upper50conf = apply(curr_fit$all_draws$xi, 2, quantile, p = 0.75, na.rm = TRUE))              
          
          model_performance = 
            validate_isotonic(x_validation = validation_data_ungrouped$x, 
                              true_prob_validation = true_prob_validation, 
                              fitted_isotonic = fitted_isotonic);
          
          summarized_performance[curr_row_performance, 
                                 c("rmse","bias","loglik","kl_div","coverage_50")] = 
            model_performance$summary_fitted_isotonic[,c("rmse","bias","loglik","kl_div","coverage")];
          summarized_performance[curr_row_performance, "loglik_intercept"] = 
            curr_loglik_intercept;
          rm(fitted_isotonic, curr_row_performance, model_performance);
          
          # Bayesian model performance
          bayesian_model_performance = 
            validate_bayesian_isotonic(x_validation_category = validation_data_ungrouped$x_cat,
                                       true_prob_validation = true_prob_validation, 
                                       draws_xi = curr_fit$all_draws$xi);
          
          
          summarized_bayesian_performance[curr_row_bayesian_performance, 
                                          c("rmse","bias","loglik","kl_div","coverage_50")] = 
            bayesian_model_performance$summary_fitted_isotonic[c("rmse","bias","loglik","kl_div","coverage")];
          
          
          summarized_bayesian_performance[curr_row_bayesian_performance, "number_divergences"] =
            curr_fit$number_divergences;
          summarized_bayesian_performance[curr_row_bayesian_performance, "rhat"] = curr_fit$max_rhat;
          summarized_bayesian_performance[curr_row_bayesian_performance, "any_nan"] =
            curr_fit$any_nan;
          summarized_bayesian_performance[curr_row_bayesian_performance, "run_time_secs"] = 
            curr_fit$total_run_time_secs;
          summarized_bayesian_performance[curr_row_bayesian_performance, "chain_relative_diff"] = 
            min(curr_fit$chain_run_times_secs) / max(curr_fit$chain_run_times_secs);
          
          #Model summary ----
          if(return_summarized_models) {
            model_summary = 
              validate_bayesian_isotonic(
                x = prespecified_validation_data_ungrouped$x_cat, 
                true_prob = true_prob_prespecified_validation, 
                draws_xi = curr_fit$all_draws$xi);
            
            summarized_models[curr_rows_model,
                              c("x","true_prob","fitted_prob","lower50","upper50")] = 
              bind_cols(model_summary$returned_predictions[,c("x_cat","true_prob","fitted_prob","lower50","upper50")],
                        x = prespecified_validation_data_ungrouped$x) %>%
              select(x, true_prob, fitted_prob, lower50, upper50) %>%
              arrange(x);
            rm(curr_rows_model, model_summary);
          }
          rm(bayesian_model_performance, curr_row_bayesian_performance);
          rm(curr_fit, stan_args);
        }
        rm(curr_prior);
      }
      # BRMS with monotonic effects ----
      if(include_comparators) {
        curr_prior = "brms_mo"
        cat("\n######################################\n");
        cat("# iteration:",i, "/", n_sim,"\n");
        cat("# the current data seed is",data_seeds[i],"\n");
        cat("#", curr_prior,"\n");
        cat("######################################\n\n");
        
        curr_row_performance = 
          with(summarized_performance, 
               which(sim_id == i & 
                       priors == curr_prior));
        
        curr_row_bayesian_performance = 
          with(summarized_bayesian_performance, 
               which(sim_id == i & 
                       priors == curr_prior));
        
        if(return_summarized_models) {
          curr_rows_model = 
            with(summarized_models, 
                 which(sim_id == i & 
                         priors == curr_prior));
        }
        
        if(length(c(hs_names, ga_names)) == 0) {
          data_grouped = 
            tibble(x = x, 
                   y = y) %>%
            group_by(x) %>%
            summarize(n = length(x),
                      y = sum(y)) %>%
            arrange(x) %>%
            mutate(x_cat = 1:n())
          
          validation_data_ungrouped <- 
            left_join(
              tibble(x = x_validation) %>%
                mutate(row_num = 1:n()),
              data_grouped %>%
                select(x, x_cat)) %>%
            arrange(x) %>%
            fill(x_cat, 
                 .direction = "downup") %>%
            arrange(row_num)
          
          prespecified_validation_data_ungrouped <- 
            left_join(
              tibble(x = prespecified_x_validation) %>%
                mutate(row_num = 1:n()),
              data_grouped %>%
                select(x, x_cat), 
              by = "x") %>%
            arrange(x) %>%
            fill(x_cat, 
                 .direction = "downup") %>%
            arrange(row_num) 
        }
        summarized_performance[curr_row_performance, "obs_sparse_group"] = 
          summarized_bayesian_performance[curr_row_bayesian_performance, "obs_sparse_group"] =
          any(pull(data_grouped, y) == 0 | 
                pull(data_grouped, y) == pull(data_grouped, n))
        
        summarized_performance[curr_row_performance, "prob_sparse_group"] =  
          summarized_bayesian_performance[curr_row_bayesian_performance, "prob_sparse_group"] =
          tibble(x = x, 
                 true_prob = true_prob) %>%
          group_by(x) %>%
          summarize(prob_homogenous = 
                      prod(true_prob) + prod(1 - true_prob)) %>%
          summarize(prob_any_homogenous = 1 - prod(1 - prob_homogenous)) %>%
          pull(prob_any_homogenous)
        
        #Fit Stan model ----
        curr_fit = 
          update(brm_fit, newdata = data_grouped, seed = stan_seeds[i]);
        
        posterior_xi = 
          (fitted(curr_fit, summary = F)  / matrix(data_grouped$n, nrow = 2 * mc_samps, ncol = nrow(data_grouped), byrow = T)) 
        posterior_median_xi = 
          posterior_xi  %>% apply(2, median)
        
        #Non-Bayesian Model performance ----
        
        fitted_isotonic =
          bind_cols(x = pull(data_grouped, x), 
                    y = posterior_median_xi,  
                    lower50conf = apply(posterior_xi, 2, quantile, p = 0.25, na.rm = TRUE), 
                    upper50conf = apply(posterior_xi, 2, quantile, p = 0.75, na.rm = TRUE))              
        
        model_performance = 
          validate_isotonic(x_validation = validation_data_ungrouped$x, 
                            true_prob_validation = true_prob_validation, 
                            fitted_isotonic = fitted_isotonic);
        
        summarized_performance[curr_row_performance, 
                               c("rmse","bias","loglik","kl_div","coverage_50")] = 
          model_performance$summary_fitted_isotonic[,c("rmse","bias","loglik","kl_div","coverage")];
        summarized_performance[curr_row_performance, "loglik_intercept"] = 
          curr_loglik_intercept;
        rm(fitted_isotonic, curr_row_performance, model_performance);
        
        # Bayesian model performance
        bayesian_model_performance = 
          validate_bayesian_isotonic(x_validation_category = validation_data_ungrouped$x_cat,
                                     true_prob_validation = true_prob_validation, 
                                     draws_xi = posterior_xi);
        
        
        summarized_bayesian_performance[curr_row_bayesian_performance, 
                                        c("rmse","bias","loglik","kl_div","coverage_50")] = 
          bayesian_model_performance$summary_fitted_isotonic[c("rmse","bias","loglik","kl_div","coverage")];
        
        
        summarized_bayesian_performance[curr_row_bayesian_performance, "number_divergences"] =
          sum(get_divergent_iterations(curr_fit$fit));
        summarized_bayesian_performance[curr_row_bayesian_performance, "rhat"] = max(rhat(curr_fit));
        summarized_bayesian_performance[curr_row_bayesian_performance, "any_nan"] =
          any(is.nan(posterior_xi));
        curr_fit_runtimes <- rowSums(get_elapsed_time(curr_fit$fit))
        summarized_bayesian_performance[curr_row_bayesian_performance, "run_time_secs"] = 
          max(curr_fit_runtimes);
        summarized_bayesian_performance[curr_row_bayesian_performance, "chain_relative_diff"] = 
          min(curr_fit_runtimes) / max(curr_fit_runtimes);
        
        #Model summary ----
        if(return_summarized_models) {
          model_summary = 
            validate_bayesian_isotonic(
              x = prespecified_validation_data_ungrouped$x_cat, 
              true_prob = true_prob_prespecified_validation, 
              draws_xi = posterior_xi);
          
          summarized_models[curr_rows_model,
                            c("x","true_prob","fitted_prob","lower50","upper50")] = 
            bind_cols(model_summary$returned_predictions[,c("x_cat","true_prob","fitted_prob","lower50","upper50")],
                      x = prespecified_validation_data_ungrouped$x) %>%
            select(x, true_prob, fitted_prob, lower50, upper50) %>%
            arrange(x);
          rm(curr_rows_model, model_summary);
        }
        rm(bayesian_model_performance, curr_row_bayesian_performance);
        rm(model_performance)
        rm(curr_fit_runtimes,posterior_xi, posterior_median_xi )
        rm(curr_fit);
        
      }
      
      rm(data_grouped, x, true_prob, y, mean_y, 
         prespecified_validation_data_ungrouped,
         x_validation, validation_data_ungrouped, 
         true_prob_validation, curr_loglik_intercept);
      
    } else if(i == 1) {
      cat("Conducting dry run only. No simulations were run.\n");
    }
  }
  
  # +++ Report out results ----  
  
  list(summarized_performance = summarized_performance,
       summarized_bayesian_performance = summarized_bayesian_performance,
       summarized_models = summarized_models, 
       nonbayes_added_weight = nonbayes_added_weight,
       hs_scales = hs_scales,
       hs_prior_types = hs_prior_types,
       ga_lower_bounds = ga_lower_bounds,
       ga_shapes = ga_shapes,
       n_training = n_training,
       n_validation = n_validation,
       true_prob_curve = true_prob_curve,
       predictor_dist = predictor_dist, 
       random_seed = random_seed,
       data_seeds = data_seeds,
       stan_seeds = stan_seeds,
       total_run_time_secs = difftime(Sys.time(), begin_all, units = "secs"));
} 


