
################################################################################
# Creates a dataset with lagged variable values
# covariates_df: the dataset to be modified
# num_lags: number of lagged values to add to the dataset
# vars_to_lag: the variables whose values to lag

# ex: to include the values of 'T2D' for the past 3 intervals including the current interval (the default) then set vars_to_lag = c('T2D') and num_lags = 3

################################################################################

create_covariates_df <- function(covariates_df, num_lags, vars_to_lag) {
  
  if (num_lags == 0 | is.null(vars_to_lag)) {return(covariates_df)}
  
  covariates_df <- covariates_df %>% group_by(id)
  
  for (v in 1:length(vars_to_lag)) {
    
    for (i in 1:num_lags) {
      
      varname = paste0(vars_to_lag[v], '_', i)
      
      covariates_df <- covariates_df %>% mutate(., !!varname := lag(!!as.symbol(vars_to_lag[v]), i))
      
    }
  }
  
  return(covariates_df %>% ungroup)
  
}

################################################################################
# Computes cross-entropy / binary log loss
################################################################################

calc_logloss <- function(true, pred){
  
  pred <- pmax(pmin(1-1e-15, pred), 1e-15)
  
  logloss <- -sum(true*log(pred) + (1-true)*log(1-pred)) / length(pred)
  
  return(logloss)
}

################################################################################
# Computes squared error loss
################################################################################

calc_mse <- function(true, pred){
  
  mse <- sum((true-pred)^2) / length(pred)
  
  return(mse)
}

################################################################################
# Returns a vector which contains the cross validation `k`-fold
# identifier for each row in `data`
################################################################################

get_CV_ids <- function(data, k){
  
  folds_ids <- sample(cut(1:nrow(data), breaks = k, labels = FALSE),
                      size = nrow(data),
                      replace = FALSE)
  
  return(folds_ids)
  
}

################################################################################
# Returns the cross-validation loss for combinations of hyperparameter values for GBM
# covariates_df: the dataset containing `predictor` variables
# label_vector: the response/label variable
# nfold: number of cross-validation folds
# tree_depth: a vector of values for boosted tree depth
# shrinkage_factor: a vector of shrinkage values  
# num_trees: a vector of values for number of trees (iterations) 
# num_cores: number of cores to use for parallelization
# family: the loss function ('bernoulli' for binary log loss, otherwise squared error loss)
################################################################################

find_params_boosted_tree_model <- function(covariates_df, label_vector, nfold, tree_depth, 
                                           shrinkage_factor, num_trees, num_cores, family='bernoulli'){
  
  params <- expand.grid(n.trees = num_trees, interaction.depth = tree_depth, shrinkage = shrinkage_factor)
  
  cv_test_logloss <- vector(length=nrow(params))
  
  test_preds <- vector(length=nrow(covariates_df))
  
  cl <- makeCluster(num_cores)
  
  clusterEvalQ(cl, library(gbm))
  
  clusterExport(cl, c('params', 'nfold', 'label_vector', 
                      'covariates_df', 'get_CV_ids', 'calc_logloss', 'calc_mse', 'test_preds'),
                envir=environment())
  
  gbm_out <- 
    pbsapply(1:nrow(params), function(iter){
      
      cv_ids <- get_CV_ids(covariates_df, k=nfold)
      
      for (k in 1:nfold) {
        
        train_indices <- which(cv_ids != k)
        
        test_indices <- which(cv_ids == k)
        
        model <- 
          gbm(formula = Y~., data = data.frame(Y=label_vector[train_indices], covariates_df[train_indices, ,drop=FALSE]), 
              n.trees = params$n.trees[iter], interaction.depth = params$interaction.depth[iter], shrinkage = params$shrinkage[iter],
              distribution = family)
        
        test_preds[which(cv_ids == k)] <- predict(model, type='response', n.trees = model$n.trees,
                                                  newdata = covariates_df[test_indices, ,drop=FALSE])
        
      }
      
      if(family=='bernoulli'){calc_logloss(label_vector, test_preds)}else{calc_mse(label_vector, test_preds)}
      
    }, cl=cl)
  
  stopCluster(cl)
  
  best_params <- as.list(c(params[resample(which(gbm_out == min(gbm_out)), 1), ]))
  best_params$attempted_params <- params
  best_params$logloss <- gbm_out
  
  return(best_params)
  
}

################################################################################
# Returns predictions from a gbm model
# covariates_df: the dataset containing 'predictor' variables for fitting the gbm model
# label_vector: the response/label variable for fitting the gbm model
# params: a list containing a value for n.trees (number of gbm trees), interaction.depth (boosted tree depth), and shrinkage (the shrinkage value)
# predict_data: a new dataset for which predictions are desired
# boosted_tree_family: the loss function ('bernoulli' for binary log loss, otherwise squared error loss)
################################################################################

estimate_prob_boosted_tree_model <- function(covariates_df, label_vector,
                                             params, predict_data, boosted_tree_family = 'bernoulli'){
  
  boosted_tree_model <- 
    gbm(
      formula = Y~., 
      data = data.frame(Y=label_vector, covariates_df),
      interaction.depth = params$interaction.depth,
      shrinkage = params$shrinkage,
      n.trees = params$n.trees,
      distribution = boosted_tree_family
    )
  
  out <- predict(boosted_tree_model, newdata=predict_data, type="response", n.trees=boosted_tree_model$n.trees)
  
  
  return(out)
  
}

################################################################################
# Creates an indicator for whether a patient has followed the 
# medication initiation and grace period strategy up to the beginning of 
# each month.
# initiated_treatment: the medication initiated, ACE or thiazide
# grace_period_length: the number of months, m, by which a grace period is defined
# min_interval: defined as min(data$month) by resultsfnc()
# data: the cohort dataset in person-month format
################################################################################

create_following_strategy_var <- function(initiated_treatment,
                                          grace_period_length, 
                                          min_interval,
                                          data) {
  
  data <- data %>%
    group_by(id) %>%
    mutate(following_strategy = cumprod(
      case_when(num_intervals_without_trt <= grace_period_length & 
                  (initiated_med == initiated_treatment | month == min_interval) ~ 1,
                TRUE ~ 0)
    )) %>%
    ungroup()
  
  return(data)
  
}

################################################################################
# Returns inverse probability weights for the initiation of a medication at baseline.
# Can be used to construct weights for the initiation of one out of many possible medications
# min_interval: the index of the first time-interval of the study (the first month of treatment initiation in this data example)
# initiated_treatment: the medication to be initiated at baseline in the treatment strategy of interest
# model_formula_vars: a vector of the names of confounders to be included in the model for treatment initiation
# parametric_model: should a logistic GLM be used to compute probabilities of treatment initiation (T/F)
# data: the dataset
# nfold: if sample splitting is to be used with a non-parametric model for probability of treatment initiation, the number of folds to use for cross-validation
# tree_depth: if sample splitting is to be used with a non-parametric model for probability of treatment initiation, a vector of values for boosted tree depth
# shrinkage_factor: if sample splitting is to be used with a non-parametric model for probability of treatment initiation, a vector of shrinkage values
# num_trees: if sample splitting is to be used with a non-parametric model for probability of treatment initiation, a vector of values for number of trees (iterations)
# num_cores: if sample splitting is to be used with a non-parametric model for probability of treatment initiation, number of cores to use for parallelization
################################################################################

initiation_weightsfnc <- function(min_interval,
                                  initiated_treatment,
                                  model_formula_vars,
                                  parametric_model,
                                  data,
                                  nfold,
                                  tree_depth,
                                  shrinkage_factor,
                                  num_trees,
                                  num_cores){
  
  data <- data %>% 
    dplyr::select(c('id','month','has_med','dispensed','initiated_med',
                    'num_intervals_without_trt', 'following_strategy',
                    'censored', 'outcome','split_ids',
                    all_of(model_formula_vars))) %>%
    mutate(prob_trt = 1,
           initiated_med = case_when(initiated_med == initiated_treatment ~ 1,
                                     initiated_med != initiated_treatment ~ 0))
  
  data_k <- data %>% dplyr::select(-prob_trt) %>% 
    filter(month == min_interval)
  
  for (split in 1:length(unique(data$split_ids))) {
    
    if(parametric_model) {
      
      model_formula_RHS <- paste(
        names(data_k)[which(!names(data_k) %in% c('id','month','has_med','dispensed','initiated_med',
                                                  'num_intervals_without_trt', 'following_strategy',
                                                  'censored', 'outcome','split_ids'))], 
        collapse=" + "
      )
      
      if (length(unique(data$split_ids)) == 1) {
        
        model <- 
          glm(formula = as.formula(paste("initiated_med", model_formula_RHS, sep=" ~ ")),
              data = data_k,
              family = quasibinomial())
        
      } else {
        
        model <- 
          glm(formula = as.formula(paste("initiated_med", model_formula_RHS, sep=" ~ ")),
              data = data_k %>% filter(split_ids != split),
              family = quasibinomial())
        
      }
      
      prob_trt_k <- predict(model, newdata = data_k %>% filter(split_ids == split), type='response')
      
    } else {
      
      params_trt <- find_params_boosted_tree_model(covariates_df = data_k %>% 
                                                     filter(split_ids != split) %>% 
                                                     dplyr::select(-c('id','month','has_med','dispensed','initiated_med',
                                                                      'num_intervals_without_trt', 'following_strategy',
                                                                      'censored','outcome','split_ids')), 
                                                   label_vector = data_k %>% 
                                                     filter(split_ids != split) %>% 
                                                     dplyr::select(c(initiated_med)) %>% {.[[1]]},
                                                   nfold = nfold,
                                                   tree_depth = tree_depth,
                                                   shrinkage_factor = shrinkage_factor,
                                                   num_trees = num_trees,
                                                   family = 'bernoulli',
                                                   num_cores = num_cores)
      
      prob_trt_k <- estimate_prob_boosted_tree_model(covariates_df = data_k %>% 
                                                       filter(split_ids != split) %>% 
                                                       dplyr::select(-c('id','month','has_med','dispensed','initiated_med',
                                                                        'num_intervals_without_trt', 'following_strategy',
                                                                        'censored','outcome','split_ids')), 
                                                     label_vector = data_k %>% 
                                                       filter(split_ids != split) %>% 
                                                       dplyr::select(c(initiated_med)) %>% {.[[1]]},
                                                     params = params_trt,
                                                     boosted_tree_family = 'bernoulli',
                                                     predict_data = data_k %>% 
                                                       filter(split_ids == split) %>% 
                                                       dplyr::select(-c('id','month','has_med','dispensed','initiated_med',
                                                                        'num_intervals_without_trt', 'following_strategy',
                                                                        'censored','outcome','split_ids')))
      
    }
    
    data <- data %>% 
      group_by(id) %>%
      mutate(model_prob = case_when(month == min_interval & split_ids == split ~ 1,
                                    TRUE ~ 0)) %>%
      ungroup()
    
    data[which(data$model_prob == 1), ]$prob_trt <- prob_trt_k
    
    data <- data %>% 
      dplyr::select(-model_prob)
    
  }
  
  data <- data %>%
    mutate(initiation_weight = case_when(month == min_interval & initiated_med == 1 ~ 1/prob_trt,
                                         month == min_interval & initiated_med == 0 ~ 0,
                                         month != min_interval ~ 1))
  
  return(data$initiation_weight)
  
}

################################################################################
# Returns inverse probability weights for censoring due to loss to follow-up (and death if death is treated as a censoring event).
# min_interval: the index of the first time-interval of the study (the first month of treatment initiation in this data example)
# initiated_treatment: the medication to be initiated at baseline in the treatment strategy of interest
# model_formula_vars: a vector of the names of confounders to be included in the model for censoring
# parametric_model: should a logistic GLM be used to compute probabilities of censoring (T/F)
# data: the dataset
# nfold: if sample splitting is to be used with a non-parametric model for probability of censoring, the number of folds to use for cross-validation
# tree_depth: if sample splitting is to be used with a non-parametric model for probability of censoring, a vector of values for boosted tree depth
# shrinkage_factor: if sample splitting is to be used with a non-parametric model for probability of censoring, a vector of shrinkage values
# num_trees: if sample splitting is to be used with a non-parametric model for probability of censoring, a vector of values for number of trees (iterations)
# num_cores: if sample splitting is to be used with a non-parametric model for probability of censoring, number of cores to use for parallelization
################################################################################

censoring_weightsfnc <- function(min_interval,
                                 initiated_treatment,
                                 model_formula_vars,
                                 parametric_model,
                                 data,
                                 nfold,
                                 tree_depth,
                                 shrinkage_factor,
                                 num_trees,
                                 num_cores){
  
  data <- data %>% 
    dplyr::select(c('id','month','has_med','dispensed','initiated_med',
                    'num_intervals_without_trt', 'following_strategy',
                    'censored', 'outcome','split_ids',
                    all_of(model_formula_vars))) %>%
    mutate(prob_censored = 1)
  
  data_k <- data %>% dplyr::select(-prob_censored) %>% 
    filter(month >= min_interval)
  
  for (split in 1:length(unique(data$split_ids))) {
    
    if(parametric_model) {
      
      model_formula_RHS <- paste(
        names(data_k)[which(!names(data_k) %in% c('id','dispensed','initiated_med',
                                                  'following_strategy',
                                                  'outcome','split_ids'))], 
        collapse=" + "
      )
      
      if (length(unique(data$split_ids)) == 1) {
        
        model <- 
          glm(formula = as.formula(paste("censored", model_formula_RHS, sep=" ~ ")),
              data = data_k %>% filter(initiated_med == initiated_treatment),
              family = quasibinomial())
        
      } else {
        
        model <- 
          glm(formula = as.formula(paste("censored", model_formula_RHS, sep=" ~ ")),
              data = data_k %>% filter(split_ids != split & initiated_med == initiated_treatment),
              family = quasibinomial())
        
      }
      
      prob_censored_k <- predict(model, newdata = data_k %>% filter(split_ids == split), type='response')
      
    } else {
      
      params_censored <- find_params_boosted_tree_model(covariates_df = data_k %>% 
                                                          filter(split_ids != split & 
                                                                   initiated_med == initiated_treatment) %>% 
                                                          dplyr::select(-c('id','dispensed','initiated_med',
                                                                           'following_strategy',
                                                                           'outcome','split_ids')), 
                                                        label_vector = data_k %>% 
                                                          filter(split_ids != split &
                                                                   initiated_med == initiated_treatment) %>% 
                                                          dplyr::select(c(censored)) %>% {.[[1]]},
                                                        nfold = nfold,
                                                        tree_depth = tree_depth,
                                                        shrinkage_factor = shrinkage_factor,
                                                        num_trees = num_trees,
                                                        family = 'bernoulli',
                                                        num_cores = num_cores)
      
      prob_censored_k <- estimate_prob_boosted_tree_model(covariates_df = data_k %>% 
                                                            filter(split_ids != split &
                                                                     initiated_med == initiated_treatment) %>% 
                                                            dplyr::select(-c('id','dispensed','initiated_med',
                                                                             'following_strategy',
                                                                             'outcome','split_ids')), 
                                                          label_vector = data_k %>% 
                                                            filter(split_ids != split &
                                                                     initiated_med == initiated_treatment) %>% 
                                                            dplyr::select(c(censored)) %>% {.[[1]]},
                                                          params = params_censored,
                                                          boosted_tree_family = 'bernoulli',
                                                          predict_data = data_k %>% 
                                                            filter(split_ids == split) %>% 
                                                            dplyr::select(-c('id','dispensed','initiated_med',
                                                                             'following_strategy',
                                                                             'outcome','split_ids')))
      
    }
    
    data[which(data$split_ids == split), ]$prob_censored <- prob_censored_k
    
  }
  
  data <- data %>%
    mutate(censoring_weight = case_when(censored == 1 ~ 0,
                                        censored == 0 ~ 1/(1-prob_censored)))
  
  return(data$censoring_weight)
  
}

################################################################################
# Returns inverse probability weights for adherence according to an investigator specified grace period strategy.
# Can be used to construct weights for the natural or stochastic grace period strategies.
# skip_first_interval: whether the first interval of follow-up should be included (T/F). Generally, all individuals will receive treatment in the first interval in the observed data so this should be set to FALSE.
# f_int: a vector of probabilities defining the marginal probability of adherence in each grace period interval, should have length equal to grace period length
# f_cond_int: a vector of probabilities defining the conditional (on number of intervals without treatment) probability of adherence in each grace period interval, should have length equal to grace period length
# grace_period_length: the investigator specified grace period length. When the number of intervals without treatment equals this value, the conditional probability of treatment will be set to 1.

## To specify a natural grace period strategy, set f_int and f_cond_int to NULL. 
## To specify a stochastic grace period strategy, set values for either f_int or f_cond_int. f_cond_int will override f_int.
## ex: to define a stochastic grace period strategy with a uniform distribution of adherence with a grace period length of 3, set f_int = c(0.25, 0.25, 0.25) or f_cond_int = c(0.25, 0.33, 0.5), and grace_period_length = 3.

# min_interval: the index of the first time-interval of the study (the first month of treatment initiation in this data example)
# max_interval: the investigator specified final-interval of intervention
# initiated_treatment: the medication to be initiated at baseline in the treatment strategy of interest
# pooled_time_treatment_model: whether separate models should be fit at each time interval for the probability of adherence, or whether the model should pool over time intervals (TRUE to pool over time intervals, FALSE for separate models)
# model_formula_vars: a vector of the names of confounders to be included in the model for treatment adherence
# vars_to_lag: if pooled_time_treatment_model == F, then vars_to_lag can be used to specify, as a subset of models_formula_vars, the names of the variables to include histories of in the model for adherence
# num_lags: if pooled_time_treatment_model == F, then num_lags specifies the number of lagged histories for the variables in 'vars_to_lag' to include in the model for adherence

# ex: to include the values of 'T2D' for the past 3 intervals including the current interval (the default) then set vars_to_lag = c('T2D') and num_lags = 3

# parametric_model: should a logistic GLM be used to compute probabilities of adherence (T/F)
# data: the dataset
# nfold: if sample splitting is to be used with a non-parametric model for probability of adherence, the number of folds to use for cross-validation
# tree_depth: if sample splitting is to be used with a non-parametric model for probability of adherence, a vector of values for boosted tree depth
# shrinkage_factor: if sample splitting is to be used with a non-parametric model for probability of adherence, a vector of shrinkage values
# num_trees: if sample splitting is to be used with a non-parametric model for probability of adherence, a vector of values for number of trees (iterations)
# num_cores: if sample splitting is to be used with a non-parametric model for probability of adherence, number of cores to use for parallelization
################################################################################

adherence_weightsfnc <- function(skip_first_interval = F,
                                 f_int = NULL,
                                 f_cond_int = NULL,
                                 grace_period_length,
                                 min_interval,
                                 max_interval,
                                 initiated_treatment,
                                 pooled_time_treatment_model = F,
                                 model_formula_vars,
                                 vars_to_lag,
                                 num_lags,
                                 parametric_model,
                                 data,
                                 nfold,
                                 tree_depth,
                                 shrinkage_factor,
                                 num_trees,
                                 num_cores){
  
  data <- data %>% 
    dplyr::select(c('id','month','has_med','dispensed','initiated_med',
                    'num_intervals_without_trt', 'following_strategy',
                    'censored','outcome','split_ids',
                    all_of(model_formula_vars))) %>%
    mutate(prob_trt = 1)
  
  if (skip_first_interval) {min_interval <- min_interval+1}
  
  if (pooled_time_treatment_model) {
    
    data_k <- data %>% dplyr::select(-prob_trt) %>% 
      filter(month >= min_interval)
    
    for (split in 1:length(unique(data$split_ids))) {
      
      if(parametric_model) {
        
        model_formula_RHS <- paste(
          names(data_k)[which(!names(data_k) %in% c('id','has_med','dispensed','initiated_med',
                                                    'following_strategy',
                                                    'censored','outcome','split_ids'))], 
          collapse=" + "
        )
        
        if (length(unique(data$split_ids)) == 1) {
          
          model <- 
            glm(formula = as.formula(paste("has_med", model_formula_RHS, sep=" ~ ")),
                data = data_k %>% filter(initiated_med == initiated_treatment & 
                                           following_strategy == 1),
                family = quasibinomial())
          
        } else {
          
          model <- 
            glm(formula = as.formula(paste("has_med", model_formula_RHS, sep=" ~ ")),
                data = data_k %>% filter(split_ids != split & 
                                           initiated_med == initiated_treatment & 
                                           following_strategy == 1),
                family = quasibinomial())
          
        }
        
        prob_trt_k <- predict(model, newdata = data_k %>% filter(split_ids == split), type='response')
        
      } else {
        
        params_trt <- find_params_boosted_tree_model(covariates_df = data_k %>% 
                                                       filter(split_ids != split &
                                                                initiated_med == initiated_treatment & 
                                                                following_strategy == 1) %>% 
                                                       dplyr::select(-c('id','has_med','dispensed','initiated_med',
                                                                        'following_strategy',
                                                                        'censored','outcome','split_ids')), 
                                                     label_vector = data_k %>% 
                                                       filter(split_ids != split &
                                                                initiated_med == initiated_treatment & 
                                                                following_strategy == 1) %>% 
                                                       dplyr::select(c(has_med)) %>% {.[[1]]},
                                                     nfold = nfold,
                                                     tree_depth = tree_depth,
                                                     shrinkage_factor = shrinkage_factor,
                                                     num_trees = num_trees,
                                                     family = 'bernoulli',
                                                     num_cores = num_cores)
        
        prob_trt_k <- estimate_prob_boosted_tree_model(covariates_df = data_k %>% 
                                                         filter(split_ids != split &
                                                                  initiated_med == initiated_treatment & 
                                                                  following_strategy == 1) %>% 
                                                         dplyr::select(-c('id','has_med','dispensed','initiated_med',
                                                                          'following_strategy',
                                                                          'censored','outcome','split_ids')), 
                                                       label_vector = data_k %>% 
                                                         filter(split_ids != split &
                                                                  initiated_med == initiated_treatment & 
                                                                  following_strategy == 1) %>% 
                                                         dplyr::select(c(has_med)) %>% {.[[1]]},
                                                       params = params_trt,
                                                       boosted_tree_family = 'bernoulli',
                                                       predict_data = data_k %>% 
                                                         filter(split_ids == split) %>% 
                                                         dplyr::select(-c('id','has_med','dispensed','initiated_med',
                                                                          'following_strategy',
                                                                          'censored','outcome','split_ids')))
        
      }
      
      data <- data %>% 
        mutate(model_prob = case_when(month >= min_interval & split_ids == split ~ 1,
                                      TRUE ~ 0))
      
      data[which(data$model_prob == 1), ]$prob_trt <- prob_trt_k
      
      data <- data %>% 
        dplyr::select(-model_prob)
      
    }
    
  } else {
    
    for (k in min_interval:max_interval) {
      
      data_k <- create_covariates_df(covariates_df = data %>% dplyr::select(-prob_trt), 
                                     num_lags = min(c(num_lags,k-min_interval)), 
                                     vars_to_lag = vars_to_lag)
      
      subset_data_k <- data_k %>% filter(month == k)
      
      for (split in 1:length(unique(data$split_ids))) {
        
        if(parametric_model) {
          
          model_formula_RHS <- paste(
            names(subset_data_k)[which(!names(subset_data_k) %in% c('id','month','has_med','dispensed','initiated_med',
                                                                    'following_strategy',
                                                                    'censored','outcome','split_ids'))], 
            collapse=" + "
          )
          
          if (length(unique(data$split_ids)) == 1) {
            
            model_k <- 
              glm(formula = as.formula(paste("has_med", model_formula_RHS, sep=" ~ ")),
                  data = subset_data_k %>% filter(initiated_med == initiated_treatment & 
                                                    following_strategy == 1),
                  family = quasibinomial())
            
          } else {
            
            model_k <- 
              glm(formula = as.formula(paste("has_med", model_formula_RHS, sep=" ~ ")),
                  data = subset_data_k %>% filter(split_ids != split &
                                                    initiated_med == initiated_treatment & 
                                                    following_strategy == 1),
                  family = quasibinomial())
            
          }
          
          prob_trt_k <- predict(model_k, newdata = subset_data_k %>% filter(split_ids == split), type='response')
          
        } else {
          
          params_trt_k <- find_params_boosted_tree_model(covariates_df = subset_data_k %>% 
                                                           filter(split_ids != split &
                                                                    initiated_med == initiated_treatment & 
                                                                    following_strategy == 1) %>% 
                                                           dplyr::select(-c('id','month','has_med','dispensed','initiated_med',
                                                                            'following_strategy',
                                                                            'censored','outcome','split_ids')), 
                                                         label_vector = subset_data_k %>% 
                                                           filter(split_ids != split &
                                                                    initiated_med == initiated_treatment & 
                                                                    following_strategy == 1) %>% 
                                                           dplyr::select(c(has_med)) %>% {.[[1]]},
                                                         nfold = nfold,
                                                         tree_depth = tree_depth,
                                                         shrinkage_factor = shrinkage_factor,
                                                         num_trees = num_trees,
                                                         family = 'bernoulli',
                                                         num_cores = num_cores)
          
          prob_trt_k <- estimate_prob_boosted_tree_model(covariates_df = subset_data_k %>% 
                                                           filter(split_ids != split &
                                                                    initiated_med == initiated_treatment & 
                                                                    following_strategy == 1) %>% 
                                                           dplyr::select(-c('id','month','has_med','dispensed','initiated_med',
                                                                            'following_strategy',
                                                                            'censored','outcome','split_ids')), 
                                                         label_vector = subset_data_k %>% 
                                                           filter(split_ids != split &
                                                                    initiated_med == initiated_treatment & 
                                                                    following_strategy == 1) %>% 
                                                           dplyr::select(c(has_med)) %>% {.[[1]]},
                                                         params = params_trt_k,
                                                         boosted_tree_family = 'bernoulli',
                                                         predict_data = subset_data_k %>% 
                                                           filter(split_ids == split) %>% 
                                                           dplyr::select(-c('id','month','has_med','dispensed','initiated_med',
                                                                            'following_strategy',
                                                                            'censored','outcome','split_ids')))
          
        }
        
        data <- data %>% 
          mutate(model_prob = case_when(month == k & month >= min_interval & split_ids == split ~ 1,
                                        TRUE ~ 0))
          
          data[which(data$model_prob == 1), ]$prob_trt <- prob_trt_k
        
        data <- data %>% 
          dplyr::select(-model_prob)
        
      }
      
    }
    
  }
  
  
  if (is.null(f_int) & is.null(f_cond_int)) {
    
    data <- data %>%
      group_by(id) %>% 
      mutate(num_wt = case_when(num_intervals_without_trt == grace_period_length & has_med == 0 ~ 0,
                                TRUE ~ 1),
             denom_wt = case_when(num_intervals_without_trt == grace_period_length ~ prob_trt,
                                  TRUE ~ 1),
             Wt = num_wt/denom_wt) %>%
      ungroup()
    
  } else {
    
    if (is.null(f_cond_int)) {
      
      f_int_cond_adh <- f_int / head(c(1, 1-cumsum(f_int)), -1)
      
    } else {
      
      f_int_cond_adh <- f_cond_int
      
    }
    
    data <- data %>% 
      group_by(id) %>% 
      mutate(
        num_wt = case_when(num_intervals_without_trt < grace_period_length & has_med == 1 ~ 
                             f_int_cond_adh[num_intervals_without_trt+1],
                           num_intervals_without_trt < grace_period_length & has_med == 0 ~ 
                             1 - f_int_cond_adh[num_intervals_without_trt+1],
                           num_intervals_without_trt == grace_period_length & has_med == 1 ~ 1,
                           num_intervals_without_trt == grace_period_length & has_med == 0 ~ 0,
                           num_intervals_without_trt > grace_period_length ~ 0),
        denom_wt = case_when(has_med == 1 & dispensed == 1 ~ prob_trt,
                             has_med == 1 & dispensed == 0 ~ 1,
                             has_med == 0 ~ 1 - prob_trt),
        Wt = num_wt/denom_wt) %>%
      ungroup()
    
  }
  
  if (skip_first_interval) {data[which(data$month == min_interval-1), ]$Wt <- 1}
  
  return(data %>% dplyr::select(-c(num_wt, denom_wt)) %>% rename(Wt_adherence = Wt, prob_adherence = prob_trt))
  
}

################################################################################
# Computes iterative outcome models for either the iterative conditional expectation g-formula estimator or doubly-robust AIPW estimator.
# skip_first_interval: whether the first interval of follow-up should be included in the grace period intervention (T/F). Generally, all individuals will receive treatment in the first interval in the observed data so this should be set to FALSE.
# min_interval: the index of the first time-interval of the study (the first month of treatment initiation in this data example)
# max_interval: the investigator specified final-interval of intervention
# initiated_treatment: the medication to be initiated at baseline in the treatment strategy of interest
# f_int: a vector of probabilities defining the marginal probability of adherence in each grace period interval, should have length equal to grace period length
# f_cond_int: a vector of probabilities defining the conditional (on number of intervals without treatment) probability of adherence in each grace period interval, should have length equal to grace period length
# grace_period_length: the investigator specified grace period length. When the number of intervals without treatment equals this value, the conditional probability of treatment will be set to 1.

## To specify a natural grace period strategy, set f_int and f_cond_int to NULL. 
## To specify a stochastic grace period strategy, set values for either f_int or f_cond_int. f_cond_int will override f_int.
## ex: to define a stochastic grace period strategy with a uniform distribution of adherence with a grace period length of 3, set f_int = c(0.25, 0.25, 0.25) or f_cond_int = c(0.25, 0.33, 0.5), and grace_period_length = 3.

# model_formula_vars: a vector of the names of confounders to be included in the model for the outcome
# vars_to_lag: as a subset of models_formula_vars, the names of the variables to include histories of in the model for the outcome
# num_lags: the number of lagged histories for the variables in 'vars_to_lag' to include in the model for the outcome

# ex: to include the values of 'T2D' for the past 3 intervals including the current interval (the default) then set vars_to_lag = c('T2D') and num_lags = 3

# parametric_model: should a logistic GLM be used to compute the conditional outcome risk (T/F)
# data: the dataset
# nfold: if sample splitting is to be used with a non-parametric model for the outcome, the number of folds to use for cross-validation
# tree_depth: if sample splitting is to be used with a non-parametric model for the outcome, a vector of values for boosted tree depth
# shrinkage_factor: if sample splitting is to be used with a non-parametric model for the outcome, a vector of shrinkage values
# num_trees: if sample splitting is to be used with a non-parametric model for the outcome, a vector of values for number of trees (iterations)
# num_cores: if sample splitting is to be used with a non-parametric model for the outcome, number of cores to use for parallelization
################################################################################

iterative_outcomefnc <- function(skip_first_interval,
                                 max_interval,
                                 min_interval,
                                 initiated_treatment,
                                 f_int = NULL,
                                 f_cond_int = NULL,
                                 grace_period_length,
                                 model_formula_vars,
                                 vars_to_lag,
                                 num_lags,
                                 parametric_model,
                                 data,
                                 nfold,
                                 tree_depth,
                                 shrinkage_factor,
                                 num_trees,
                                 num_cores) {
  
  data <- data %>% 
    dplyr::select(c('id','month','has_med','dispensed','initiated_med',
                    'num_intervals_without_trt', 'following_strategy',
                    'censored','outcome','split_ids',
                    all_of(model_formula_vars))) %>%
    mutate(h_k_final = 0,
           outcome_resid = 0) %>%
    group_by(id) %>% 
    mutate(individual_max_month = max(month)) %>% 
    ungroup()
  
  for (split in 1:length(unique(data$split_ids))) {
    
    data <- data %>% 
      mutate(h_k_1 = outcome,
             Eh_k_1 = 0,
             h_k = 0)
    
    for (k in max_interval:min_interval) {
      
      data_k <- create_covariates_df(covariates_df = data %>% dplyr::select(-c(h_k, h_k_final, Eh_k_1, outcome, outcome_resid, individual_max_month)), 
                                     num_lags = min(c(num_lags,k-min_interval)), 
                                     vars_to_lag = vars_to_lag)
      
      subset_data_k <- filter(data_k, month == k & 
                                following_strategy == 1)
      subset_data_k_following_strategy <- filter(subset_data_k, 
                                                 (num_intervals_without_trt < grace_period_length | has_med == 1) &
                                                   initiated_med == initiated_treatment &
                                                   censored == 0)
      
      if (skip_first_interval & k == min_interval) {
        
        subset_data_k <- subset_data_k %>% dplyr::select(-c(has_med, num_intervals_without_trt))
        subset_data_k_following_strategy <- subset_data_k_following_strategy %>% dplyr::select(-c(has_med, num_intervals_without_trt))
        
      }
      
      if(parametric_model) {
        
        model_formula_RHS <- paste(
          names(subset_data_k)[which(!names(subset_data_k) %in% c('id','month','dispensed','initiated_med',
                                                                  'following_strategy',
                                                                  'censored','h_k_1','split_ids'))], 
          collapse=" + "
        )
        
        model_k <- 
          glm(formula = as.formula(paste("h_k_1", model_formula_RHS, sep=" ~ ")),
              data = subset_data_k_following_strategy %>% filter(split_ids != split),
              family = quasibinomial())
        
        prob_outcome_k <- predict(model_k, newdata = subset_data_k, type='response')
        
        if (!(skip_first_interval & k == min_interval)){
          
          prob_outcome_1_k <- predict(model_k, newdata = subset_data_k %>% mutate(has_med=1), type='response')
          prob_outcome_0_k <- predict(model_k, newdata = subset_data_k %>% mutate(has_med=0), type='response')
          prob_outcome_N_k <- predict(model_k, newdata = subset_data_k %>% 
                                        mutate(has_med = case_when(num_intervals_without_trt == grace_period_length ~ 1,
                                                                   TRUE ~ has_med)), type='response')
          
        }
        
      } else {
        
        params_outcome_k <- find_params_boosted_tree_model(covariates_df = subset_data_k_following_strategy %>% 
                                                             filter(split_ids != split) %>% 
                                                             dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                              'following_strategy',
                                                                              'censored','h_k_1','split_ids')), 
                                                           label_vector = subset_data_k_following_strategy %>% 
                                                             filter(split_ids != split) %>% 
                                                             dplyr::select(c(h_k_1)) %>% {.[[1]]},
                                                           nfold = nfold,
                                                           tree_depth = tree_depth,
                                                           shrinkage_factor = shrinkage_factor,
                                                           num_trees = num_trees,
                                                           family = 'gaussian',
                                                           num_cores = num_cores)
        
        prob_outcome_k <- estimate_prob_boosted_tree_model(covariates_df = subset_data_k_following_strategy %>%
                                                             filter(split_ids != split) %>% 
                                                             dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                              'following_strategy',
                                                                              'censored','h_k_1','split_ids')), 
                                                           label_vector = subset_data_k_following_strategy %>% 
                                                             filter(split_ids != split) %>% 
                                                             dplyr::select(c(h_k_1)) %>% {.[[1]]},
                                                           params = params_outcome_k,
                                                           boosted_tree_family = 'gaussian',
                                                           predict_data = subset_data_k %>%
                                                             dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                              'following_strategy',
                                                                              'censored','h_k_1','split_ids')))
        
        if (!(skip_first_interval & k == min_interval)){
          
          prob_outcome_1_k <- estimate_prob_boosted_tree_model(covariates_df = subset_data_k_following_strategy %>%
                                                                 filter(split_ids != split) %>% 
                                                                 dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                                  'following_strategy',
                                                                                  'censored','h_k_1','split_ids')), 
                                                               label_vector = subset_data_k_following_strategy %>% 
                                                                 filter(split_ids != split) %>% 
                                                                 dplyr::select(c(h_k_1)) %>% {.[[1]]},
                                                               params = params_outcome_k,
                                                               boosted_tree_family = 'gaussian',
                                                               predict_data = subset_data_k %>% mutate(has_med=1) %>%
                                                                 dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                                  'following_strategy',
                                                                                  'censored','h_k_1','split_ids')))
          
          prob_outcome_0_k <- estimate_prob_boosted_tree_model(covariates_df = subset_data_k_following_strategy %>%
                                                                 filter(split_ids != split) %>% 
                                                                 dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                                  'following_strategy',
                                                                                  'censored','h_k_1','split_ids')), 
                                                               label_vector = subset_data_k_following_strategy %>% 
                                                                 filter(split_ids != split) %>% 
                                                                 dplyr::select(c(h_k_1)) %>% {.[[1]]},
                                                               params = params_outcome_k,
                                                               boosted_tree_family = 'gaussian',
                                                               predict_data = subset_data_k %>% mutate(has_med=0) %>%
                                                                 dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                                  'following_strategy',
                                                                                  'censored','h_k_1','split_ids')))
          
          prob_outcome_N_k <- estimate_prob_boosted_tree_model(covariates_df = subset_data_k_following_strategy %>%
                                                                 filter(split_ids != split) %>% 
                                                                 dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                                  'following_strategy',
                                                                                  'censored','h_k_1','split_ids')), 
                                                               label_vector = subset_data_k_following_strategy %>% 
                                                                 filter(split_ids != split) %>% 
                                                                 dplyr::select(c(h_k_1)) %>% {.[[1]]},
                                                               params = params_outcome_k,
                                                               boosted_tree_family = 'gaussian',
                                                               predict_data = subset_data_k %>% 
                                                                 mutate(has_med = case_when(num_intervals_without_trt == grace_period_length ~ 1,
                                                                                            TRUE ~ has_med)) %>%
                                                                 dplyr::select(-c('id','month','dispensed','initiated_med',
                                                                                  'following_strategy',
                                                                                  'censored','h_k_1','split_ids')))
          
        }
        
      }
      
      data[which(data$month == k & 
                   data$following_strategy == 1), ]$Eh_k_1 <- prob_outcome_k
      
      if (skip_first_interval & k == min_interval) {
        
        data[which(data$month == k & 
                     data$following_strategy == 1), ]$h_k <- prob_outcome_k
        
      } else {
        
        if (is.null(f_int) & is.null(f_cond_int)) {
          
          data[which(data$month == k & 
                       data$following_strategy == 1), ]$h_k <- prob_outcome_N_k
          
        } else {
          
          if (is.null(f_cond_int)) {
            
            f_int_cond_adh <- f_int / head(c(1, 1-cumsum(f_int)), -1)
            
          } else {
            
            f_int_cond_adh <- f_cond_int
            
          }
          
          f_int_cond_adh <- c(f_int_cond_adh, rep(1, (max(data$month)) - length(f_int_cond_adh)))
          
          data[which(data$month == k & 
                       data$following_strategy == 1), ]$h_k <- 
            f_int_cond_adh[data[which(data$month == k & 
                                        data$following_strategy == 1), ]$num_intervals_without_trt + 1]*prob_outcome_1_k +
            (1-f_int_cond_adh[data[which(data$month == k & 
                                           data$following_strategy == 1), ]$num_intervals_without_trt + 1])*prob_outcome_0_k
          
        }
        
      }
      
      if (k != min_interval) {
        data <- data %>% 
          mutate(
            h_k_1 = case_when(month == k-1 & h_k_1 == 0 & month != individual_max_month ~ lead(h_k,1),
                              TRUE ~ as.numeric(h_k_1))
          )
      }
      
    }
    
    data[which(data$split_ids == split), ]$outcome_resid <- 
      data[which(data$split_ids == split), ]$h_k_1 - data[which(data$split_ids == split), ]$Eh_k_1
    
    data[which(data$split_ids == split), ]$h_k_final <- data[which(data$split_ids == split), ]$h_k
    
  }
  
  return(data %>% dplyr::select(-c(individual_max_month, h_k, h_k_1, Eh_k_1)) %>% rename(h_k = h_k_final))
  
}

################################################################################
# Computes the outcome risk under a medication initiation, loss to follow-up abolished, and natural or stochastic grace period adherence strategy
# initiated_treatment: the medication to be initiated at baseline in the treatment strategy of interest
# data: the dataset in person-month format. Should contain the following variables,
## 'id': a unique identifier for each individual in the study,
## 'month': a value for the interval number (month),
## 'has_med': an indicator of whether the individual has medication in the current interval based on their last dispensation and the amount of dispensed medication (1 if they have medication, 0 otherwise),
## 'dispensed': an indicator of whether the inidividual had a medication dispensation in the current interval (1 if yes, 0 if no),
## 'initiated_med': the name of the medication initiated at baseline,
## 'num_intervals_without_trt': an integer specifying the number of consecutive past intervals with has_med==0,
## 'censored': an indicator taking the value 1 (0 otherwise) at the interval at which an individual has been defined as lost to follow-up, can also be defined as the interval at which death occurs if death is to be treated as a censoring event,
## 'outcome': an indicator taking the value 1 (0 otherwise) at the interval at which the outcome has occurred,
## time-fixed and time-varying covariates

# pooled_time_treatment_model: whether separate models should be fit at each time interval for the probability of adherence, or whether the model should pool over time intervals (TRUE to pool over time intervals, FALSE for separate models)
# skip_first_interval: whether the first interval of follow-up should be included in the grace period intervention (T/F). Generally, all individuals will receive treatment in the first interval in the observed data so this should be set to FALSE.
# max_interval: the investigator specified final-interval of intervention
# initiated_treatment: the medication to be initiated at baseline in the treatment strategy of interest
# f_int: a vector of probabilities defining the marginal probability of adherence in each grace period interval, should have length equal to grace period length
# f_cond_int: a vector of probabilities defining the conditional (on number of intervals without treatment) probability of adherence in each grace period interval, should have length equal to grace period length
# grace_period_length: the investigator specified grace period length. When the number of intervals without treatment equals this value, the conditional probability of treatment will be set to 1.

## To specify a natural grace period strategy, set f_int and f_cond_int to NULL. 
## To specify a stochastic grace period strategy, set values for either f_int or f_cond_int. f_cond_int will override f_int.
## ex: to define a stochastic grace period strategy with a uniform distribution of adherence with a grace period length of 3, set f_int = c(0.25, 0.25, 0.25) or f_cond_int = c(0.25, 0.33, 0.5), and grace_period_length = 3.

# vars_to_lag_trt: if pooled_time_treatment_model == F, then vars_to_lag_trt can be used to specify, as a subset of models_formula_vars_trt, the names of the variables to include histories of in the model for treatment adherence
# num_lags_trt: if pooled_time_treatment_model == F, then num_lags_trt specifies the number of lagged histories for the variables in 'vars_to_lag_trt' to include in the model for treatment adherence
# model_formula_vars_trt: a vector of the names of confounders to be included in the model for treatment initiation and adherence
# parametric_model_trt: should a logistic GLM be used to compute probabilities of treatment initiation and adherence (T/F)
# nfold_trt: if sample splitting is to be used with a non-parametric model for probability of treatment initiation and adherence, the number of folds to use for cross-validation
# tree_depth_trt: if sample splitting is to be used with a non-parametric model for probability of treatment initiation and adherence, a vector of values for boosted tree depth
# shrinkage_factor_trt: if sample splitting is to be used with a non-parametric model for probability of treatment initiation and adherence, a vector of shrinkage values
# num_trees_trt: if sample splitting is to be used with a non-parametric model for probability of treatment initiation and adherence, a vector of values for number of trees (iterations)
# num_cores_trt: if sample splitting is to be used with a non-parametric model for probability of treatment initiation adherence, number of cores to use for parallelization
# parametric_outcome_non_DR: if IPW (rather than a doubly robust estimator) is to be used, whether a marginal structural model for the outcome should be fit (T/F)
# model_formula_vars_censoring: a vector of the names of confounders to be included in the model for censoring
# parametric_model_censoring: should a logistic GLM be used to compute probabilities of censoring (T/F)
# nfold_censoring: if sample splitting is to be used with a non-parametric model for probability of censoring, the number of folds to use for cross-validation
# tree_depth_censoring: if sample splitting is to be used with a non-parametric model for probability of censoring, a vector of values for boosted tree depth
# shrinkage_factor_censoring: if sample splitting is to be used with a non-parametric model for probability of censoring, a vector of shrinkage values
# num_trees_censoring: if sample splitting is to be used with a non-parametric model for probability of censoring, a vector of values for number of trees (iterations)
# num_cores_censoring: if sample splitting is to be used with a non-parametric model for probability of censoring, number of cores to use for parallelization
# trim_weights_quantile: to trim extreme weights, choose a numeric value > 0 and < 1. For no trimming, choose 1.
# vars_to_lag_outcome: as a subset of models_formula_vars_outcome, the names of the variables to include histories of in the model for the outcome
# num_lags_outcome: the number of lagged histories for the variables in 'vars_to_lag_outcome' to include in the model for the outcome
# model_formula_vars_outcome: a vector of the names of confounders to be included in the model for the outcome
# parametric_model_outcome: should a logistic GLM be used to compute the conditional outcome risk (T/F)
# nfold_outcome: if sample splitting is to be used with a non-parametric model for the outcome, the number of folds to use for cross-validation
# tree_depth_outcome: if sample splitting is to be used with a non-parametric model for the outcome, a vector of values for boosted tree depth
# shrinkage_factor_outcome: if sample splitting is to be used with a non-parametric model for the outcome, a vector of shrinkage values
# num_trees_outcome: if sample splitting is to be used with a non-parametric model for the outcome, a vector of values for number of trees (iterations)
# num_cores_outcome: if sample splitting is to be used with a non-parametric model for the outcome, number of cores to use for parallelization
# nsplits: if sample splitting it to be used with a non-parametric model for treatment initiation, treatment adherence, censoring, and outcome, the number of sample splits
# DR: whether a doubly robust AIPW estimator should be used rather than IPW (T/F)
################################################################################

resultsfnc <- function(initiated_treatment,
                       data, 
                       pooled_time_treatment_model = F,
                       skip_first_interval = F,
                       max_interval = NULL, 
                       f_int = NULL,
                       f_cond_int = NULL,
                       grace_period_length,
                       vars_to_lag_trt = NULL,
                       num_lags_trt = NULL,
                       model_formula_vars_trt,
                       parametric_model_trt = F,
                       nfold_trt = NULL,
                       tree_depth_trt = NULL,
                       shrinkage_factor_trt = NULL,
                       num_trees_trt = NULL,
                       num_cores_trt = NULL,
                       parametric_outcome_non_DR = NULL,
                       model_formula_vars_censoring = NULL,
                       parametric_model_censoring = NULL,
                       nfold_censoring = NULL,
                       tree_depth_censoring = NULL,
                       shrinkage_factor_censoring = NULL,
                       num_trees_censoring = NULL,
                       num_cores_censoring = NULL,
                       trim_weights_quantile = 1,
                       vars_to_lag_outcome = NULL,
                       num_lags_outcome = NULL,
                       model_formula_vars_outcome = NULL,
                       parametric_model_outcome = NULL,
                       nfold_outcome = NULL,
                       tree_depth_outcome = NULL,
                       shrinkage_factor_outcome = NULL,
                       num_trees_outcome = NULL,
                       num_cores_outcome = NULL,
                       nsplits = NULL,
                       DR = T){
  
  min_interval <- min(data$month)
  
  if (is.null(max_interval)) {max_interval <- max(data$month)} else {
    data <- data %>% filter(month <= max_interval)
  }
  
  data <- create_following_strategy_var(initiated_treatment = initiated_treatment,
                                        grace_period_length = grace_period_length,
                                        min_interval = min_interval,
                                        data = data)
  
  if (!is.null(nsplits)) {
    split_ids <- sample(rep(1:nsplits,ceiling(length(unique(data$id))/nsplits))[1:length(unique(data$id))])
  } else {
    split_ids <- rep(1,length(unique(data$id)))
  }
  
  data <- left_join(data, data.frame(id = unique(data$id), split_ids))
  
  data_wts <- adherence_weightsfnc(skip_first_interval = skip_first_interval,
                                   f_int = f_int,
                                   f_cond_int = f_cond_int,
                                   grace_period_length = grace_period_length,
                                   min_interval = min_interval,
                                   max_interval = max_interval,
                                   initiated_treatment = initiated_treatment,
                                   pooled_time_treatment_model = pooled_time_treatment_model,
                                   model_formula_vars = model_formula_vars_trt,
                                   vars_to_lag = vars_to_lag_trt,
                                   num_lags = num_lags_trt,
                                   parametric_model = parametric_model_trt,
                                   data = data,
                                   nfold = nfold_trt,
                                   tree_depth = tree_depth_trt,
                                   shrinkage_factor = shrinkage_factor_trt,
                                   num_trees = num_trees_trt,
                                   num_cores = num_cores_trt)
  
  data_wts$Wt_initiation <- initiation_weightsfnc(min_interval = min_interval,
                                                  initiated_treatment = initiated_treatment,
                                                  model_formula_vars = model_formula_vars_trt,
                                                  parametric_model = parametric_model_trt,
                                                  data = data,
                                                  nfold = nfold_trt,
                                                  tree_depth = tree_depth_trt,
                                                  shrinkage_factor = shrinkage_factor_trt,
                                                  num_trees = num_trees_trt,
                                                  num_cores = num_cores_trt)
  
  data_wts$Wt <- data_wts$Wt_adherence*data_wts$Wt_initiation
  
  if (!is.null(model_formula_vars_censoring)) {
    
    data_wts$Wt_censoring <- censoring_weightsfnc(min_interval = min_interval,
                                                  initiated_treatment = initiated_treatment,
                                                  model_formula_vars = model_formula_vars_censoring,
                                                  parametric_model = parametric_model_censoring,
                                                  data = data,
                                                  nfold = nfold_censoring,
                                                  tree_depth = tree_depth_censoring,
                                                  shrinkage_factor = shrinkage_factor_censoring,
                                                  num_trees = num_trees_censoring,
                                                  num_cores = num_cores_censoring)
    
    data_wts$Wt <- data_wts$Wt*data_wts$Wt_censoring
    
  }
  
  data_wts <- data_wts %>% 
    group_by(id) %>%
    mutate(Wt = cumprod(Wt)) %>%
    ungroup() %>%
    mutate(Wt = case_when(Wt > quantile(Wt, trim_weights_quantile) ~ quantile(Wt, trim_weights_quantile), TRUE ~ Wt))
  
  if (!DR){
    
    if (parametric_outcome_non_DR) {
      
      standard_data <- data.frame(month = min_interval:max_interval)
      
      outcome_model <- glm(outcome ~ 
                             ns(month, knots = seq(min_interval+1,max_interval-1,2), Boundary.knots = c(min_interval,max_interval)), 
                           data = data_wts,
                           family = binomial(),
                           weights = weight
      )
      
      predicted <- standard_data %>%
        mutate(cond_surv = 1-predict(outcome_model, ., type="response")) %>% mutate(CI = 1-cumprod(cond_surv)) 
      
    } else {
      
      predicted <- data_wts %>% 
        group_by(month) %>% 
        summarise(cond_surv = 1- (sum(outcome*Wt) / sum(Wt))) %>%
        {data.frame(month = min_interval:max_interval, CI = 1-cumprod(.$cond_surv))}
      
    }
    
    return(list(data = data_wts,
                outcome = data.frame(CI=c(0,predicted[["CI"]])) %>% mutate(month = min_interval:(max_interval+1))))
    
  } else {
    
    data_outcome <- iterative_outcomefnc(skip_first_interval = skip_first_interval,
                                         max_interval = max_interval,
                                         min_interval = min_interval,
                                         initiated_treatment = initiated_treatment,
                                         f_int = f_int,
                                         f_cond_int = f_cond_int,
                                         grace_period_length = grace_period_length,
                                         model_formula_vars = model_formula_vars_outcome,
                                         vars_to_lag = vars_to_lag_outcome,
                                         num_lags = num_lags_outcome,
                                         parametric_model = parametric_model_outcome,
                                         data = data,
                                         nfold = nfold_outcome,
                                         tree_depth = tree_depth_outcome,
                                         shrinkage_factor = shrinkage_factor_outcome,
                                         num_trees = num_trees_outcome,
                                         num_cores = num_cores_outcome)
    
    data <- data_outcome %>% mutate(prob_adherence = data_wts$prob_adherence,
                                    Wt_adherence = data_wts$Wt_adherence,
                                    Wt_initiation = data_wts$Wt_initiation,
                                    Wt_censoring = data_wts$Wt_censoring,
                                    Wt = data_wts$Wt)
    
    data <- data %>% 
      mutate(D_k = Wt*outcome_resid)
    
    data <- data %>% 
      group_by(id) %>% 
      mutate(cumsum_Dk = cumsum(D_k)) %>% 
      ungroup()
    
    EIF_data <- data.frame(
      EIF = (data %>% group_by(id) %>% filter(month == min(max_interval, max(month))))$cumsum_Dk +
        (data %>% group_by(id) %>% filter(month == min_interval))$h_k,
      split_ids = data %>% filter(month==min_interval) %>% {.$split_ids}
    )
    
    DR_outcome <- 
      sapply(1:nsplits, function(split) {
        
        split_data <- EIF_data %>% filter(split_ids == split)
        
        mean(split_data$EIF)
        
      })
    
    DR_var <- 
      sapply(1:nsplits, function(split) {
        
        split_data <- EIF_data %>% filter(split_ids == split)
        
        sum((split_data$EIF - mean(split_data$EIF))^2) / (length(split_data$EIF)^2)
        
      })
    
    return(list(data = data,
                EIF_data = EIF_data,
                DR_outcome = DR_outcome,
                DR_var = DR_var))
    
  }
}

################################################################################
# Taking a 'results' list produced by 'resultsfnc' with DR=T, returns the doubly robust sample split estimate, its standard error, and confidence interval by taking the average of the estimates over the splits
# results: the results list produced by 'resultsfnc'
# CI_level: the nominal confidence interval coverage probability
################################################################################

combine_DR_splits <- function(results, CI_level) {
  
  outcome = mean(results$DR_outcome)
  
  se = sqrt((1/(length(results$DR_var)^2))*sum(results$DR_var))
  
  CI = c(outcome + abs(qnorm((1-CI_level)/2))*se, 
         outcome - abs(qnorm((1-CI_level)/2))*se)
  
  return(list(outcome = outcome,
              se = se,
              CI = CI))
  
}

################################################################################
# Taking a 'results' list produced by 'resultsfnc' with DR=T, returns the doubly robust sample split estimate, its standard error, and confidence interval
# results: the results list produced by 'resultsfnc'
# CI_level: the nominal confidence interval coverage probability
################################################################################

compute_risk <- function(results, CI_level) {
  
  outcome = mean(results$EIF_data$EIF)
  
  se <- sqrt(
    sum((results$EIF_data$EIF - outcome)^2) / (length(results$EIF_data$EIF)^2)
  )
  
  CI = c(outcome + abs(qnorm((1-CI_level)/2))*se, 
         outcome - abs(qnorm((1-CI_level)/2))*se)
  
  return(list(outcome = outcome,
              se = se,
              CI = CI))
  
}

################################################################################
# Taking a 'results' list produced by 'resultsfnc' with DR=T, returns the doubly robust sample split estimate, its standard error, and confidence interval after trimming the IP weights
# trim_weights_quantile: to trim extreme weights, choose a numeric value > 0 and < 1. For no trimming, choose 1
# results: the results list produced by 'resultsfnc'
# CI_level: the nominal confidence interval coverage probability
# max_interval: the investigator specified final-interval of intervention
################################################################################

trim_weights_compute_risk <- function(trim_weights_quantile, results, CI_level, max_interval=NULL) {
  
  min_interval <- min(results$data$month)
  
  if (is.null(max_interval)) {max_interval <- max(results$data$month)} else {
    results$data <- results$data %>% filter(month <= max_interval)
  }
  
  results$data <- results$data %>% 
    mutate(Wt = case_when(Wt > quantile(Wt, trim_weights_quantile) ~ quantile(Wt, trim_weights_quantile), TRUE ~ Wt))
  
  results$data <- results$data %>% 
    mutate(D_k = Wt*outcome_resid)
  
  results$data <- results$data %>% 
    group_by(id) %>% 
    mutate(cumsum_Dk = cumsum(D_k)) %>% 
    ungroup()
  
  results$EIF_data <- data.frame(
    EIF = (results$data %>% group_by(id) %>% filter(month == min(max_interval, max(month))))$cumsum_Dk +
      (results$data %>% group_by(id) %>% filter(month == min_interval))$h_k,
    split_ids = results$data %>% filter(month==min_interval) %>% {.$split_ids}
  )
  
  outcome = mean(results$EIF_data$EIF)
  
  se <- sqrt(
    sum((results$EIF_data$EIF - outcome)^2) / (length(results$EIF_data$EIF)^2)
  )
  
  CI = c(outcome + abs(qnorm((1-CI_level)/2))*se, 
         outcome - abs(qnorm((1-CI_level)/2))*se)
  
  return(list(outcome = outcome,
              se = se,
              CI = CI))
  
}

################################################################################
# Taking two 'results' lists produced by 'resultsfnc' with DR=T, returns the doubly robust sample split estimate of the risk difference, its standard error, and confidence interval
# results1: a results list produced by 'resultsfnc'
# results2: a second results list produced by 'resultsfnc'
# CI_level: the nominal confidence interval coverage probability
################################################################################

compute_risk_difference <- function(results1, results2, CI_level) {
  
  outcome <- mean(results1$EIF_data$EIF - results2$EIF_data$EIF)
  
  se <- sqrt(
    sum(((results1$EIF_data$EIF - results2$EIF_data$EIF) - outcome)^2) / (length(results1$EIF_data$EIF)^2)
  )
  
  CI = c(outcome + abs(qnorm((1-CI_level)/2))*se, 
         outcome - abs(qnorm((1-CI_level)/2))*se)
  
  return(list(outcome = outcome,
              se = se,
              CI = CI))
  
}

################################################################################
# Taking two 'results' lists produced by 'resultsfnc' with DR=T, returns the doubly robust sample split estimate of the risk difference, its standard error, and confidence interval after trimming the IP weights
# trim_weights_quantile: to trim extreme weights, choose a numeric value > 0 and < 1. For no trimming, choose 1
# results1: a results list produced by 'resultsfnc'
# results2: a second results list produced by 'resultsfnc'
# results: the results list produced by 'resultsfnc'
# CI_level: the nominal confidence interval coverage probability
# max_interval: the investigator specified final-interval of intervention
################################################################################

trim_weights_compute_risk_difference <- function(trim_weights_quantile, results1, results2, CI_level, max_interval=NULL) {
  
  min_interval <- min(results1$data$month)
  
  if (is.null(max_interval)) {max_interval <- max(results1$data$month)} else {
    results1$data <- results1$data %>% filter(month <= max_interval)
    results2$data <- results2$data %>% filter(month <= max_interval)
  }
  
  results1$data <- results1$data %>% 
    mutate(Wt = case_when(Wt > quantile(Wt, trim_weights_quantile) ~ quantile(Wt, trim_weights_quantile), TRUE ~ Wt))
  
  results2$data <- results2$data %>% 
    mutate(Wt = case_when(Wt > quantile(Wt, trim_weights_quantile) ~ quantile(Wt, trim_weights_quantile), TRUE ~ Wt))
  
  results1$data <- results1$data %>% 
    mutate(D_k = Wt*outcome_resid)
  
  results2$data <- results2$data %>% 
    mutate(D_k = Wt*outcome_resid)
  
  results1$data <- results1$data %>% 
    group_by(id) %>% 
    mutate(cumsum_Dk = cumsum(D_k)) %>% 
    ungroup()
  
  results2$data <- results2$data %>% 
    group_by(id) %>% 
    mutate(cumsum_Dk = cumsum(D_k)) %>% 
    ungroup()
  
  results1$EIF_data <- data.frame(
    EIF = (results1$data %>% group_by(id) %>% filter(month == min(max_interval, max(month))))$cumsum_Dk +
      (results1$data %>% group_by(id) %>% filter(month == min_interval))$h_k,
    split_ids = results1$data %>% filter(month==min_interval) %>% {.$split_ids}
  )
  
  results2$EIF_data <- data.frame(
    EIF = (results2$data %>% group_by(id) %>% filter(month == min(max_interval, max(month))))$cumsum_Dk +
      (results2$data %>% group_by(id) %>% filter(month == min_interval))$h_k,
    split_ids = results2$data %>% filter(month==min_interval) %>% {.$split_ids}
  )
  
  outcome <- mean(results1$EIF_data$EIF - results2$EIF_data$EIF)
  
  se <- sqrt(
    sum(((results1$EIF_data$EIF - results2$EIF_data$EIF) - outcome)^2) / (length(results1$EIF_data$EIF)^2)
  )
  
  CI = c(outcome + abs(qnorm((1-CI_level)/2))*se, 
         outcome - abs(qnorm((1-CI_level)/2))*se)
  
  return(list(outcome = outcome,
              se = se,
              CI = CI))
  
}

################################################################################
# Taking a dataset with weights, computes the IPW estimate of risk
# data_wts: dataset of the same format as specified in 'resultsfnc' with a weight variable 'Wt'
################################################################################

compute_IPW_outcome <- function(data_wts){
  
  predicted <- data_wts %>% 
    group_by(month) %>% 
    summarise(cond_surv = 1- (sum(outcome*Wt) / sum(Wt))) %>%
    {data.frame(month = min(data_wts$month):max(data_wts$month), CI = 1-cumprod(.$cond_surv))}
  
  
  return(list(data = data_wts,
              outcome = data.frame(CI=c(0,predicted[["CI"]])) %>% mutate(month = min(data_wts$month):(max(data_wts$month)+1))))
  
}
