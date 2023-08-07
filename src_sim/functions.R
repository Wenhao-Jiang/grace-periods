
################################################################################
# Creates the initial simulated dataset in person-time format, assigning treatment
# and creating placeholders for time-varying variables including adherence and death.
# Time-varying variables will be subsequently updated iteratively with the 'create_sim_data' function.
# params: a list with the following structure,
## list(
##   sample_size = the number of individuals to simulate,
##   num_intervals = the number of intervals for each individual,
##   p_Y = the probability of death in an interval in the absence of treatment adherence,
##   p_AY = the probability of death in an interval if adherent (in that interval) to medication A,
##   p_BY = the probability of death in an interval if adherent (in that interval) to medication B,
##   p_nonadherence_A = the probability of nonadherence in the current interval if adherent to medication A in the prior interval,
##   p_nonadherence_B = the probability of nonadherence in the current interval if adherent to medication B in the prior interval,
##   p_nonadherence_baseline_1 = the probability of nonadherence in the current interval if nonadherent in the prior interval,
##   p_nonadherence_baseline_2 = the probability of nonadherence in the current interval if nonadherent in the prior 2 consecutive intervals,
##   p_nonadherence_baseline_3 = the probability of nonadherence in the current interval if nonadherent in the prior 3 consecutive intervals,
##   p_nonadherence_baseline_4 = the probability of nonadherence in the current interval if nonadherent in the prior 4 consecutive intervals,
##   p_nonadherence_baseline_5 = the probability of nonadherence in the current interval if nonadherent in the prior 5 consecutive intervals,
##   p_nonadherence_baseline_6 = the probability of nonadherence in the current interval if nonadherent in the prior 6 consecutive intervals
## )
################################################################################

initialize_sim_data <- function(params){
  
  with(params, {
    
    trt_A <- as.numeric(c(rep(1, sample_size/2), rep(0, sample_size/2)))
    trt_B <- 1-trt_A
    
    out <- data.frame(
      id = rep(1:sample_size, each = num_intervals),
      trt_A = rep(trt_A, each = num_intervals),
      trt_B = rep(trt_B, each = num_intervals),
      interval = rep(1:num_intervals, times = sample_size),
      adherent = rep(0, times = sample_size*num_intervals),
      dead = rep(0, times = sample_size*num_intervals)
    )
    
    return(out)
    
  })
  
}

################################################################################
# Initializes the simulated data then iteratively updates the time-varying variables
# params: a list with the following structure,
## list(
##   sample_size = the number of individuals to simulate,
##   num_intervals = the number of intervals for each individual,
##   p_Y = the probability of death in an interval in the absence of treatment adherence,
##   p_AY = the probability of death in an interval if adherent (in that interval) to medication A,
##   p_BY = the probability of death in an interval if adherent (in that interval) to medication B,
##   p_nonadherence_A = the probability of nonadherence in the current interval if adherent to medication A in the prior interval,
##   p_nonadherence_B = the probability of nonadherence in the current interval if adherent to medication B in the prior interval,
##   p_nonadherence_baseline_1 = the probability of nonadherence in the current interval if nonadherent in the prior interval,
##   p_nonadherence_baseline_2 = the probability of nonadherence in the current interval if nonadherent in the prior 2 consecutive intervals,
##   p_nonadherence_baseline_3 = the probability of nonadherence in the current interval if nonadherent in the prior 3 consecutive intervals,
##   p_nonadherence_baseline_4 = the probability of nonadherence in the current interval if nonadherent in the prior 4 consecutive intervals,
##   p_nonadherence_baseline_5 = the probability of nonadherence in the current interval if nonadherent in the prior 5 consecutive intervals,
##   p_nonadherence_baseline_6 = the probability of nonadherence in the current interval if nonadherent in the prior 6 consecutive intervals
## )
################################################################################

create_sim_data <- function(params) {
  
  with(params, {
    
    sim_data <- initialize_sim_data(params)
    
    sim_data$adherent[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
      rbinom(n=sample_size, size=1,
             p = 1-p_nonadherence_baseline_1)
    
    sim_data$dead[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] <- 
      rbinom(n=sample_size, size=1,
             p = p_Y*(1-sim_data$adherent[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)]) +
               p_AY*sim_data$adherent[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)]*
                sim_data$trt_A[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)] +
               p_BY*sim_data$adherent[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)]*
                sim_data$trt_B[seq(1, 1+(num_intervals*(sample_size-1)), by=num_intervals)])
    
    for (int in 2:num_intervals) {
      
      sim_data <- sim_data %>% 
        group_by(id) %>% mutate(
          
        adherent = case_when(
          interval == int & lag(adherent, 1) == 1 & trt_A == 1 ~ 
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_A)),
          
          interval == int & lag(adherent, 1) == 1 & trt_B == 1 ~ 
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_B)),
          
          interval == int & lag(adherent, 1) == 0 & (lag(adherent, 2) == 1 | interval == 2) ~ 
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_baseline_1)),
          
          interval == int & lag(adherent, 1) == 0 & lag(adherent, 2) == 0 & (lag(adherent, 3) == 1 | interval == 3) ~ 
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_baseline_2)),
          
          interval == int & lag(adherent, 1) == 0 & lag(adherent, 2) == 0 & lag(adherent, 3) == 0 & (lag(adherent, 4) == 1 | interval == 4) ~ 
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_baseline_3)),
          
          interval == int & lag(adherent, 1) == 0 & lag(adherent, 2) == 0 & lag(adherent, 3) == 0 & lag(adherent, 4) == 0 &
            (lag(adherent, 5) == 1 | interval == 5) ~ 
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_baseline_4)),
          
          interval == int & lag(adherent, 1) == 0 & lag(adherent, 2) == 0 & lag(adherent, 3) == 0 & lag(adherent, 4) == 0 & lag(adherent, 5) == 0 &
            (lag(adherent, 6) == 1 | interval == 6) ~ 
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_baseline_5)),
          
          interval == int & lag(adherent, 1) == 0 & lag(adherent, 2) == 0 & lag(adherent, 3) == 0 & lag(adherent, 4) == 0 & lag(adherent, 5) == 0 &
            lag(adherent, 6) == 0 & (lag(adherent, 7) == 1 | interval == 7) ~  
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_baseline_6)),
          
          interval != int ~ 
            adherent
        )
      ) %>% ungroup()
      
      sim_data <- sim_data %>% 
        group_by(id) %>% mutate(
        dead = case_when(
          interval == int & lag(dead, 1) == 1 ~ as.numeric(1),
          
          interval == int & lag(dead, 1) == 0 & adherent == 1 & trt_A == 1 ~ 
            as.numeric(rbinom(size=1, n=n(), p = p_AY)),
          
          interval == int & lag(dead, 1) == 0 & adherent == 1 & trt_B == 1 ~ 
            as.numeric(rbinom(size=1, n=n(), p = p_BY)),
          
          interval == int & lag(dead, 1) == 0 & adherent == 0 ~ 
            as.numeric(rbinom(size=1, n=n(), p = p_Y)),
          
          interval != int ~ 
            dead
        )
      ) %>% ungroup()
      
    }
    
    sim_data <- sim_data %>% 
      group_by(id) %>% 
      filter(dead==0 | lag(dead==0) | interval==1) %>% 
      ungroup()
    
    sim_data$temp <- rep(1:(length(rle(sim_data$adherent)$length)), rle(sim_data$adherent)$length)
    
    sim_data <- sim_data %>% 
      group_by(id, temp) %>%
      mutate(num_intervals_nonadherent = row_number()*(1-adherent)) %>%
      group_by(id) %>%
      mutate(num_intervals_without_trt = lag(num_intervals_nonadherent,1,default=0)) %>% 
      ungroup() %>% 
      dplyr::select(-c(temp, num_intervals_nonadherent))
    
    return(sim_data)
    
  })
  
}

################################################################################
# Computes the outcome risks under a stochastic or natural grace period adherence strategy
# with initiation of medication A and medication B separately. Uses IPW with a logistic GLM.
# data: the simulated dataset
# wt_denom_formula: a formula object defining the logistic GLM formula for the inverse probability of treatment (adherence) weights
# grace_period_length: the investigator specified grace period length. When the number of intervals without treatment equals this value, the conditional probability of treatment will be set to 1.
# f_int: a vector of probabilities defining the marginal probability of adherence in each grace period interval, should have length equal to grace period length
# f_cond_int: a vector of probabilities defining the conditional (on number of intervals without treatment) probability of adherence in each grace period interval, should have length equal to grace period length

## To specify a natural grace period strategy, set f_int and f_cond_int to NULL. 
## To specify a stochastic grace period strategy, set values for either f_int or f_cond_int. f_cond_int will override f_int.
## ex: to define a stochastic grace period strategy with a uniform distribution of adherence with a grace period length of 3, set f_int = c(0.25, 0.25, 0.25) or f_cond_int = c(0.25, 0.33, 0.5), and grace_period_length = 3.

################################################################################

resultsfnc <- function(data, 
                       wt_denom_formula, 
                       grace_period_length = NULL, 
                       f_int = NULL, 
                       f_cond_int = NULL){
  
  data_wts <- data %>% mutate(Wt = estimate_weights(
    data = data, 
    wt_denom_formula = wt_denom_formula, 
    grace_period_length = grace_period_length,
    f_int = f_int,
    f_cond_int = f_cond_int
  ))
  
  estimate <- data.frame(interval = rep(1:max(data_wts$interval), each=2),
             trt_A = rep(0:1, max(data_wts$interval)),
             cond_surv = 1 - (
               (data_wts %>% 
                  group_by(interval, trt_A) %>% 
                  summarise(num_deaths = sum(Wt*dead)) %>% {.$num_deaths}) / 
               (data_wts %>% 
                  group_by(interval, trt_A) %>% 
                  summarise(sum_wts = sum(Wt)) %>% {.$sum_wts}))) %>% 
    group_by(trt_A) %>%
    mutate(CI = 1-cumprod(cond_surv))
  
  return(estimate)
}

################################################################################
# Returns inverse probability weights for adherence according to an investigator specified grace period strategy.
# Can be used to construct weights for the natural or stochastic grace period strategies.
# data: the simulated dataset
# wt_denom_formula: a formula object defining the logistic GLM formula for the inverse probability of treatment (adherence) weights
# grace_period_length: the investigator specified grace period length. When the number of intervals without treatment equals this value, the conditional probability of treatment will be set to 1.
# f_int: a vector of probabilities defining the marginal probability of adherence in each grace period interval, should have length equal to grace period length
# f_cond_int: a vector of probabilities defining the conditional (on number of intervals without treatment) probability of adherence in each grace period interval, should have length equal to grace period length

## To specify a natural grace period strategy, set f_int and f_cond_int to NULL. 
## To specify a stochastic grace period strategy, set values for either f_int or f_cond_int. f_cond_int will override f_int.
## ex: to define a stochastic grace period strategy with a uniform distribution of adherence with a grace period length of 3, set f_int = c(0.25, 0.25, 0.25) or f_cond_int = c(0.25, 0.33, 0.5), and grace_period_length = 3.

################################################################################

estimate_weights <- function(data, 
                             wt_denom_formula, 
                             grace_period_length = NULL, 
                             f_int = NULL, 
                             f_cond_int = NULL){
  
  data$prob_adherence <- 
    predict(speedglm(formula = wt_denom_formula, 
               data = data, 
               family = gaussian()), 
           newdata = data, type="response")
  
  if (is.null(f_int) & is.null(f_cond_int)) {
    
    data_with_weights <- data %>%
      group_by(id) %>% 
      mutate(denom_wt = case_when(num_intervals_without_trt == grace_period_length ~ prob_adherence,
                                  TRUE ~ 1),
             num_wt = case_when(num_intervals_without_trt == grace_period_length & adherent == 0 ~ 0,
                                TRUE ~ 1),
             Wt = num_wt/denom_wt) %>%
      mutate(Wt_cumprod = cumprod(Wt)) %>%
      ungroup()
    
  } else {
    
    if (is.null(f_cond_int)) {
      
      f_int_cond_adh <- f_int / head(c(1, 1-cumsum(f_int)), -1)
      
    } else {
      
      f_int_cond_adh <- f_cond_int
      
    }
    
    f_int_cond_adh <- c(f_int_cond_adh, rep(1, (max(data$num_intervals_without_trt+1)) - length(f_int_cond_adh)))
           
    f_int_cond_nonadh <- 1-f_int_cond_adh
    
    data$denom_wt <- data$adherent*data$prob_adherence + (1-data$adherent)*(1-data$prob_adherence)
    data$num_wt <- data$adherent*f_int_cond_adh[data$num_intervals_without_trt+1] + 
      (1-data$adherent)*f_int_cond_nonadh[data$num_intervals_without_trt+1]
    data$Wt <- data$num_wt/data$denom_wt
    
    data_with_weights <- data %>%
      group_by(id) %>% 
      mutate(Wt_cumprod = cumprod(Wt)) %>%
      ungroup()
    
  }

  return(data_with_weights$Wt_cumprod)
}

################################################################################
# Resamples (with replacement) individuals from a provided dataset and 
# computes the outcome risks under a stochastic or natural grace period adherence strategy
# with initiation of medication A and medication B separately. Uses IPW with a logistic GLM.
# Can be replicated to compute a bootstrap estimator.
# data: the simulated dataset
# wt_denom_formula: a formula object defining the logistic GLM formula for the inverse probability of treatment (adherence) weights
# grace_period_length: the investigator specified grace period length. When the number of intervals without treatment equals this value, the conditional probability of treatment will be set to 1.
# f_int: a vector of probabilities defining the marginal probability of adherence in each grace period interval, should have length equal to grace period length
# f_cond_int: a vector of probabilities defining the conditional (on number of intervals without treatment) probability of adherence in each grace period interval, should have length equal to grace period length

## To specify a natural grace period strategy, set f_int and f_cond_int to NULL. 
## To specify a stochastic grace period strategy, set values for either f_int or f_cond_int. f_cond_int will override f_int.
## ex: to define a stochastic grace period strategy with a uniform distribution of adherence with a grace period length of 3, set f_int = c(0.25, 0.25, 0.25) or f_cond_int = c(0.25, 0.33, 0.5), and grace_period_length = 3.

################################################################################

bootstrapfnc <- function(data, 
                         wt_denom_formula, 
                         grace_period_length = NULL, 
                         f_int = NULL, 
                         f_cond_int = NULL){
  
  IDs <- sample(unique(data$id), length(unique(data$id)), replace = TRUE)
  datatemp <- as.data.table(data)
  setkey(datatemp, "id")
  data_resample <- datatemp[J(IDs), allow.cartesian = TRUE]
  data_resample <- data_resample %>% mutate(id = case_when(interval==1 ~ 1, TRUE ~ 0)) %>% mutate(id = cumsum(id))
  
  return(
    resultsfnc(
      data = data_resample,
      wt_denom_formula = wt_denom_formula,
      grace_period_length = grace_period_length,
      f_int = f_int,
      f_cond_int = f_cond_int
    )
  )
}
