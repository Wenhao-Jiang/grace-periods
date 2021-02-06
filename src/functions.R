initialize_sim_data <- function(params){
  
  with(params, {
    
    trt_A <- as.numeric(rbernoulli(sample_size, p=0.5))
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
          
          interval == int & lag(adherent, 1) == 0 & lag(adherent, 2) == 0 & lag(adherent, 3) == 0 ~ 
            as.numeric(rbinom(size=1, n=n(), p = 1 - p_nonadherence_baseline_3)),
          
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
      ungroup() %>% 
      select(-temp)
    
    return(sim_data)
    
  })
  
}


resultsfnc <- function(data, wt_denom_formula, grace_period_length = NULL, f_int = NULL){
  
  data_wts <- data %>% mutate(Wt = estimate_weights(
    data = data, 
    wt_denom_formula = wt_denom_formula, 
    grace_period_length = grace_period_length,
    f_int = f_int
  ))
  
  estimate <- data_wts %>% 
    group_by(interval, trt_A) %>% 
    summarise(cond_surv = 1-mean(Wt*dead)) %>% as.data.frame() %>% 
    group_by(trt_A) %>%
    mutate(CI = 1-cumprod(cond_surv)) 
  
  return(estimate)
}


estimate_weights <- function(data, wt_denom_formula, grace_period_length = NULL, f_int = NULL){
  
  data$prob_adherence <- 
    predict(glm(formula = wt_denom_formula, 
                data = data, 
                family = binomial()), 
            newdata = data, type="response")
  
  if (is.null(f_int)) {
    
    data_with_weights <- data %>%
      group_by(id) %>% 
      mutate(denom_wt = case_when(lag(num_intervals_nonadherent, 1) == grace_period_length ~ 1/prob_adherence,
                                  TRUE ~ 1),
             num_wt = case_when(lag(num_intervals_nonadherent, 1) == grace_period_length & adherent == 0 ~ 0,
                                TRUE ~ 1),
             Wt = num_wt/denom_wt) %>%
      mutate(Wt_cumprod = cumprod(Wt)) %>%
      ungroup()
    
  } else {
    
    f_int_cond_adh <- f_int / head(c(1, 1-cumsum(f_int)), -1)
    
    f_int_cond_nonadh <- 1 - f_int_cond_adh
    
    data_with_weights <- data %>%
      group_by(id) %>% 
      mutate(denom_wt = 1/prob_adherence,
             num_wt = case_when(adherent == 1 ~ f_int_cond_adh[num_intervals_nonadherent+1],
                                adherent == 0 ~ f_int_cond_adh[num_intervals_nonadherent+1]),
             Wt = num_wt/denom_wt) %>%
      mutate(Wt_cumprod = cumprod(Wt)) %>%
      ungroup()
    
  }
  
  #print(data_with_weights$Wt_cumprod %>% quantile(seq(0,1,0.1)))
  
  return(data_with_weights$Wt_cumprod)
}


bootstrapfnc <- function(data, wt_denom_formula, grace_period_length = NULL, f_int = NULL){
  
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
      f_int = f_int
    )
  )
}