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
    
    return(sim_data)
    
  })
  
}
