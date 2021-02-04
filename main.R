source('./src/dependencies.R')
source('./src/functions.R')

#Set sim parameters
params <- list(
  sample_size = 1000,
  num_intervals = 120,
  p_Y = 0.01,
  p_AY = 0.005,
  p_BY = 0.005,
  p_nonadherence_A = 0.1,
  p_nonadherence_B = 0.2,
  p_nonadherence_baseline_1 = 0.05,
  p_nonadherence_baseline_2 = 0.025,
  p_nonadherence_baseline_3 = 0
)

num_cores <- 10

#Produce sim data

params$sample_size <- params$sample_size/num_cores

cl <- makeCluster(num_cores)
clusterEvalQ(cl, library(tidyverse))
clusterExport(cl, c('params'),
              envir=environment())

registerDoParallel(cl)

sim_data <- 
  foreach(i=1:num_cores, .combine=rbind)%dopar%{
    create_sim_data(params) %>% group_by(id) %>% filter(dead==0 | lag(dead==0)) %>% ungroup()
  } %>% mutate(id = case_when(interval==1 ~ 1, TRUE ~ 0)) %>% mutate(id = cumsum(id))

stopCluster(cl)

save(sim_data, file='./output/sim_data.Rdata')

#Estimate effect of treatment A versus treatment B 
#(all data is compatible with per-protocol regimes - no one is censored)
glm(dead ~ interval + trt_A, data=sim_data, family=binomial())$coef %>% exp

#Plot proportion of adherent individuals over time, stratified by treatment
ggplot(
  data = sim_data %>% group_by(interval, trt_A) %>% summarise(adherence=sum(adherent)/n()) %>% as.data.frame(),
  aes(y=adherence, x=interval, col=as.factor(trt_A))
  ) + geom_point()

#histogram of non-adherence over time, stratified by treatment
ggplot(
  data = 
    sim_data %>% 
      group_by(id) %>% 
      filter(adherent == 0 & (lag(adherent,1) == 1 | lag(adherent,2) == 1 | lag(adherent, 3) == 1)) %>% 
      ungroup() %>%
      select(c(interval, trt_A)) %>%
      as.data.frame(),
  aes(x=interval, fill=as.factor(trt_A))
) + geom_histogram(alpha=0.5)
