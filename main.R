source('./src/dependencies.R')
source('./src/functions.R')


######################################################
############## Generate simulated data  ##############
######################################################

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
clusterExport(cl, c('params'))

registerDoParallel(cl)

sim_data <- 
  foreach(i=1:num_cores, .combine=rbind)%dopar%{
    create_sim_data(params)
  } %>% mutate(id = case_when(interval==1 ~ 1, TRUE ~ 0)) %>% mutate(id = cumsum(id))

stopCluster(cl)

save(sim_data, file='./output/sim_data.Rdata')


######################################################
############## View/describe adherence  ##############
######################################################

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


######################################################
############# Estimate treatment effect  #############
######################################################

#Estimate effect of treatment A versus treatment B using PH model
#(all data is compatible with per-protocol regimes - no one is censored)
glm(dead ~ interval + trt_A, data=sim_data, family=binomial())$coef %>% exp()

#Estimate effect of treatment A versus treatment B using IPW
estimate <- 
  resultsfnc(data = sim_data %>% 
             group_by(id) %>% 
             mutate(lag_adherent_1 = lag(adherent,1), lag_adherent_2 = lag(adherent,2), lag_adherent_3 = lag(adherent,3)) %>%
             replace_na(list(lag_adherent_1 = 1, lag_adherent_2 = 1, lag_adherent_3 = 1)) %>%
             ungroup(), 
           wt_denom_formula = as.formula('adherent ~ lag_adherent_1*lag_adherent_2*lag_adherent_3*trt_A'),
           grace_period_length = 3)

#Plot CI curves
ggplot(data = estimate, 
       aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_line(size=1.25) + 
  xlab("Month")  + ylab("Cumulative incidence of death") + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  scale_linetype_manual(values = c("solid","dashed"), 
                        labels = c("Treatment A",
                                   "Treatment B")) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(2,"line")) + guides(linetype = guide_legend(override.aes = list(size=1.2)))

#bootstrap
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {library(tidyverse); library(data.table)})
clusterExport(cl, c('bootstrapfnc', 'resultsfnc', 'estimate_weights', 'sim_data'))

bs <- 
  pbreplicate(500,
    bootstrapfnc(data = sim_data %>% 
                 group_by(id) %>% 
                 mutate(lag_adherent_1 = lag(adherent,1), lag_adherent_2 = lag(adherent,2), lag_adherent_3 = lag(adherent,3)) %>%
                 replace_na(list(lag_adherent_1 = 1, lag_adherent_2 = 1, lag_adherent_3 = 1)) %>%
                 ungroup(), 
               wt_denom_formula = as.formula('adherent ~ lag_adherent_1*lag_adherent_2*lag_adherent_3*trt_A'),
               grace_period_length = 3), cl=cl)
stopCluster(cl)

sapply(1:500, function(i){ bs[,i]$CI[length(bs[,i]$CI)] - bs[,i]$CI[length(bs[,i]$CI)-1]}) %>% quantile(c(0.025,0.975))
