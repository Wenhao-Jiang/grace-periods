source('./src/dependencies.R')
source('./src/functions.R')


######################################################
############## Generate simulated data  ##############
######################################################

#Set sim parameters
params <- list(
  sample_size = 100000,
  num_intervals = 60,
  p_Y = 0.02,
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


######################################################
############# Estimate treatment effect  #############
######################################################

#Estimate effect of treatment A versus treatment B using PH model
#(all data is compatible with per-protocol regimes - no one is censored)
glm(dead ~ interval + trt_A, data=sim_data, family=binomial())$coef %>% exp()

#Estimate effect of treatment A versus treatment B using IPW for the natural regime
estimate_natural <- 
  resultsfnc(data = sim_data %>% 
               group_by(id) %>% 
               mutate(lag_adherent_1 = lag(adherent,1,default=1), 
                      lag_adherent_2 = lag(adherent,2,default=1), 
                      lag_adherent_3 = lag(adherent,3,default=1)) %>%
               ungroup(), 
           wt_denom_formula = as.formula('adherent ~ interval*lag_adherent_1*lag_adherent_2*lag_adherent_3*trt_A'),
           grace_period_length = 3)

#Plot CI curves
ggplot(data = estimate_natural, 
       aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_line(size=1.25) + 
  xlab("Month")  + ylab("Cumulative incidence of death") + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  scale_linetype_manual(values = c("solid","dashed"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(2,"line")) + guides(linetype = guide_legend(override.aes = list(size=1.2)))

#bootstrap for the natural regime
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {library(tidyverse); library(data.table)})
clusterExport(cl, c('bootstrapfnc', 'resultsfnc', 'estimate_weights', 'sim_data'))

bs_natural <- 
  pbreplicate(500,
    bootstrapfnc(data = sim_data %>% 
                 group_by(id) %>% 
                 mutate(lag_adherent_1 = lag(adherent,1,default=1), 
                        lag_adherent_2 = lag(adherent,2,default=1), 
                        lag_adherent_3 = lag(adherent,3,default=1)) %>%
                 ungroup(), 
               wt_denom_formula = as.formula('adherent ~ interval*lag_adherent_1*lag_adherent_2*lag_adherent_3*trt_A'),
               grace_period_length = 3), cl=cl)
stopCluster(cl)
save(bs_natural, file='./output/bs_natural.Rdata')

sapply(1:500, function(i){ bs_natural[,i]$CI[length(bs_natural[,i]$CI)] - 
    bs_natural[,i]$CI[length(bs_natural[,i]$CI)-1]}) %>% quantile(c(0.025,0.975))


#Estimate effect of treatment A versus treatment B using IPW for a particular intervention distribution
estimate_f_int <- 
  resultsfnc(data = sim_data %>% 
               group_by(id) %>% 
               mutate(lag_adherent_1 = lag(adherent,1,default=1), 
                      lag_adherent_2 = lag(adherent,2,default=1), 
                      lag_adherent_3 = lag(adherent,3,default=1)) %>%
               ungroup(), 
             wt_denom_formula = as.formula('adherent ~ interval*lag_adherent_1*lag_adherent_2*lag_adherent_3*trt_A'),
             f_cond_int = c(0.8,0.8,0.8,1))

#Plot CI curves
ggplot(data = estimate_f_int, 
       aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_line(size=1.25) + 
  xlab("Month")  + ylab("Cumulative incidence of death") + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  scale_linetype_manual(values = c("solid","dashed"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(2,"line")) + guides(linetype = guide_legend(override.aes = list(size=1.2)))

#bootstrap for an f_int
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {library(tidyverse); library(data.table)})
clusterExport(cl, c('bootstrapfnc', 'resultsfnc', 'estimate_weights', 'sim_data'))

bs_f_int <- 
  pbreplicate(500,
              bootstrapfnc(data = sim_data %>% 
                             group_by(id) %>% 
                             mutate(lag_adherent_1 = lag(adherent,1,default=1), 
                                    lag_adherent_2 = lag(adherent,2,default=1), 
                                    lag_adherent_3 = lag(adherent,3,default=1)) %>%
                             ungroup(), 
                           wt_denom_formula = as.formula('adherent ~ interval*lag_adherent_1*lag_adherent_2*lag_adherent_3*trt_A'),
                           f_cond_int = c(0.9,0.9,0.9,1)), cl=cl)
stopCluster(cl)
save(bs_f_int, file='./output/bs_f_int.Rdata')

sapply(1:500, function(i){ bs_f_int[,i]$CI[length(bs_f_int[,i]$CI)] - 
    bs_f_int[,i]$CI[length(bs_f_int[,i]$CI)-1]}) %>% quantile(c(0.025,0.975))
