source('./src_sim/dependencies.R')
source('./src_sim/functions.R')


######################################################
############## Generate simulated data  ##############
######################################################

#Set sim parameters
params <- list(
  sample_size = 50000,
  num_intervals = 24,
  p_Y = 0.03,
  p_AY = 0.005,
  p_BY = 0.005,
  p_nonadherence_A = 0.25,
  p_nonadherence_B = 0.50,
  p_nonadherence_baseline_1 = 0.20,
  p_nonadherence_baseline_2 = 0.20,
  p_nonadherence_baseline_3 = 0,
  p_nonadherence_baseline_4 = 0,
  p_nonadherence_baseline_5 = 0,
  p_nonadherence_baseline_6 = 0
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

save(sim_data, file='./output_sim/sim_data.Rdata')


######################################################
############## View adherence  #######################
######################################################

#Plot proportion of adherent individuals over time, stratified by treatment
ggplot(
  data = sim_data %>% group_by(interval, trt_A) %>% summarise(adherence=sum(adherent)/n()) %>% as.data.frame(),
  aes(y=adherence, x=interval, col=as.factor(trt_A))
) + geom_point()


######################################################
############# Estimate treatment effect  #############
######################################################

weight_formula = as.formula('adherent ~ I(num_intervals_without_trt>0)*trt_A + trt_A*I(interval==1) + I(num_intervals_without_trt>=3)')

#Estimate effect of treatment A versus treatment B using PH model
glm(dead ~ interval + trt_A, data=sim_data, family=binomial())$coef %>% exp()

#Estimate effect of treatment A versus treatment B using IPW for the natural regime
estimate_natural <- 
  resultsfnc(data = sim_data, 
             wt_denom_formula = weight_formula,
             grace_period_length = 3)

#Plot CI curves
g1 <- ggplot(data = estimate_natural, 
             aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_smooth(size=1.25, col='black', se=F) + 
  xlab("Month")  + ylab("Cumulative incidence of death") + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  scale_linetype_manual(values = c("solid","dotted"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(2,"line")) + guides(linetype = guide_legend(override.aes = list(size=1.2)))


#Estimate effect of treatment A versus treatment B using IPW for a particular intervention distribution
estimate_f_int <- 
  resultsfnc(data = sim_data, 
             wt_denom_formula = weight_formula,
             f_cond_int = c(1,1,1)) # static always treat
             #f_cond_int = c(0.9179567, 0.918423, 0.9241418)) # truncated exponential on (0,1) with rate paramter = 10
             #f_cond_int = c(0.5552792, 0.5897977, 0.6791787)) # truncated exponential on (0,1) with rate paramter = 3
             #f_cond_int = c(0.4063653, 0.6845376, 0.8639574)) # truncated normal on (0,1) with mean=0.25, sd=0.25
             #f_cond_int = c(0.5817036, 0.8006517, 0.917605)) # truncated normal on (0,1) with mean=0.1, sd=0.25
             #f_cond_int = c(0.75, 0.3333333, 0.5)) 
             #f_cond_int = c(0.9, 0.3333333, 0.5)) 

#Plot CI curves
g2 <- ggplot(data = estimate_f_int, 
             aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_smooth(size=1.25, col='black', se=F) + 
  xlab("Month")  + ylab("Cumulative incidence of death") + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  scale_linetype_manual(values = c("solid","dotted"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(2,"line")) + guides(linetype = guide_legend(override.aes = list(size=1.2)))


Fig1 <- ggarrange(g1 + rremove("ylab") + rremove("xlab"), g2 + rremove("ylab") + rremove("xlab"),
                  common.legend = T, legend = 'bottom',
                  labels=c('A) Natural grace period strategy', 'B) Stochastic grace period strategy'),
                  label.x = -0.1, label.y = 0.9,
                  font.label = list(size=11, face='italic'))
annotate_figure(Fig1,
                bottom = text_grob("Month", color = "black", size = 14),
                left = text_grob("Cumulative incidence of death", color = "black", rot = 90, size=14)
)

ggsave("simple_example.png", dpi=600)

#bootstrap for the natural regime
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {library(tidyverse); library(data.table); library(speedglm)})
clusterExport(cl, c('bootstrapfnc', 'resultsfnc', 'estimate_weights', 'sim_data', 'weight_formula'))

bs_natural <- 
  pbreplicate(500,
              bootstrapfnc(data = sim_data, 
                           wt_denom_formula = weight_formula,
                           grace_period_length = 3), cl=cl)
stopCluster(cl)
save(bs_natural, file='./output_sim/bs_natural.Rdata')

#compute bootstrap 95% confidence interval
sapply(1:500, function(i){ bs_natural[,i]$CI[length(bs_natural[,i]$CI)] - 
    bs_natural[,i]$CI[length(bs_natural[,i]$CI)-1]}) %>% quantile(c(0.025,0.975))

#bootstrap for an f_int
cl <- makeCluster(num_cores)
clusterEvalQ(cl, {library(tidyverse); library(data.table); library(speedglm)})
clusterExport(cl, c('bootstrapfnc', 'resultsfnc', 'estimate_weights', 'sim_data', 'weight_formula'))

bs_f_int <- 
  pbreplicate(500,
              bootstrapfnc(data = sim_data, 
                           wt_denom_formula = weight_formula,
                           f_cond_int = c(0.5,0.75,0.9,0.9,0.9,0.9)), cl=cl)
stopCluster(cl)
save(bs_f_int, file='./output_sim/bs_f_int.Rdata')

#compute bootstrap 95% confidence interval
sapply(1:500, function(i){ bs_f_int[,i]$CI[length(bs_f_int[,i]$CI)] - 
    bs_f_int[,i]$CI[length(bs_f_int[,i]$CI)-1]}) %>% quantile(c(0.025,0.975))
