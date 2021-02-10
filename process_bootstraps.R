source('./src/dependencies.R')
source('./src/functions.R')

#differences

sapply(1:500, function(i){ bs_f_int[,i]$CI[length(bs_f_int[,i]$CI)] - 
    bs_f_int[,i]$CI[length(bs_f_int[,i]$CI)-1]}) %>% quantile(c(0.025,0.5,0.975))

sapply(1:500, function(i){ bs_f_int_A[,i]$CI[length(bs_f_int_A[,i]$CI)] - 
    bs_f_int_A[,i]$CI[length(bs_f_int_A[,i]$CI)-1]}) %>% quantile(c(0.025,0.5,0.975))

sapply(1:500, function(i){ bs_f_int_B[,i]$CI[length(bs_f_int_B[,i]$CI)] - 
    bs_f_int_B[,i]$CI[length(bs_f_int_B[,i]$CI)-1]}) %>% quantile(c(0.025,0.5,0.975))

sapply(1:500, function(i){ bs_f_int_50[,i]$CI[length(bs_f_int_50[,i]$CI)] - 
    bs_f_int_50[,i]$CI[length(bs_f_int_50[,i]$CI)-1]}) %>% quantile(c(0.025,0.5,0.975))

sapply(1:500, function(i){ bs_natural[,i]$CI[length(bs_natural[,i]$CI)] - 
    bs_natural[,i]$CI[length(bs_natural[,i]$CI)-1]}) %>% quantile(c(0.025,0.5,0.975))


#trt A means

sapply(1:500, function(i){ bs_f_int[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans() 

sapply(1:500, function(i){ bs_f_int_A[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans()

sapply(1:500, function(i){ bs_f_int_B[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans()

sapply(1:500, function(i){ bs_f_int_50[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans()

sapply(1:500, function(i){ bs_natural[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans()


#trt B means

sapply(1:500, function(i){ bs_f_int[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans() 

sapply(1:500, function(i){ bs_f_int_A[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans()

sapply(1:500, function(i){ bs_f_int_B[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans()

sapply(1:500, function(i){ bs_f_int_50[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans()

sapply(1:500, function(i){ bs_natural[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans()


#plots
#natural
p1 <- ggplot(data = data.frame(trt_A = c(rep(1,60), rep(0,60)), 
                         interval = rep(1:60, 2),
                         CI = c(sapply(1:500, function(i){ bs_natural[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans(),
                                sapply(1:500, function(i){ bs_natural[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans())), 
       aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_line(size=1.25) + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  xlab("") + ylab("") + 
  scale_linetype_manual(values = c("solid","dotted"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.key.size = unit(1.5,"line")) + guides(linetype = guide_legend(override.aes = list(size=0.75), nrow=1))

#trt A
p2 <- ggplot(data = data.frame(trt_A = c(rep(1,60), rep(0,60)), 
                               interval = rep(1:60, 2),
                               CI = c(sapply(1:500, function(i){ bs_f_int_A[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans(),
                                      sapply(1:500, function(i){ bs_f_int_A[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans())), 
             aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_line(size=1.25) + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  xlab("") + ylab("") + 
  scale_linetype_manual(values = c("solid","dotted"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(1.5,"line")) + guides(linetype = guide_legend(override.aes = list(size=0.75), nrow=1))

#trt B
p3 <- ggplot(data = data.frame(trt_A = c(rep(1,60), rep(0,60)), 
                               interval = rep(1:60, 2),
                               CI = c(sapply(1:500, function(i){ bs_f_int_B[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans(),
                                      sapply(1:500, function(i){ bs_f_int_B[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans())), 
             aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_line(size=1.25) + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  xlab("") + ylab("") + 
  scale_linetype_manual(values = c("solid","dotted"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(1.5,"line")) + guides(linetype = guide_legend(override.aes = list(size=0.75), nrow=1))

#90%
p4 <- ggplot(data = data.frame(trt_A = c(rep(1,60), rep(0,60)), 
                               interval = rep(1:60, 2),
                               CI = c(sapply(1:500, function(i){ bs_f_int[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans(),
                                      sapply(1:500, function(i){ bs_f_int[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans())), 
             aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_line(size=1.25) + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  xlab("") + ylab("") + 
  scale_linetype_manual(values = c("solid","dotted"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(1.5,"line")) + guides(linetype = guide_legend(override.aes = list(size=0.75), nrow=1))

#50%
p5 <- ggplot(data = data.frame(trt_A = c(rep(1,60), rep(0,60)), 
                               interval = rep(1:60, 2),
                               CI = c(sapply(1:500, function(i){ bs_f_int_50[,i]$CI[seq(2, 120, 2)]}) %>% rowMeans(),
                                      sapply(1:500, function(i){ bs_f_int_50[,i]$CI[seq(1, 119, 2)]}) %>% rowMeans())), 
             aes(x=interval, y=CI, linetype=as.factor(trt_A))) + geom_line(size=1.25) + 
  theme_tufte() + scale_x_continuous(breaks=seq(0, 120, 12)) +
  xlab("") + ylab("") + 
  scale_linetype_manual(values = c("solid","dotted"), 
                        labels = c("Medication B",
                                   "Medication A")) +
  ylim(c(0,0.4)) +
  theme(legend.title=element_blank(), text=element_text(size=18),
        axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        legend.position=c(.55,.25),
        legend.key.size = unit(1.5,"line")) + guides(linetype = guide_legend(override.aes = list(size=0.75), nrow=1))

mylegend<-g_legend(p1)


Fig1 <- ggarrange(p1,p2,p3,p4, common.legend = T, legend = 'bottom',
                  labels=c('A) Natural', 'B) f_int equal to trt A', 'C) f_int 90%'),
                  label.x = 0.2, label.y = 0.9,
                  font.label = list(size=11, face='italic'))
annotate_figure(Fig1,
                bottom = text_grob("Month", color = "black", size = 14),
                left = text_grob("Cumulative incidence of death", color = "black", rot = 90, size=14)
)

ggsave("1.png", dpi=600)
