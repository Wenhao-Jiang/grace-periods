# Reproduces all tables in the paper

source('./src/dependencies.R')
source('./src/functions.R')

load('cohort.RData')

pacman::p_load(tableone,tables,pbapply,xtable)

#################################
#################################
# Supplementary table summarizing baseline characteristics
#################################
#################################

myVars <- c("female",
            "bariatric_pre",
            "Chronic_Kidney_Disease_pre",
            "Liver_Cirrhosis_pre",
            "T2D_pre",
            "T1D_pre",
            "age")
binaryVars <- c("female",
                "bariatric_pre",
                "Chronic_Kidney_Disease_pre",
                "Liver_Cirrhosis_pre",
                "T2D_pre",
                "T1D_pre")
continuous <- c("age")
tab0 <- CreateTableOne(vars = myVars, data = cohort %>% 
                         mutate(female = as.factor(female),
                                bariatric_pre = as.factor(bariatric_pre),
                                Chronic_Kidney_Disease_pre = as.factor(Chronic_Kidney_Disease_pre),
                                Liver_Cirrhosis_pre = as.factor(Liver_Cirrhosis_pre),
                                T2D_pre = as.factor(T2D_pre),
                                T1D_pre = as.factor(T1D_pre)) %>%
                         filter(month==1) %>% {mutate(., initiated_med=ifelse(.$initiated_med=='ACE',
                                                                              paste0("\\specialcell{Initiated ACEI \\\\ N=",table(.$initiated_med)[[1]],"}"),
                                                                              paste0("\\specialcell{Initiated thiazide \\\\ N=",table(.$initiated_med)[[2]],"}")))},
                       strata = "initiated_med")
tab0.p <- print(tab0, nonnormal=continuous, quote=FALSE, test=FALSE, noSpaces=TRUE, printToggle = FALSE, contDigits=1, catDigits=1)

tab1.format <- t(t(as.matrix(tab0.p)) %>% as.data.frame %>% rename(
  "Female gender"="female = 1 (%)",
  "Bariatric diagnosis"="bariatric_pre = 1 (%)",
  "Chronic kidney disease"="Chronic_Kidney_Disease_pre = 1 (%)",
  "Liver Cirrhosis"="Liver_Cirrhosis_pre = 1 (%)",
  "Type 2 diabetes"="T2D_pre = 1 (%)",
  "Type 1 diabetes"="T1D_pre = 1 (%)",
  "Age (years)"="age (median [IQR])"
))

named1 <- rownames(tab1.format)
tags1 <- grepl("^ ", rownames(tab1.format))
rownames(tab1.format) <- c(ifelse(tags1==FALSE, named1, paste("\\hskip .5cm", named1, sep=' ')))
print(xtable(tab1.format, align=c("lrr")),
      type="latex", sanitize.text.function = function(x){x}, file="table1.tex",
      floating=FALSE, tabular.environment="longtable")

#################################
#################################
#Table 1
#################################
#################################

load("./save_results_ACE_stochastic_0.95_constant_m_3.Rdata")
load("./save_results_thiazide_stochastic_0.95_constant_m_3.Rdata")
load("./save_results_ACE_natural_m_3.Rdata")
load("./save_results_thiazide_natural_m_3.Rdata")
load("./save_results_ACE_stochastic_0.9_constant_m_3.Rdata")
load("./save_results_thiazide_stochastic_0.9_constant_m_3.Rdata")
load("./save_results_ACE_stochastic_0.75_constant_m_3.Rdata")
load("./save_results_thiazide_stochastic_0.75_constant_m_3.Rdata")
load("./save_results_ACE_stochastic_0.5_constant_m_3.Rdata")
load("./save_results_thiazide_stochastic_0.5_constant_m_3.Rdata")

estimates_m_3 <- data.frame(
  protocol = c("ACE stochastic 0.95 constant", 
               "thiazide stochastic 0.95 constant",
               "ACE natural",
               "thiazide natural", 
               "ACE stochastic 0.9 constant", 
               "thiazide stochastic 0.9 constant",
               "ACE stochastic 0.75 constant", 
               "thiazide stochastic 0.75 constant", 
               "ACE stochastic 0.5 constant", 
               "thiazide stochastic 0.5 constant"),
  estimate = c(
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.95_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.95_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_natural_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_natural_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.9_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.9_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.75_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.75_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.5_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.5_constant_m_3, 0.95)$outcome
  ),
  se = c(
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.95_constant_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.95_constant_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_natural_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_natural_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.9_constant_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.9_constant_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.75_constant_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.75_constant_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.5_constant_m_3, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.5_constant_m_3, 0.95)$se
  )
) %>% mutate(lower = estimate-1.96*se,
             upper = estimate+1.96*se) %>%
  mutate(non_null = case_when(sign(lower) == sign(upper) ~ TRUE,
                              TRUE ~ FALSE))

estimates_m_3 %>% arrange(estimate) %>% mutate_if(is.numeric, round, 3)

#################################
#################################
#Table S3
#################################
#################################

load("./save_results_ACE_stochastic_0.95_constant_m_2.Rdata")
load("./save_results_thiazide_stochastic_0.95_constant_m_2.Rdata")
load("./save_results_ACE_natural_m_2.Rdata")
load("./save_results_thiazide_natural_m_2.Rdata")
load("./save_results_ACE_stochastic_0.9_constant_m_2.Rdata")
load("./save_results_thiazide_stochastic_0.9_constant_m_2.Rdata")
load("./save_results_ACE_stochastic_0.75_constant_m_2.Rdata")
load("./save_results_thiazide_stochastic_0.75_constant_m_2.Rdata")
load("./save_results_ACE_stochastic_0.5_constant_m_2.Rdata")
load("./save_results_thiazide_stochastic_0.5_constant_m_2.Rdata")

estimates_m_2 <- data.frame(
  protocol = c("ACE stochastic 0.95 constant", 
               "thiazide stochastic 0.95 constant",
               "ACE natural",
               "thiazide natural", 
               "ACE stochastic 0.9 constant", 
               "thiazide stochastic 0.9 constant",
               "ACE stochastic 0.75 constant", 
               "thiazide stochastic 0.75 constant", 
               "ACE stochastic 0.5 constant", 
               "thiazide stochastic 0.5 constant"),
  estimate = c(
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.95_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.95_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_natural_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_natural_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.9_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.9_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.75_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.75_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.5_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.5_constant_m_2, 0.95)$outcome
  ),
  se = c(
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.95_constant_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.95_constant_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_natural_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_natural_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.9_constant_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.9_constant_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.75_constant_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.75_constant_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.5_constant_m_2, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.5_constant_m_2, 0.95)$se
  )
) %>% mutate(lower = estimate-1.96*se,
             upper = estimate+1.96*se) %>%
  mutate(non_null = case_when(sign(lower) == sign(upper) ~ TRUE,
                              TRUE ~ FALSE))

estimates_m_2 %>% arrange(estimate) %>% mutate_if(is.numeric, round, 3)

#################################
#################################
#Table S5
#################################
#################################

load("./save_results_ACE_stochastic_0.95_constant_m_4.Rdata")
load("./save_results_thiazide_stochastic_0.95_constant_m_4.Rdata")
load("./save_results_ACE_natural_m_4.Rdata")
load("./save_results_thiazide_natural_m_4.Rdata")
load("./save_results_ACE_stochastic_0.9_constant_m_4.Rdata")
load("./save_results_thiazide_stochastic_0.9_constant_m_4.Rdata")
load("./save_results_ACE_stochastic_0.75_constant_m_4.Rdata")
load("./save_results_thiazide_stochastic_0.75_constant_m_4.Rdata")
load("./save_results_ACE_stochastic_0.5_constant_m_4.Rdata")
load("./save_results_thiazide_stochastic_0.5_constant_m_4.Rdata")

estimates_m_4 <- data.frame(
  protocol = c("ACE stochastic 0.95 constant", 
               "thiazide stochastic 0.95 constant",
               "ACE natural",
               "thiazide natural", 
               "ACE stochastic 0.9 constant", 
               "thiazide stochastic 0.9 constant",
               "ACE stochastic 0.75 constant", 
               "thiazide stochastic 0.75 constant", 
               "ACE stochastic 0.5 constant", 
               "thiazide stochastic 0.5 constant"),
  estimate = c(
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.95_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.95_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_natural_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_natural_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.9_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.9_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.75_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.75_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.5_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.5_constant_m_4, 0.95)$outcome
  ),
  se = c(
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.95_constant_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.95_constant_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_natural_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_natural_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.9_constant_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.9_constant_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.75_constant_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.75_constant_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_ACE_stochastic_0.5_constant_m_4, 0.95)$se,
    trim_weights_compute_risk(0.99, save_results_thiazide_stochastic_0.5_constant_m_4, 0.95)$se
  )
) %>% mutate(lower = estimate-1.96*se,
             upper = estimate+1.96*se) %>%
  mutate(non_null = case_when(sign(lower) == sign(upper) ~ TRUE,
                              TRUE ~ FALSE))

estimates_m_4 %>% arrange(estimate) %>% mutate_if(is.numeric, round, 3)

#################################
#################################
#Table 2
#################################
#################################

estimates_diff_m_3 <- data.frame(
  protocol = c("stochastic 0.95 constant",  
               "natural",
               "stochastic 0.9 constant", 
               "stochastic 0.75 constant", 
               "stochastic 0.5 constant"),
  estimate = c(
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.95_constant_m_3, save_results_thiazide_stochastic_0.95_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_natural_m_3, save_results_thiazide_natural_m_3, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.9_constant_m_3, save_results_thiazide_stochastic_0.9_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.75_constant_m_3, save_results_thiazide_stochastic_0.75_constant_m_3, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.5_constant_m_3, save_results_thiazide_stochastic_0.5_constant_m_3, 0.95)$outcome
  ),
  se = c(
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.95_constant_m_3, save_results_thiazide_stochastic_0.95_constant_m_3, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_natural_m_3, save_results_thiazide_natural_m_3, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.9_constant_m_3, save_results_thiazide_stochastic_0.9_constant_m_3, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.75_constant_m_3, save_results_thiazide_stochastic_0.75_constant_m_3, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.5_constant_m_3, save_results_thiazide_stochastic_0.5_constant_m_3, 0.95)$se
  )) %>% mutate(lower = estimate-1.96*se,
                upper = estimate+1.96*se) %>%
  mutate(non_null = case_when(sign(lower) == sign(upper) ~ TRUE,
                              TRUE ~ FALSE))

estimates_diff_m_3 %>% arrange(estimate) %>% mutate_if(is.numeric, round, 3)

#################################
#################################
#Table S4
#################################
#################################

estimates_diff_m_2 <- data.frame(
  protocol = c("stochastic 0.95 constant",  
               "natural",
               "stochastic 0.9 constant", 
               "stochastic 0.75 constant", 
               "stochastic 0.5 constant"),
  estimate = c(
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.95_constant_m_2, save_results_thiazide_stochastic_0.95_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_natural_m_2, save_results_thiazide_natural_m_2, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.9_constant_m_2, save_results_thiazide_stochastic_0.9_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.75_constant_m_2, save_results_thiazide_stochastic_0.75_constant_m_2, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.5_constant_m_2, save_results_thiazide_stochastic_0.5_constant_m_2, 0.95)$outcome
  ),
  se = c(
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.95_constant_m_2, save_results_thiazide_stochastic_0.95_constant_m_2, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_natural_m_2, save_results_thiazide_natural_m_2, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.9_constant_m_2, save_results_thiazide_stochastic_0.9_constant_m_2, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.75_constant_m_2, save_results_thiazide_stochastic_0.75_constant_m_2, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.5_constant_m_2, save_results_thiazide_stochastic_0.5_constant_m_2, 0.95)$se
  )) %>% mutate(lower = estimate-1.96*se,
                upper = estimate+1.96*se) %>%
  mutate(non_null = case_when(sign(lower) == sign(upper) ~ TRUE,
                              TRUE ~ FALSE))

estimates_diff_m_2 %>% arrange(estimate) %>% mutate_if(is.numeric, round, 3)

#################################
#################################
#Table S6
#################################
#################################

estimates_diff_m_4 <- data.frame(
  protocol = c("stochastic 0.95 constant",  
               "natural",
               "stochastic 0.9 constant", 
               "stochastic 0.75 constant", 
               "stochastic 0.5 constant"),
  estimate = c(
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.95_constant_m_4, save_results_thiazide_stochastic_0.95_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_natural_m_4, save_results_thiazide_natural_m_4, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.9_constant_m_4, save_results_thiazide_stochastic_0.9_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.75_constant_m_4, save_results_thiazide_stochastic_0.75_constant_m_4, 0.95)$outcome,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.5_constant_m_4, save_results_thiazide_stochastic_0.5_constant_m_4, 0.95)$outcome
  ),
  se = c(
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.95_constant_m_4, save_results_thiazide_stochastic_0.95_constant_m_4, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_natural_m_4, save_results_thiazide_natural_m_4, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.9_constant_m_4, save_results_thiazide_stochastic_0.9_constant_m_4, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.75_constant_m_4, save_results_thiazide_stochastic_0.75_constant_m_4, 0.95)$se,
    trim_weights_compute_risk_difference(0.99, save_results_ACE_stochastic_0.5_constant_m_4, save_results_thiazide_stochastic_0.5_constant_m_4, 0.95)$se
  )) %>% mutate(lower = estimate-1.96*se,
                upper = estimate+1.96*se) %>%
  mutate(non_null = case_when(sign(lower) == sign(upper) ~ TRUE,
                              TRUE ~ FALSE))

estimates_diff_m_4 %>% arrange(estimate) %>% mutate_if(is.numeric, round, 3)
