requiredPackages = c(
  'tidyverse',
  'ggthemes',
  'data.table',
  'pbapply',
  'foreach',
  'doParallel',
  'ggpubr',
  'truncnorm',
  'ReIns',
  'speedglm'
)
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}
