requiredPackages = c(
  'tidyverse',
  'ggthemes',
  'data.table',
  'pbapply',
  'foreach',
  'doParallel'
)
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}