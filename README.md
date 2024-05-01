# Grace periods in comparative effectiveness studies of sustained treatments

# Introduction
Here we provide the code to reproduce the analysis described in: Wanis KN, Sarvet AL, Wen L, Block JP, Rifas-Shiman SL, Robins JM, Young JG. Grace periods in comparative effectiveness studies of sustained treatments. Journal of the Royal Statistical Society Series A: Statistics in Society. https://doi.org/10.1093/jrsssa/qnae002.

# Organization
- `main_sim.R` — R file which creates a simulated dataset and then applies an inverse probability weighted estimator to compute risk under natural and stochastic grace period treatment strategies.
- `analysis_m_2.R, analysis_m_3.R, analysis_m_4.R` — R files which contain code to reproduce the analysis presented in the main text for grace periods of length m=2, m=3, and m=4.
- `src_sim` — Folder containing scripts including main functions for producing simulated data and for computing risk under grace period treatment strategies.
- `src`  — Folder containing scripts including main functions for computing risk under grace period treatment strategies.
- `produce_paper_tables.R` — R file which uses outputs from the files 'analysis_m_2.R', 'analysis_m_3.R', and 'analysis_m_4.R' to reproduce the tables shown in the paper.

# Correspondence
If you have any questions, comments, or discover an error, please contact Kerollos Wanis at knwanis@gmail.com.
