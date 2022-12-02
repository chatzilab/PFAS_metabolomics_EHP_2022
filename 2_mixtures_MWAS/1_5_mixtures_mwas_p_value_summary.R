# P-value Summary
library("jag2")
library("tidyverse")
source(here::here("0_project_setup", "!directories.r"))
# Read in data and calculate summary stats ----------------------------------
mwas_results <- read_csv(
  file = fs::path(dir_results_mixtures, 
                  "all_pfas_mixtures_results_hyper_g_v3.csv")) %>% 
  as_tibble()

# Calculate summaries
(mwas_summary <- mwas_results %>% 
    mutate(associated_bci = if_else(lcl_beta > 0 | ucl_beta < 0,
                                    "Associated", 
                                    "Not Associated")) %>% 
    group_by(cohort, mixture, exposure) %>% 
    summarise(n_features = length(feature), 
              percent_significant_p05 = jag2::npct(significance, 
                                                   "p < 0.05", 
                                                   n.digits = 2), 
              percent_significant_q05 = jag2::npct(significancefdr, 
                                                   "q < 0.05", 
                                                   n.digits = 2), 
              pct_associated = jag2::npct(associated_bci, 
                                          "Associated", 
                                          n.digits = 2)))

# pivot summary wider
mwas_summary_w <- pivot_wider(mwas_summary, 
                              id_cols = c(exposure, mixture, n_features),
                              names_from = cohort, 
                              values_from = c(percent_significant_p05, 
                                              percent_significant_q05, 
                                              pct_associated)) %>% 
  select(exposure, mixture, contains("SOL"), contains("CHS"), contains("Pooled"))

# Save Summary Data to csv -----------------------------------
write_excel_csv(mwas_summary_w, 
                file = fs::path(
                  dir_reports, 
                  "Summary of Mixtures MWAS Results_hyper_g_v3.csv"))
