library(tidyverse)
library(janitor)

# Source setup scripts
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!set_exposure_outcome_vars.R"))

#Key for superpathways
superpathwaykey <- readxl::read_xlsx(
  fs::path(dir_data, 
           "Supporting files",  
           "superpathway_key_sept_21.xlsx")) %>% 
  rename(path = pathway)

# Get list of all results folders ------------------------
dir_results_exposures <- fs::path(dir_results_mum_mixtures, 
                                  "Mixture effect hyper_g") 

cohort_mixture = c("solar/sol_all_pfas_p05", 
                   "solar/sol_pfsas_p05",
                   "solar/sol_pfcas_p05",
                   "chs/chs_all_pfas_p05",
                   "chs/chs_pfsas_p05",
                   "chs/chs_pfcas_p05", 
                   "pooled/pooled_all_pfas_p05")

cohort_mixture2 = c("sol_all_pfas", "sol_pfsas", "sol_pfcas",
                    "chs_all_pfas", "chs_pfsas", "chs_pfcas", 
                    "pooled_all_pfas")


dir_results_exposures_chrt_mode <- map(dir_results_exposures,
                                       ~fs::path(.x, 
                                                 cohort_mixture)) %>% 
  unlist() 

# 0) Load Mummichog RDS files --------------------------------------------------
pthwy_pvals <- map(dir_results_exposures_chrt_mode, 
                   ~read_csv(fs::path(.x,
                                      "mummichog_integ_pathway_enrichment.csv"))) %>% 
  modify(~clean_names(.x))

pthwy_ecs <-  map(dir_results_exposures_chrt_mode, 
                  ~read_csv(fs::path(.x,
                                     "mummichog_pathway_enrichment.csv"))) %>% 
  modify(~clean_names(.x))



names(pthwy_pvals) <- cohort_mixture2
names(pthwy_ecs) <- cohort_mixture2

# Combine dataframes
mum_res_lst <- map2(pthwy_pvals, pthwy_ecs, ~tidylog::full_join(.x, .y))


# 1) Combine cohorts ---------------------
mum_pw <- mum_res_lst %>% 
  bind_rows(., .id = "cohort_mixture") %>% 
  janitor::clean_names() %>% 
  rename(path = x1) %>%
  select(cohort_mixture, everything())

# Get columns for PFAS, cohort,
mum_pw1 <- mum_pw %>% 
  mutate( 
    cohort = str_split_fixed(cohort_mixture, "_", 2)[,1], 
    mixture = str_split_fixed(cohort_mixture, "_", 2)[,2],
    enrichment = hits_sig/hits_total, 
    neg_logp = -log10(combined_pvals),
    path_2  = str_replace(path, "metabolism", "met.") %>% 
      str_replace("Metabolism", "met.") %>% 
      str_replace(., " pathway", "")) %>% 
  select(cohort, mixture, everything(), -pathway_number, 
         -c(total_size:sig_hits), -fet, -ease, -gamma)

# Pivot wider on cohort
mum_pw_w1 <- mum_pw1 %>% 
  group_by(mixture) %>% 
  nest() %>% 
  mutate(data_w = map(data, 
                      ~pivot_wider(., 
                                   id_cols = c(path, 
                                               path_2,
                                               pathway_total,
                                               hits_total), 
                                   names_from = cohort, 
                                   values_from = c(mummichog_pvals:combined_pvals, 
                                                   hits_sig:neg_logp)))) %>% 
  select(-data) %>% 
  unnest(data_w) %>% 
  ungroup()

# 2) Perform meta analysis of p values ----------------------------------
wgt_sol = sqrt(312)
wgt_chs = sqrt(137)

mum_pw_w1 <- mum_pw_w1 %>% 
  mutate(pval_sol = if_else(is.na(combined_pvals_sol), .99, 
                            combined_pvals_sol),
         pval_chs = if_else(is.na(combined_pvals_chs), .99, 
                            combined_pvals_chs), 
         hits_sig_chs = replace_na(hits_sig_chs, 0), 
         hits_sig_sol = replace_na(hits_sig_sol, 0), 
         hits_sig_pooled = replace_na(hits_sig_pooled, 0)) %>%
  rowwise() %>% 
  mutate(pval_meta = metap::sumz(p = c_across(pval_sol:pval_chs), 
                                 weights = c(wgt_sol, wgt_chs))$p[[1]] %>% 
           if_else(is.na(combined_pvals_sol) | is.na(combined_pvals_chs), 
                   NA_real_, .),
         
         enrichment_meta = 
           ((enrichment_sol*wgt_sol)+(enrichment_chs*wgt_chs))/(wgt_sol+wgt_chs), 
         neg_logp_meta = -log10(pval_meta), 
         q_meta = p.adjust(pval_meta), 
         sig_meta = if_else(pval_meta < 0.05, "Sig.", "Not. Sig"), 
         sig_fdr = if_else(q_meta < 0.2, "Sig. FDR < 0.2", "Not. Sig"), 
         hits_sig_meta = (hits_sig_sol + hits_sig_chs)/2
  ) %>% 
  ungroup()

# select only pathways which were reported in both cohorts
mum_pw_w_reduced <- mum_pw_w1 %>% 
  tidylog::filter(pathway_total > 2, 
                  hits_sig_sol > 3 |
                    hits_sig_chs > 3) %>%
  mutate(sig = case_when(combined_pvals_sol < 0.05 & 
                           combined_pvals_chs < 0.05 ~ "Sig. Both Cohorts", 
                         combined_pvals_sol < 0.05 ~ "Sig. SOLAR Only", 
                         combined_pvals_chs < 0.05 ~ "Sig. CHS Only", 
                         TRUE ~ "Not Significant"), 
         sig_overall_p = sig_meta) 


# Combine with superpathway metadata
mum_pw_final <- mum_pw_w_reduced %>% 
  tidylog::left_join(superpathwaykey)


# Clean Environment
# rm(mum_pw, mum_pw1, mum_pw_w1, wgt_chs, wgt_sol)

# Save Data 
write_rds(mum_pw_final,
          fs::path(dir_results_mum_mixtures,
                   "Mixture effect hyper_g",
                   "SOL CHS PFAS Mummichog wide sig PW_v3.RDS"))

# Save csv 
write_csv(mum_pw_final %>% filter(mixture == "all_pfas"),
          fs::path(dir_results_mum_mixtures,
                   "Mixture effect hyper_g",
                   "SOL CHS PFAS Mummichog wide sig PW_v3.csv"))




# 3) Create long dataframe ----------------------
sol_only <- mum_pw_final %>% 
  select(everything(),
         -contains("chs"), 
         -contains("meta"),
         -contains("pooled"),
         sig_meta) %>% 
  rename_all(~str_remove(., "_sol")) %>%
  rename(sig_cohort = sig) %>%
  mutate(cohort = "solar")

chs_only <- mum_pw_final %>% 
  select(everything(), 
         -contains("sol"),
         -contains("meta"),
         -contains("pooled"),
         sig_meta) %>% 
  rename_all(~str_remove(., "_chs")) %>% 
  rename(sig_cohort = sig) %>%
  mutate(cohort = "chs")


meta_only <- mum_pw_final %>% 
  select(-contains("sol"), 
         -contains("chs"), 
         -contains("pooled"),
         -q_meta) %>% 
  rename(meta_p_sig = sig_meta) %>%
  rename_all(~str_remove(., "_meta")) %>% 
  rename(sig_meta = meta_p_sig, 
         sig_cohort = sig) %>%
  mutate(cohort = "meta")


pooled_only <- mum_pw_final %>% 
  select(everything(), 
         -contains("sol"),
         -contains("chs"),
         -contains("meta"),
         sig_meta) %>% 
  rename_all(~str_remove(., "_pooled")) %>% 
  rename(sig_cohort = sig) %>%
  mutate(cohort = "pooled")

# Bind cohorts
mum_pw_final_long <- bind_rows(sol_only, chs_only, meta_only, pooled_only) %>% 
  select(cohort, everything())

table(mum_pw_final_long$sig_cohort, 
      mum_pw_final_long$cohort,
      mum_pw_final_long$mixture)



# Save Data 
write_rds(mum_pw_final_long,
          fs::path(dir_results_mum_mixtures,
                   "Mixture effect hyper_g",
                   "SOL CHS PFAS Mummichog long sig PW_v3.RDS"))

