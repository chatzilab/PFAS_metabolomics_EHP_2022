# Get cpds from mum

library(tidyverse)
library(ggplot2)
library(cowplot)
library(here)
library(ggrepel)
library(fs)
library(janitor)
ggplot2::theme_set(cowplot::theme_cowplot())

# Source setup scripts
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!set_exposure_outcome_vars.R"))
# source(here::here("2_3_Combine_mum_pw_between_chrts_mixtures.R"))

#1) Get Empirical compound to Pathway Key from metaboanalyst -----------------------
mum_pws <- read_rds(fs::path(dir_results_mum_mixtures,
                             "Mixture effect hyper_g",
                             "SOL CHS PFAS Mummichog long sig PW_v3.RDS")) %>% 
  filter(!is.na(ec_hits))

#2) Read in mz/rt key ------------------------------------------
# (Note: doesn't matter which cohort you use, this file is the same)
mzrtkey  <- read_csv(fs::path(dir_results_mum_mixtures,
                              "Mixture effect hyper_g",
                              "solar",
                              "sol_all_pfas_p05",
                              "mummichog_matched_compound_all.csv")) %>% 
  clean_names() %>%
  mutate(feature = str_c(query_mass, 
                         retention_time, 
                         sep = "_"))


# 3) Create pw to ec dataset -------------------------------------
# Split ecs based on ";", turn into dataframe 
pw_ec_df <- str_split(mum_pws$ec_hits, ";") %>% 
  enframe() %>% 
  rename(empirical_compound = value)

# Bind pw_ec_df with full data
pw_ec_key <- bind_cols(mum_pws, pw_ec_df)

# Unnest data to get a unique row for each empirical compound
pw_ec_key2 <- pw_ec_key %>% 
  unnest(empirical_compound) %>% 
  select(cohort, path, path_2, empirical_compound)

# Check for distinct values
ecd_pw_key_final <- pw_ec_key2 %>% 
  tidylog::distinct(path, path_2, empirical_compound) %>% 
  tidylog::group_by(empirical_compound) %>% 
  tidylog::summarise(path = str_c(path, collapse = "; "), 
                     path_2 = str_c(path_2, collapse = "; ")) %>% 
  ungroup()

# 4) Combine mzrt key and ecd_pw_key -----------------------------------------
ecd_pw_key <- tidylog::left_join(ecd_pw_key_final, 
                                 mzrtkey, 
                                 by = c("empirical_compound")) %>% 
  tidylog::filter(!is.na(feature)) %>% 
  rename(name = feature) %>% 
  tidylog::distinct()

# 4) Merge with Compound Names ----------------------------------------------------------
# Read in hand curated list of molecule names (this was from previous version 
# of data analysis, but the cpd/name key is correct). This data was created in "other_scripts/Convert Mummichog cpd ids to cpd names.R"
cpd_name_key <- read_rds(
  fs::path(dir_data_local,
           "Supporting Files", 
           "kegg_cpd_name_key.rds")) %>% 
  select(matched_compound, 
         met_name) %>% 
  distinct()

# left join with key from current analysis
ecd_pw_key_2 <- ecd_pw_key %>% 
  tidylog::left_join(cpd_name_key) %>% 
  select(path, path_2, met_name, everything()) %>% 
  rename(pathway = path) 


# Create final empirical compound pathway key
ecd_pw_key_final <- ecd_pw_key_2 %>% 
  tidylog::group_by(pathway, 
                    path_2,
                    empirical_compound,
                    matched_form,
                    query_mass,
                    retention_time,
                    mass_diff,
                    name) %>% 
  tidylog::summarise(met_name = str_c(met_name, 
                                      collapse = "; "), 
                     matched_compound = str_c(matched_compound,
                                              collapse = "; ")) %>% 
  select(empirical_compound, met_name, everything())


# Save final Key
write_rds(ecd_pw_key_final,
          fs::path(
            dir_data_local,
            "Supporting Files",
            "mummichog_pw_ec_feature_key_with_cpd_names.rds"))

