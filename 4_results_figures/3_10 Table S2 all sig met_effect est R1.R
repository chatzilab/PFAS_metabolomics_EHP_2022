# Table S2: All effect estimates 
rm(list = ls())
source(here::here("0_project_setup", "!libraries.R"))
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))


# aromatic AAs
aro_aa_met <- read_rds(
  fs::path(
    dir_results_mixtures, 
    "ind_met_effect_ests", 
    "aeromatic amino acid pooled effect estimates.rds")) %>% 
  bind_rows(.id = "cohort") %>% 
  select(-met_name) %>% 
  rename(met_name = met_name_tyr_pw) %>% 
  mutate(superpath = "Aromatic Amino Acid Metabolism", 
         pathway = "Tyrosine Metabolism")


nonaro_met <- read_rds(
  fs::path(dir_results_mixtures, 
           "ind_met_effect_ests", 
           "nonaeromatic aa effect estimates.rds")) %>% 
  bind_rows(.id = "cohort") %>% 
  mutate(superpath = "Non-aromatic Amino Acid Metabolism")


fa_met <- read_rds(
  fs::path(dir_results_mixtures, 
           "ind_met_effect_ests", 
           "lipid effect estimates.rds")) %>% 
  bind_rows(.id = "cohort") %>% 
  mutate(superpath = "Lipid Metabolism")


other_met <- read_rds(
  fs::path(dir_results_mixtures, 
           "ind_met_effect_ests", 
           "other_met_pathway effect estimates.rds"))  %>% 
  bind_rows(.id = "cohort") %>% 
  mutate(superpath = "Other met. pathways")


# Join all
# Bind rows 
all_ee <- bind_rows(aro_aa_met, nonaro_met, fa_met, other_met)

# four duplicates, but they are all the same metabolite
dup_mzrt(all_ee)
length(unique(all_ee$name))

# Select Key columns
sol_all <- all_ee %>% 
  tidylog::filter(bci_sig_sol == "Sig Sol", 
                  cohort == "solar") %>% 
  select(cohort, superpath, pathway, met_name, 
         matched_compound, mode, matched_form, query_mass,
         beta_ci_sol_mixture, 
         p_value_sol_mixture,
         q_value_sol_mixture) %>% 
  rename_all(~str_remove(., "_sol"))


chs_all <- all_ee %>% 
  tidylog::filter(bci_sig_chs == "Sig CHS", 
                  cohort == "chs") %>% 
  select(cohort, superpath, pathway, met_name, 
         matched_compound, mode, matched_form, query_mass,
         beta_ci_chs_mixture, 
         p_value_chs_mixture, 
         q_value_chs_mixture) %>% 
  rename_all(~str_remove(., "_chs"))

pooled_all <- all_ee %>% 
  tidylog::filter(bci_sig_pooled == "Sig Pooled", 
                  cohort == "pooled") %>% 
  select(cohort, superpath, pathway, met_name, 
         matched_compound, mode, matched_form, query_mass,
         beta_ci_pooled_mixture,
         p_value_pooled_mixture,
         q_value_pooled_mixture)  %>% 
  rename_all(~str_remove(., "_pooled"))


# Bind back together, count number of matched compounds for each mzrt
all_ee2 <- bind_rows(sol_all, chs_all, pooled_all) %>% 
  mutate(query_mass = as.character(query_mass), 
         across(where(is.numeric), 
                ~formatC(., digits = 2) %>% 
                  str_replace(" 1", ">0.99")), 
         n_met_5ppm = str_count(matched_compound, ";")+1) %>% 
  select(-matched_compound, -mode, -matched_form, -query_mass, -n_met_5ppm) %>% 
  arrange(cohort, superpath, pathway, met_name)


colnames(all_ee2) <- c("Cohort", "Superpathway", 
                       "Pathway", "Metabolite Name",
                       "Psi (95% BCI)", "P-value", "Q-value")
# "LC-MS Analytical Mode"
# Dimensions should be 66 rows by 10 coulmns
dim(all_ee2)
length(unique(all_ee2$`Metabolite Name`))


write_csv(all_ee2 %>% filter(Cohort != "pooled"), 
          fs::path(dir_reports, "Table S2 Metabolite.csv"))

write_csv(data.frame(met_name = unique(all_ee2$`Metabolite Name`)), 
          fs::path(dir_reports, "All Unique Metabolite Names.csv"))


# Calculate number of unique metabolites in each cohort

all_ee2 %>% 
  group_by(Cohort, Superpathway) %>% 
  summarise(n_mets = length(unique(`Metabolite Name`))) %>% 
  filter(Cohort != "pooled")


# Pivot wider 
all_ee2_wide <- all_ee2 %>% 
  janitor::clean_names() %>% 
  pivot_wider(id_cols = c(superpathway,  pathway, metabolite_name),
              values_from = c("psi_95_percent_bci", "p_value","q_value"), 
              names_from = cohort) %>% 
  mutate(nsig_all = is.na(p_value_solar) + 
           is.na(p_value_chs) + 
           is.na(p_value_pooled),
         nsig_cohorts = is.na(p_value_solar) + 
           is.na(p_value_chs)) #%>% 



table(all_ee2_wide$nsig_cohorts)  
length(unique(all_ee$met_name))




# Look at the PIPs ----------------------
## SOLAR --------------------
sol_all_pips <- all_ee %>% 
  tidylog::filter(bci_sig_sol == "Sig Sol", 
                  cohort == "solar") %>% 
  select(met_name, superpath, contains("estimate_pip_sol")) %>% 
  pivot_longer(cols = contains("pip"), values_to = "pip") %>% 
  mutate(pfas = str_remove(name, "estimate_pip_sol_"), 
         pip = as.numeric(pip))

# Minimum and maximum PIPs
summary(sol_all_pips$pip)

# Percent of ALL PIPs > 1/6
table(sol_all_pips$pip>1/6)/length(sol_all_pips$pip)

# Number of metabolites in each pathway (needed later)
nmet_superpath_sol <-  sol_all_pips  %>% 
  group_by(superpath) %>% 
  summarise(nmet_superpath = length(unique(met_name))) %>%
  ungroup()


# Percent of PIPs >1/6 for each metabolite
sol_pip_summary_by_met <- sol_all_pips  %>% 
  group_by(superpath, met_name) %>% 
  tidylog::summarise(n_over_0.16_metabolite = sum(pip>(1/6)), 
            pct_over_0.16_metabolite = sum(pip>(1/6))/length(pip)) %>% 
  ungroup() 

table(sol_pip_summary_by_met$n_over_0.16_metabolite >1)/nrow(sol_pip_summary_by_met)

# Percent of PIPs >1/6 for each PFAS, overall
sol_pip_summary_by_pfas <- sol_all_pips  %>% 
  group_by(pfas) %>% 
  tidylog::summarise(n_over_0.16_pfas = sum(pip>(1/6)), 
                     pct_over_0.16_pfas = sum(pip>(1/6))/length(pip)) %>% 
  ungroup() %>% 
  arrange(-pct_over_0.16_pfas)

# Percent of PIPs >1/6 for each PFAS, by pathway
sol_pip_summary_by_pfas <- sol_all_pips  %>% 
  group_by(superpath, pfas) %>% 
  tidylog::summarise(n_over_0.16_pfas = sum(pip>(1/6)), 
                     pct_over_0.16_pfas = sum(pip>(1/6))/length(pip)) %>% 
  ungroup() %>% 
  arrange(superpath, -pct_over_0.16_pfas)

table(sol_pip_summary_by_pfas$n_over_0.16_pfas >1)/nrow(sol_pip_summary_by_pfas)



# Dataset for Figure
sol_pip_summary_total <- sol_pip_summary_by_met %>% 
  group_by(superpath, n_over_0.16_metabolite) %>% 
  tidylog::summarise(n_in_superpath=length(n_over_0.16_metabolite)) %>% 
  ungroup() %>%
  tidylog::left_join(nmet_superpath_sol) %>% 
  mutate(pct = n_in_superpath/nmet_superpath) %>% 
  select(-n_in_superpath, -nmet_superpath) 
  
# Summary of PIPs
sol_pip_summary_total %>% 
  tidylog::pivot_wider(names_from = n_over_0.16_metabolite, 
                       values_from = pct) %>% 
  janitor::clean_names() %>% 
  select(superpath, x1, x2, x3, x4, x5, x6) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)),
                pct_mets_with_2_pfas_over_.16 = sum(c(x2, x3, x4, x5, x6)))


# Plot
ggplot(sol_pip_summary_total, 
       aes(y = pct, x = n_over_0.16_metabolite)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~superpath, 
             ncol = 1)




## CHS ---------------------
chs_all_pips <- all_ee %>% 
  tidylog::filter(bci_sig_chs == "Sig CHS", 
                  cohort == "chs") %>%
  select(met_name, superpath, contains("estimate_pip_chs")) %>% 
  pivot_longer(cols = contains("pip"), values_to = "pip") %>% 
  mutate(pfas = str_remove(name, "estimate_pip_chs_"), 
         pip = as.numeric(pip))


# Minimum and maximum PIPs
summary(chs_all_pips$pip)

# Percent of ALL PIPs > 1/6
table(chs_all_pips$pip>1/6)/length(chs_all_pips$pip)

# Number of metabolites in each pathway (needed later)
nmet_superpath_chs <-  chs_all_pips  %>% 
  group_by(superpath) %>% 
  summarise(nmet_superpath = length(unique(met_name))) %>%
  ungroup()


# Percent of PIPs >1/6 for each metabolite
chs_pip_summary_by_met <- chs_all_pips  %>% 
  group_by(superpath, met_name) %>% 
  tidylog::summarise(n_over_0.16_metabolite = sum(pip>(1/6)), 
                     pct_over_0.16_metabolite = sum(pip>(1/6))/length(pip)) %>% 
  ungroup() 

table(chs_pip_summary_by_met$n_over_0.16_metabolite >1)/nrow(chs_pip_summary_by_met)
sum(chs_pip_summary_by_met$n_over_0.16_metabolite >1)

# Percent of PIPs >1/6 for each PFAS, overall
chs_pip_summary_by_pfas <- chs_all_pips  %>% 
  group_by(pfas) %>% 
  tidylog::summarise(n_over_0.16_pfas = sum(pip>(1/6)), 
                     pct_over_0.16_pfas = sum(pip>(1/6))/length(pip)) %>% 
  ungroup() %>% 
  arrange(-pct_over_0.16_pfas)

# Percent of PIPs >1/6 for each PFAS, by pathway
chs_pip_summary_by_pfas <- chs_all_pips  %>% 
  group_by(superpath, pfas) %>% 
  tidylog::summarise(n_over_0.16_pfas = sum(pip>(1/6)), 
                     pct_over_0.16_pfas = sum(pip>(1/6))/length(pip)) %>% 
  ungroup() %>% 
  arrange(superpath, -pct_over_0.16_pfas)


table(chs_pip_summary_by_pfas$n_over_0.16_metabolite >1)/nrow(chs_pip_summary_by_pfas)



# Dataset for Figure
chs_pip_summary_total <- chs_pip_summary_by_met %>% 
  group_by(superpath, n_over_0.16_metabolite) %>% 
  tidylog::summarise(n_in_superpath=length(n_over_0.16_metabolite)) %>% 
  ungroup() %>%
  tidylog::left_join(nmet_superpath_chs) %>% 
  mutate(pct = n_in_superpath/nmet_superpath) %>% 
  select(-n_in_superpath, -nmet_superpath) 

# Summary of PIPs
chs_pip_summary_total %>% 
  tidylog::pivot_wider(names_from = n_over_0.16_metabolite, 
                       values_from = pct) %>% 
  janitor::clean_names() %>% 
  select(superpath, x1, x2, x3, x4, x5, x6) %>% 
  rowwise() %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)),
         pct_mets_with_2_pfas_over_.16 = sum(c(x2, x3, x4, x5, x6)))


# Plot
ggplot(chs_pip_summary_total, 
       aes(y = pct, x = n_over_0.16_metabolite)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~superpath, 
             ncol = 1)

