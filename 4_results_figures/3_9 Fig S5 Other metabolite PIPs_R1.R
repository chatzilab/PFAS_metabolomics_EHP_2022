# plot PIPs for lipid metabolism
rm(list = ls())
source(here::here("0_project_setup", "!libraries.R"))
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))

# Read in Data
other_met <- read_rds(
  fs::path(dir_results_mixtures, 
           "ind_met_effect_ests", 
           "other_met_pathway effect estimates.rds"))

# Subset SOLAR and CHS
sol <- other_met$solar %>% 
  tidylog::select(met_name:mass_diff, pathway_reduced, contains("sol")) %>% 
  mutate(cohort = "SOLAR") %>% 
  tidylog::rename_all(~str_remove(., "_sol"))


chs <- other_met$chs %>% 
  select(met_name:mass_diff,pathway_reduced, contains("chs")) %>% 
  mutate(cohort = "CHS") %>% 
  tidylog::rename_all(~str_remove(., "_chs"))


# Combine
sig_ecs_l <- bind_rows(sol, chs) %>% 
  tidylog::select(cohort, everything())


sig_ecs_l_2 <- sig_ecs_l %>% 
  select(-contains("beta_ci_"), -contains("neg_log_p"),
         -contains("sd_"), -contains("wald_")) %>% 
  mutate(across(estimate_beta_mixture:q_value_pfos, as.numeric)) %>%
  pivot_longer(cols = estimate_beta_mixture:q_value_pfos, 
               names_to = "effect_name")

# separate effect and term
sig_ecs_l_3 <- sig_ecs_l_2 %>% 
  mutate(effect_1 = str_split_fixed(effect_name, "_", 3)[,1], 
         effect_2 = str_split_fixed(effect_name, "_", 3)[,2],
         exposure = str_split_fixed(effect_name, "_", 3)[,3], 
         effect = str_c(effect_1, effect_2, sep = "_")) %>% 
  select(-effect_1, -effect_2, effect_name)

# Pivot wider
fa_mets <- sig_ecs_l_3 %>% 
  pivot_wider(id_cols = c(cohort:mass_diff, exposure, pathway_reduced), 
              names_from = effect, 
              values_from = value)  %>% 
  select(cohort, exposure, everything()) %>% 
  mutate(estimate_b = estimate_beta*estimate_pip, 
         sig_beta = if_else(lcl_beta >= 0 |
                              ucl_beta <= 0, estimate_beta, 0), 
         pfas = rename_pfas(exposure,arrange_by_class = TRUE))


sol_data <- fa_mets %>% 
  tidylog::filter(cohort == "SOLAR", 
                  exposure != "mixture")

chs_data <-  fa_mets %>% 
  tidylog::filter(cohort == "CHS", 
                  exposure != "mixture")

# Plot -----------------------------------------
(sol <- sol_data %>% 
   mutate(pip_over_50 = if_else(estimate_pip > 1/6, 
                                formatC(estimate_pip,
                                        digits = 2,
                                        format = 'g', 
                                        flag = "#") , 
                                ""), 
          color_pips = if_else(estimate_pip>0.7, "1","2")) %>% 
   ggplot(aes(x = pfas,
              y = met_name,
              fill = estimate_pip)) + 
   geom_tile() +
   geom_text(aes(label = pip_over_50, color = color_pips)) +
   theme(axis.title.y = element_blank(), 
         axis.title.x = element_blank(), 
         strip.text.x = element_blank(), 
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
         strip.background = element_rect(fill = "white"),
         legend.position = "bottom",
         legend.justification = c("center","bottom"),
         strip.text.y = element_text(angle = 0, hjust = 0), 
         legend.key.size = unit(.6, "cm")) +
   facet_grid(pathway_reduced~cohort, scales = "free",
              space = "free") +
   scale_color_manual(values = c("black", "white"), 
                      guide = NULL) +
   scale_fill_viridis_c(option = "plasma", 
                        breaks = seq(from = 1/6, to = 5/6,by = 1/6),
                        labels = c("0.17", "",
                                   "0.5", "",
                                   "0.83"), 
                        name = "Posterior inclusion probability"))  

dim(chs_data)


ggsave(sol,
       filename = fs::path(dir_reports,
                           "Figure S5 SOLAR other metabolite pips_R1.jpg"),
       width = 7, height = 2.5)


# Summary of PIPs ------------------------------
quantile(sol_data$estimate_pip, 0.9, na.rm = TRUE)

sol_data %>% 
  group_by(pfas) %>% 
  summarise(n_pip_above_50 = sum(estimate_pip>.88)/length(estimate_pip))



summary(chs_data$estimate_pip)
quantile(chs_data$estimate_pip, 0.9, na.rm = TRUE)


chs_data %>% 
  group_by(pfas) %>% 
  summarise(n_pip_above_50 = sum(estimate_pip>.63)/length(estimate_pip))