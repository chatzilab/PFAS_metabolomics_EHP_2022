# plot PIPs for figure S3 (R1)
# have to run  prior to this script
rm(list = ls())
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))



# Read data from script 3_3
nonaro_met <- read_rds(
  fs::path(dir_results_mixtures, 
           "ind_met_effect_ests", 
           "nonaeromatic aa effect estimates.rds"))

# Subset SOLAR and CHS
sol <- nonaro_met$solar %>% 
  tidylog::select(met_name:mass_diff, pathway_reduced, contains("sol")) %>% 
  mutate(cohort = "SOLAR") %>% 
  tidylog::rename_all(~str_remove(., "_sol"))


chs <- nonaro_met$chs %>% 
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
nonaro_mets <- sig_ecs_l_3 %>% 
  pivot_wider(id_cols = c(cohort:mass_diff, exposure, pathway_reduced), 
              names_from = effect, 
              values_from = value)  %>% 
  select(cohort, exposure, everything()) %>% 
  mutate(estimate_b = estimate_beta*estimate_pip, 
         sig_beta = if_else(lcl_beta >= 0 |
                              ucl_beta <= 0, estimate_beta, 0), 
         pfas = rename_pfas(exposure,arrange_by_class = TRUE))


sol_data <- nonaro_mets %>% 
  filter(cohort == "SOLAR", 
         exposure != "mixture")

chs_data <-  nonaro_mets %>% 
  filter(cohort == "CHS", 
         exposure != "mixture")


# Plot -----------------------------------------
## SOLAR --------
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
         axis.text.x = element_blank(), 
         strip.background = element_rect(fill = "white"),
         legend.position = "none",
         strip.text.y = element_text(angle = 0, hjust = 0)) +
   facet_grid(pathway_reduced~cohort, scales = "free",
              space = "free") +
   scale_color_manual(values = c("black", "white"), 
                      guide = NULL) +
   scale_fill_viridis_c(option = "plasma", 
                        breaks = seq(from = 1/6, to = 5/6,by = 1/6),
                        labels = c("0.17", "",
                                   "0.5", "",
                                   "0.83"))) 


## CHS --------------
(chs <-  chs_data %>% 
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
         axis.text.x = element_text(angle = 90, 
                                    hjust = 1,
                                    vjust = 0.5), 
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


# Combine Plots -----------------------
(non_aa_met_pips <- plot_grid(NULL, sol, 
                           NULL, chs, 
                           ncol = 1, align = "v", 
                           rel_heights = c(.05, 1,.05, .75),
                           labels = c("A) SOLAR","",
                                      "B) CHS", "")))


ggsave(non_aa_met_pips,
       filename = fs::path(
         dir_reports,
         "Figure S3 PFAS and nonaro aa metabolite pips_R1.tiff"),
       width = 9, height = 9)


# Summary of PIPs ------------------------------
summary(sol_data$estimate_pip)
quantile(sol_data$estimate_pip, 0.9, na.rm = TRUE)


sol_data %>% 
  group_by(pfas) %>% 
  summarise(n_pip_above_50 = sum(estimate_pip>.88)/length(estimate_pip))



length(unique(chs_data$met_name))
summary(chs_data$estimate_pip)
quantile(chs_data$estimate_pip, 0.9, na.rm = TRUE)


chs_data %>% 
  group_by(pfas) %>% 
  summarise(n_pip_above_50 = sum(estimate_pip>.632)/length(estimate_pip))
