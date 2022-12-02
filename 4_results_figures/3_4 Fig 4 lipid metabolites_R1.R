# Make coefficient plots for PFAS and FFA metabolites
rm(list = ls())
source(here::here("0_project_setup", "!libraries.R"))
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))

# Read data
annotated_sig_ecs_ee <- read_rds(
  here::here(dir_results_mixtures,
             "Sig annotated metabolite effect estimates sol chs_v4.RDS") ) %>% 
  select(-contains("blank")) %>% 
  tidylog::filter(super_pathway == "Lipid metabolism")

table(annotated_sig_ecs_ee$super_pathway)

# Change metabolite names based on lipid metabolism --------------
annotated_sig_ecs_ee_lipid <- annotated_sig_ecs_ee %>% 
  tidylog::mutate(met_name = case_when(
    met_name == "Oleic acid; Elaidic acid" ~ "Elaidic acid", 
    met_name == "CE5707" ~ "11-hydroxyeicosatetraenoate glyceryl ester", # from metabolicatlas.org
    met_name == "CE5528" ~ "12,13-epoxy-9-alkoxy-10E-octadecenoate", # from metabolicatlas.org
    str_detect(met_name, "13\\(S\\)-HPOT") ~ "13(S)-HPOT",
    # str_detect(met_name, "Dodecanoic acid") ~ "Prostaglandin G2",
    str_detect(met_name, "Prostaglandin G2") ~ "Prostaglandin G2",
    str_detect(met_name, "S-\\(PGA2\\)-glutathione") ~ "S-(PGA2)-glutathione",
    str_detect(met_name, "Prostaglandin A2") ~ "Prostaglandin E2",  # These are all the same metabolite, noted in discussion
    str_detect(met_name, "Prostaglandin H2") ~ "Prostaglandin E2", # These are all the same metabolite, noted in discussion
    str_detect(met_name, "13-OxoODE") ~ "13-OxoODE",
    str_detect(met_name, "15-Keto-prostaglandin E2") ~ "15-Keto-prostaglandin E2",
    str_detect(met_name, "9\\(S\\)-HPODE") ~ "9(S)-HPODE", 
    str_detect(met_name, "\\(E\\)-4-hydroxynon-2-enal") ~ "(E)-4-hydroxynon-2-enal", 
    TRUE ~ met_name), 
    pathway_reduced = 
      case_when(str_detect(pathway, "Prostaglandin formation from arachidonate") ~ 
                  "Prostaglandin formation\nfrom arachidonate",
                str_detect(pathway, "Linoleate") ~ "Linoleate metabolism", 
                str_detect(pathway, "De novo fatty acid") ~ "De novo fatty acid biosynthesis", 
                str_detect(pathway, "Putative anti") ~ 
                  "Putative anti-Inflammatory\nmetabolites formation from EPA",
                TRUE ~ pathway))



# Remove duplicate features
annotated_sig_ecs_ee_lipid <- annotated_sig_ecs_ee_lipid %>% 
  tidylog::filter(
    
    # 9(S)-HPODE and (E)-4-hydroxynon-2-enal are the same mz/rt": use M+H[1+] matched form
    !(met_name == "9(S)-HPODE" & matched_form == "M+2H[2+]"),
    
    # 13-OxoODE and 9(S)-HPODE are same mzrt: choose M-H2O+H[1+] matched form
    !(name == "277.216423836957_65.2041530713161" & matched_form == "M-H4O2+H[1+]"), 
    # Chose the M+H[1+] matched form:
    !(name == "335.2218739_181.3322705" & matched_form == "M-H2O+H[1+]"),
    # Chose the M+H[1+] matched form:
    !(name == "351.215567931475_66.3057204318113" & matched_form == "M-H2O+H[1+]"),
    # Chose the M-H2O+H[1+] matched form:
    !(name == "604.2696107_315.0790795" & matched_form == "M-H4O2+H[1+]") )



# Create sig. text for fts qith q<0.05
annotated_sig_ecs_ee_lipid <- annotated_sig_ecs_ee_lipid %>% 
  mutate(qsig_sol = if_else(q_value_sol_mixture<0.05, 
                            paste0("q=", 
                                   formatC(q_value_sol_mixture, digits = 1)), 
                            ""), 
         qsig_chs = if_else(q_value_chs_mixture<0.05, 
                            paste0("q=", 
                                   formatC(q_value_chs_mixture, digits = 1)), 
                            ""), 
         qsig_pooled = if_else(q_value_pooled_mixture<0.05, 
                               paste0("q=", 
                                      formatC(q_value_pooled_mixture, digits = 1)), 
                               ""))

# Subset SOL and CHS ---------------------------------------
sol <- annotated_sig_ecs_ee_lipid %>% 
  select(met_name, pathway:mass_diff, 
         contains("estimate"),
         contains("beta"),
         everything()) %>%
  tidylog::filter(ucl_beta_sol_mixture <= 0 | 
                    lcl_beta_sol_mixture >= 0) %>% 
  group_by(met_name) %>% 
  arrange(-abs(estimate_beta_sol_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name = fct_reorder(met_name, 
                                estimate_beta_sol_mixture))



chs <- annotated_sig_ecs_ee_lipid %>% 
  select(met_name:mass_diff,empirical_compound,pathway, everything()) %>% 
  tidylog::filter(ucl_beta_chs_mixture <= 0 | 
                    lcl_beta_chs_mixture >= 0) %>%
  # Group by metabolite name and select the metabolite with the largest abs(effect)
  group_by(met_name) %>% 
  arrange(-abs(estimate_beta_chs_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name = fct_reorder(met_name, 
                                estimate_beta_chs_mixture))


pooled <- annotated_sig_ecs_ee_lipid %>% 
  select(met_name:mass_diff,empirical_compound,pathway,
         contains("pooled"), everything()) %>% 
  tidylog::filter(ucl_beta_pooled_mixture <= 0 | 
                    lcl_beta_pooled_mixture >= 0) %>% 
  # Group by metabolite name and select the metabolite with the largest abs(effect)
  group_by(met_name) %>% 
  arrange(-abs(estimate_beta_pooled_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name = fct_reorder(met_name, 
                                estimate_beta_pooled_mixture))

# Get overlap ---------------------------------
sol_chs_overlap <- tidylog::full_join(sol %>% 
                                        select(met_name, 
                                               estimate_beta_sol_mixture), 
                                      chs %>% 
                                        select(met_name, 
                                               estimate_beta_chs_mixture)) %>% 
  mutate(sig_dir = estimate_beta_sol_mixture*estimate_beta_chs_mixture, 
         sig_both = !is.na(sig_dir), 
         sig_with_direction = case_when(sig_both & sig_dir < 0 ~  
                                          "Opposite direction of association in SOLAR vs. CHS",
                                        sig_both & sig_dir > 0 ~ 
                                          "Same direction of association in SOLAR and CHS", 
                                        TRUE ~  "Only significant in one cohort") %>% 
           fct_relevel("Same direction of association in SOLAR and CHS", 
                       "Only significant in one cohort",
                       "Opposite direction of association in SOLAR vs. CHS")) %>% 
  select(-estimate_beta_chs_mixture, -estimate_beta_sol_mixture)

# Get with sol/chs data
sol_fa <- sol %>% 
  tidylog::left_join(sol_chs_overlap)
length(unique(sol$name))/length(sol$name)


chs_fa <- chs %>% 
  tidylog::left_join(sol_chs_overlap)
length(unique(sol$name))/length(sol$name)


pooled_fa <- pooled %>% 
  tidylog::left_join(sol_chs_overlap)
length(unique(sol$name))/length(sol$name)


# Make List 
fa_met <- list(solar = sol_fa, 
               chs = chs_fa,
               pooled = pooled_fa)


# Save data 
write_rds(fa_met,
          fs::path(dir_results_mixtures, 
                   "ind_met_effect_ests", 
                   "lipid effect estimates.rds"))
# Make Plots ----------------------------------------------------------

## solar------------------------------
(sol_met_plot <- ggplot(sol_fa, aes(x = met_name, 
                                    y = estimate_beta_sol_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_sol_mixture , 
                     ymax = ucl_beta_sol_mixture, 
                     linetype = sig_with_direction), 
                 width = 0) +
   geom_point(aes(shape = sig_with_direction, 
                  fill = sig_with_direction, 
                  size = sig_with_direction)) + 
   # geom_text(aes(label = qsig_sol, y = 2.4),
   #           size = 3.5, hjust = 0, vjust = .5) + 
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   # Scales
   scale_fill_manual(values = c( "black", "grey50", "white"), name = NULL) +
   scale_shape_manual(values = c(21, 21, 22), name = NULL) +
   scale_linetype_manual(values = c(1,1,1), name = NULL) +
   scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
   scale_y_continuous(limits = c(-1.75, 3)) +
   coord_flip() + 
   facet_grid(pathway_reduced ~ .,  scales = "free_y", space = "free_y") +
   
   # ylab("PFAS Exposure Effect Estimate ψ (95% BCI)") +
   theme(axis.title.y = element_blank(), 
         axis.title.x = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey95"), 
         strip.background = element_rect(fill = "white"),
         legend.position = "none",
         strip.text.y = element_text(angle = 0, hjust = 0),
         text = element_text(size = 15),
         panel.grid = element_blank(),
         axis.line.x = element_line(color = "black"),
         axis.line.y = element_line(color = "black")))



## chs ------------------------------
(chs_metplot <- ggplot(chs_fa, aes(x = met_name, 
                                   y = estimate_beta_chs_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_chs_mixture , 
                     ymax = ucl_beta_chs_mixture, 
                     linetype = sig_with_direction), 
                 width = 0) + 
   geom_point(aes(shape = sig_with_direction, 
                  fill = sig_with_direction, 
                  size = sig_with_direction)) + 
   # geom_text(aes(label = qsig_chs, y = 2.4),
   #           size = 3.5, hjust = 0, vjust = .5) + 
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   # Scales
   scale_fill_manual(values = c( "black", "grey50", "white"), name = NULL) +
   scale_shape_manual(values = c(21, 21, 22), name = NULL) +
   scale_linetype_manual(values = c(1,1,1), name = NULL) +
   scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
   scale_y_continuous(limits = c(-1.75, 3)) +
   coord_flip() + 
   facet_grid(pathway_reduced ~ .,  scales = "free_y", space = "free_y") +
   ylab("PFAS Mixture Effect ψ (95% BCI)") +
   theme(axis.title.y = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey95"), 
         strip.background = element_rect(fill = "white"),
         legend.position = "bottom",
         legend.direction = "vertical",
         strip.text.y = element_text(angle = 0, hjust = 0),
         text = element_text(size = 15),
         panel.grid = element_blank(),
         axis.line.x = element_line(color = "black"),
         axis.line.y = element_line(color = "black")))


# Combine Plots
(coefplot_lipid_met <- plot_grid(NULL, sol_met_plot, 
                                 NULL, chs_metplot, 
                                 ncol = 1, align = "v", 
                                 rel_heights = c(.07, 1,.07, .6),
                                 labels = c("A) SOLAR","",
                                            "B) CHS", "")))

# Save 
ggsave(coefplot_lipid_met, 
       filename = fs::path(
         dir_reports, 
         "Figure 4 R1 associations of PFAS and lipid metabolites.jpg"), 
       width = 9.5, height = 7.5, bg = "white")
