# Make coefficient plots based on dougs comments
rm(list = ls())
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))

# Read data
annotated_sig_ecs_ee <- read_rds(
  here::here(
    dir_results_mixtures,
    "Sig annotated metabolite effect estimates sol chs_v4.RDS") ) %>% # New
  select(-contains("blank")) %>% 
  tidylog::filter(str_detect(super_pathway, "Amino acid metabolism"),
                  str_detect(super_pathway, "Aromatic",
                             negate = TRUE))

# These are unique from the metabolites in the tyrosine pathway
length(unique(annotated_sig_ecs_ee$met_name))
length(unique(annotated_sig_ecs_ee$empirical_compound))
length(unique(annotated_sig_ecs_ee$name))
length(unique(annotated_sig_ecs_ee$name))

# Rename metabolites --------------
annotated_sig_ecs_ee_aa <- annotated_sig_ecs_ee %>% 
  tidylog::mutate(
    met_name = case_when(met_name == "L-Aspartic acid; D-Aspartic acid" ~ "Aspartic acid",
                         met_name == "CE4788; L-Proline; D-Proline" ~ "Proline", 
                         str_detect(met_name, 
                                    "5-Amino-2-oxopentanoic acid") ~  
                           "5-Amino-2-oxopentanoic acid", 
                         str_detect(met_name, 
                                    "2,3,4,5-Tetrahydro-2-pyridinecarboxylic acid") ~ 
                           "2,3,4,5-Tetrahydro-2-pyridinecarboxylic acid",
                         str_detect(met_name, 
                                    "1-Pyrroline-4-hydroxy-2-carboxylate") ~ 
                           "1-Pyrroline-4-hydroxy-2-carboxylate",
                         str_detect(met_name, 
                                    "2-Keto-6-aminocaproate") ~ 
                           "6-Amino-2-oxohexanoate", 
                         TRUE ~ met_name))

# Filter out duplicate mz/rts
annotated_sig_ecs_ee_aa_reduced <- annotated_sig_ecs_ee_aa %>% 
  tidylog::filter(
    # for 160.061586644875_256.162228074228 Choose M-H[-] matched form: 
    !(name == "160.061586644875_256.162228074228" & matched_form == "M+HCOO[-]"), 
    # for 142.051011811677_41.3053679062358 Choose M-H[-] matched form: 
    !(name == "142.051011811677_41.3053679062358" & matched_form == "M-H+O[-]"),
    # for 131.118018715004_272.694680893768 choose M+H[1+]
    !(name == "131.118018715004_272.694680893768" & matched_form == "M-CO2+H[1+]"),
   # for 175.107820519341_268.639087046289, choose M+H[1+]
   !(name == "175.107820519341_268.639087046289" & matched_form == "M-CO2+H[1+]"), 
   # for 181.035512237065_26.5065490174735, choose M+H[1+]
   !(name == "181.035512237065_26.5065490174735" & matched_form == "M+Na-2H[-]"),
   
   # This mz/rt is annotated in the aromatic amino acid metabolism pathway
   !(met_name == "1-Pyrroline-4-hydroxy-2-carboxylate" & matched_form == "M-H[-]"), 
   # This mz/rt is annotated in the aromatic amino acid metabolism pathway
   !(met_name == "Oxoadipic acid" & matched_form == "M-H+O[-]"), 
   # This mz/rt is annotated in the aromatic amino acid metabolism pathway
   !(met_name == "L-Carnitine" & matched_form == "M+K[1+]"), 
   
   # This is a weird adduct: I would not tryst
   matched_form != "M+Cl37[-]"
   
   )



dup_mzrt(annotated_sig_ecs_ee_aa_reduced)

length(unique(annotated_sig_ecs_ee_aa$name))/length(annotated_sig_ecs_ee_aa$name)


# Reorder tyrosine sub path by number of cpds within path
sig_ecs_reduced <- annotated_sig_ecs_ee_aa_reduced %>% 
  mutate(pathway_reduced = case_when(
           str_detect(pathway, "Lysine") ~ "Lysine Metabolism", 
           str_detect(pathway, "Glutathione") ~ "Glutathione Metabolism", 
           str_detect(pathway, "Arginine") ~ "Arginine and Proline\nMetabolism", 
           str_detect(pathway, "Urea cycle") ~ "Urea cycle",
           TRUE ~ pathway),
         qsig_sol = if_else(q_value_sol_mixture<0.05, 
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

## Subset SOL (filter sig features) ----------------
sol <- sig_ecs_reduced %>% 
  select(met_name:mass_diff,empirical_compound, pathway_reduced, everything()) %>% 
  tidylog::filter(ucl_beta_sol_mixture <= 0 | 
                    lcl_beta_sol_mixture >= 0) %>% 
  group_by(met_name) %>%
  arrange(-abs(estimate_beta_sol_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name = fct_reorder(met_name, 
                                estimate_beta_sol_mixture), 
         borderline_sig = if_else(lcl_beta_sol_mixture == 0 | 
                                    ucl_beta_sol_mixture == 0, 
                                  "borderline", 
                                  "sig"))

## Subset CHS (filter sig features) ----------------
chs <- sig_ecs_reduced %>% 
  select(met_name:mass_diff,empirical_compound,pathway_reduced, everything()) %>% 
  tidylog::filter(ucl_beta_chs_mixture <= 0 | 
                    lcl_beta_chs_mixture >= 0) %>% 
  group_by(met_name) %>% 
  arrange(-abs(estimate_beta_chs_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name = fct_reorder(met_name, 
                                estimate_beta_chs_mixture))


## Subset Pooled (filter sig features) ----------------
pooled <- sig_ecs_reduced %>% 
  select(met_name:mass_diff,empirical_compound,pathway_reduced,
         contains("pooled"), everything()) %>% 
  tidylog::filter(ucl_beta_pooled_mixture <= 0 | 
                    lcl_beta_pooled_mixture >= 0) %>% 
  group_by(met_name) %>% 
  arrange(-abs(estimate_beta_pooled_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name = fct_reorder(met_name, 
                                estimate_beta_pooled_mixture), 
         borderline_sig = if_else(lcl_beta_pooled_mixture == 0 | 
                                    ucl_beta_pooled_mixture == 0, 
                                  "borderline", 
                                  "sig"))

## Get data on overlap ---------------------
sol_chs_overlap <- tidylog::full_join(sol %>% 
                                        select(met_name, 
                                               estimate_beta_sol_mixture), 
                                      chs %>% 
                                        select(met_name, 
                                               estimate_beta_chs_mixture)) %>% 
  mutate(sig_dir = estimate_beta_sol_mixture*estimate_beta_chs_mixture, 
         sig_both = !is.na(sig_dir), 
         sig_with_direction = case_when(
           sig_both & sig_dir < 0 ~  
             "Opposite direction of association in SOLAR vs. CHS",
           sig_both & sig_dir > 0 ~ 
             "Same direction of association in SOLAR and CHS", 
           TRUE ~  "Only significant in one cohort") %>% 
           fct_relevel("Same direction of association in SOLAR and CHS", 
                       "Only significant in one cohort",
                       "Opposite direction of association in SOLAR vs. CHS")) %>% 
  select(-estimate_beta_chs_mixture, -estimate_beta_sol_mixture)

# Join with solar/chs data
sol_nonaro <- sol %>% 
  tidylog::left_join(sol_chs_overlap) %>% 
  mutate(qsig = if_else(q_value_pooled_mixture <0.05, "*", ""),
         met_name_sig = str_c(met_name, qsig))
length(unique(sol_nonaro$name))/length(sol_nonaro$name)

chs_nonaro <- chs %>% 
  tidylog::left_join(sol_chs_overlap) %>% 
  mutate(qsig = if_else(q_value_pooled_mixture <0.05, "*", ""),
         met_name_sig = str_c(met_name, qsig))
length(unique(chs$name))/length(chs$name)

pooled_nonaro <- pooled %>% 
  tidylog::left_join(sol_chs_overlap) %>% 
  mutate(qsig = if_else(q_value_pooled_mixture <0.05, "*", ""),
         met_name_sig = str_c(met_name, qsig))
length(unique(chs$name))/length(chs$name)


# Make List 
nonaro_met <- list(solar = sol_nonaro, 
                   chs = chs_nonaro,
                   pooled = pooled_nonaro)


# Save data 
write_rds(nonaro_met,
          fs::path(dir_results_mixtures, 
                   "ind_met_effect_ests", 
                   "nonaeromatic aa effect estimates.rds"))


# Make Plots ----------------------------------------------------------
length(unique(sol$name))

## Plot solar------------------------------
(sol_met_plot <- ggplot(sol_nonaro, 
                        aes(x = met_name, 
                            y = estimate_beta_sol_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_sol_mixture , 
                     ymax = ucl_beta_sol_mixture), 
                 width = 0, size = .75) + 
   geom_point(aes(shape = sig_with_direction, 
                  fill = sig_with_direction, 
                  size = sig_with_direction)) + 
   geom_text(aes(label = qsig_sol, y = 2.1),
             size = 3.5, hjust = 0, vjust = .5) + 
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   facet_grid(pathway_reduced ~ .,  scales = "free_y", space = "free_y") +
   # Scales
   scale_fill_manual(values = c( "black", "grey50", "white"), name = NULL) +
   scale_shape_manual(values = c(21, 21, 22), name = NULL) +
   scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
   scale_y_continuous(limits = c(-1.75, 3)) +
   coord_flip(clip = "off") +    
   # ylab("PFAS Exposure Effect Estimate ψ (95% BCI)") +
   theme(axis.title.y = element_blank(), 
         axis.title.x = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey95"), 
         strip.background = element_rect(fill = "white"),
         legend.position = "none",
         legend.direction = "vertical",
         strip.text.y = element_text(angle = 0, hjust = 0))) 

## Plot CHS ------------------------------
(chs_metplot <- ggplot(chs_nonaro, 
                       aes(x = met_name, 
                           y = estimate_beta_chs_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_chs_mixture , 
                     ymax = ucl_beta_chs_mixture, 
                     linetype = sig_with_direction), 
                 width = 0, size = .75) + 
   geom_point(aes(shape = sig_with_direction, 
                  fill = sig_with_direction, 
                  size = sig_with_direction)) + 
   geom_text(aes(label = qsig_chs, y = 2.1),
             size = 3.5, hjust = 0, vjust = .5) + 
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   facet_grid(pathway_reduced~"1",  scales = "free_y", space = "free_y") +
   # Scales
   scale_fill_manual(values = c( "black", "grey50", "white"), name = NULL) +
   scale_shape_manual(values = c(21, 21, 22), name = NULL) +
   scale_linetype_manual(values = c(1,1,1), name = NULL) +
   scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
   scale_y_continuous(limits = c(-1.75, 3)) +
   # Other
   coord_flip() + 
   ylab("PFAS Mixture Effect ψ (95% BCI)") +
   theme(axis.title.y = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey95"), 
         strip.background = element_rect(fill = "white"),
         legend.position = "bottom",
         legend.direction = "vertical",
         legend.justification = c("center","bottom"),
         strip.text.y = element_text(angle = 0, hjust = 0)))


# Combine Plots
(figure_aa_met <- plot_grid(NULL, sol_met_plot, 
                            NULL, chs_metplot, 
                            ncol = 1, align = "v", 
                            rel_heights = c(.05, 1,.05, .85),
                            labels = c("A) SOLAR","",
                                       "B) CHS", "")))

## Save ----------------
ggsave(figure_aa_met,
       filename = fs::path(dir_reports,
                           "Figure 3 R1 associations of PFAS and nonaro AA metabolites.jpg"),
       width = 7.5, height = 7)


# Plot pooled------------------------------
# (pooled_nonaro_plot <- ggplot(pooled_nonaro, 
#                         aes(x = met_name, 
#                             y = estimate_beta_pooled_mixture)) + 
#    geom_errorbar(aes(ymin = lcl_beta_pooled_mixture , 
#                      ymax = ucl_beta_pooled_mixture), 
#                  width = 0, size = .75) + 
#    geom_point(aes(shape = sig_with_direction, 
#                   fill = sig_with_direction, 
#                   size = sig_with_direction)) + 
#    geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
#    facet_grid(pathway_reduced ~ .,  scales = "free_y", space = "free_y") +
#    # Scales
#    scale_fill_manual(values = c( "black", "grey50", "white"), name = NULL) +
#    scale_shape_manual(values = c(21, 21, 22), name = NULL) +
#    scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
#    scale_y_continuous(limits = c(-1.75, 3)) +
#    coord_flip(clip = "off") +    
#    ylab("PFAS Exposure Effect Estimate ψ (95% BCI)") +
#    theme(axis.title.y = element_blank(), 
#          strip.text.x = element_blank(), 
#          panel.background = element_rect(fill="grey95"), 
#          strip.background = element_rect(fill = "white"),
#          legend.position = "bottom",
#          legend.direction = "vertical",
#          legend.justification = c("center","bottom"),
#          strip.text.y = element_text(angle = 0, hjust = 0)) ) 



