# Make coefficient plots for PFAS and FFA metabolites
rm(list = ls())
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))

# Read data
annotated_sig_ecs_ee <- read_rds(
  here::here(dir_results_mixtures,
             "Sig annotated metabolite effect estimates sol chs_v4.RDS") ) %>% 
  select(-contains("blank")) %>% 
  tidylog::filter(
    str_detect(super_pathway, "Amino acid", negate = TRUE), 
    str_detect(super_pathway, "Lipid", negate = TRUE))

# Get duplicate mzrts 
dups <- dup_mzrt(annotated_sig_ecs_ee)
# two duplicate mzrts

# Filter out duplicate mzrts 
annotated_sig_ecs_ee <- annotated_sig_ecs_ee %>% 
  tidylog::filter(
    # For both dup features, select M-H[-] matched form
    !(name == dups$name[1] & matched_form == "M-H+O[-]"),
    !(name == dups$name[2] & matched_form == "M-H+O[-]"), 
    
    # This is also annotated as a metabolite in tyrosine metabolism (Noradrenochrome)
    name != "164.035423844049_70.3920835468602",
    
    
    # This is also annotated as hippuric acid in tyrosine metabolism
    !(name =="152.070757656544_74.5306937468636" & 
        matched_form == "M-NH3+H[1+]"), 
    
    
    # This is also annotated as ascorbate in tyrosine metabolism
    !(name =="175.024852523833_313.691984199347" & 
        matched_form == "M-H2O-H[-]"), 
    
    
    !(name =="194.081389732362_69.3433932711179" & 
        matched_form == "M+H[1+]"), 
    
    str_detect(pathway, "Drug metabolism", negate = TRUE)
  )
dups <- dup_mzrt(annotated_sig_ecs_ee)

# Change metabolite names based on other pathways --------------
annotated_sig_ecs_ee_other <- annotated_sig_ecs_ee %>% 
  tidylog::mutate(
    met_name = case_when(
      str_detect(met_name,"2-Hydroxyfelbamate") ~ "2-Hydroxyfelbamate", 
      str_detect(met_name,"3-Carbamoyl-2-phenylpropionaldehyde") ~ "3-Carbamoyl-2-phenylpropionaldehyde",
      str_detect(met_name,"2-Hydroxycarbamazepine") ~ "2-Hydroxycarbamazepine", 
      str_detect(met_name,"Glucuronic acid") ~ "Glucuronic acid", 
      TRUE ~ met_name), 
    
    pathway_reduced = 
      case_when(
        str_detect(pathway, "Drug metabolism - cytochrome P450") ~ 
          "Drug metabolism-\ncytochrome P450",
        str_detect(pathway, "Vitamin B6 \\(pyridoxine\\) metabolism") ~ 
          "Vitamin B6 (pyridoxine)\nmetabolism",
        TRUE ~ pathway)
  )

unique(annotated_sig_ecs_ee_other$pathway_reduced)

# Create sig. text for fts with q<0.05
annotated_sig_ecs_ee_other <- annotated_sig_ecs_ee_other %>% 
  mutate(qsig_sol = if_else(q_value_sol_mixture<0.05, 
                            paste0("q=", 
                                   formatC(q_value_sol_mixture, 
                                           digits = 1)), 
                            ""), 
         qsig_chs = if_else(q_value_chs_mixture<0.05, 
                            paste0("q=", 
                                   formatC(q_value_chs_mixture, 
                                           digits = 1)), 
                            ""), 
         qsig_pooled = if_else(q_value_pooled_mixture<0.05, 
                               paste0("q=", 
                                      formatC(q_value_pooled_mixture, 
                                              digits = 1)), 
                               ""))

# Subset SOL and CHS ---------------------------------------
sol <- annotated_sig_ecs_ee_other %>% 
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



chs <- annotated_sig_ecs_ee_other %>% 
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


pooled <- annotated_sig_ecs_ee_other %>% 
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
sol_other <- sol %>% 
  tidylog::left_join(sol_chs_overlap)
length(unique(sol$name))/length(sol$name)


chs_other <- chs %>% 
  tidylog::left_join(sol_chs_overlap)
length(unique(sol$name))/length(sol$name)


pooled_other <- pooled %>% 
  tidylog::left_join(sol_chs_overlap)
length(unique(sol$name))/length(sol$name)


# Make List 
other_met <- list(solar = sol_other, 
                  chs   = chs_other,
                  pooled = pooled_other)


# Save data 
write_rds(other_met,
          fs::path(dir_results_mixtures, 
                   "ind_met_effect_ests", 
                   "other_met_pathway effect estimates.rds"))
# Make Plots ----------------------------------------------------------

## solar------------------------------
(sol_met_plot <- ggplot(sol_other,
                        aes(x = met_name, 
                            y = estimate_beta_sol_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_sol_mixture , 
                     ymax = ucl_beta_sol_mixture, 
                     linetype = sig_with_direction), 
                 width = 0) +
   geom_point(aes(shape = sig_with_direction, 
                  fill = sig_with_direction, 
                  size = sig_with_direction)) + 
   # geom_text(aes(label = qsig_sol, y = 2.1),
   #           size = 3.5, hjust = 0, vjust = .5) + 
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   # Scales
   scale_fill_manual(values = c( "grey50","black",  "white"), name = NULL) +
   scale_shape_manual(values = c(21, 21, 22), name = NULL) +
   scale_linetype_manual(values = c(1,1,1), name = NULL) +
   scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
   scale_y_continuous(limits = c(-1.75, 3)) +
   coord_flip() + 
   facet_grid(pathway_reduced ~ .,  scales = "free_y", space = "free_y") +
   ylab("PFAS Mixture Effect ψ (95% BCI)") +
   # ylab("PFAS Exposure Effect Estimate ψ (95% BCI)") +
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

# Save 
ggsave(sol_met_plot, 
       filename = fs::path(
         dir_reports, 
         "Figure 5 R1 associations of PFAS and other metabolites SOL only.jpg"), 
       width = 8, height = 2.5, bg = "white")

