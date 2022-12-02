# Make coefficient plots based on dougs comments
rm(list = ls())
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))


# Read data
annotated_sig_ecs_ee_tyrosine <- read_rds(
  here::here(
    dir_results_mixtures,
    "Sig annotated metabolite effect estimates sol chs_v4.RDS") ) %>% # New
  select(-contains("blank")) %>% 
  tidylog::filter(str_detect(pathway, "Tyrosine"))


length(unique(annotated_sig_ecs_ee_tyrosine$met_name))
# table(annotated_sig_ecs_ee_tyrosine$q_value_sol_mixture<0.05)
# table(annotated_sig_ecs_ee_tyrosine$q_value_sol_mixture<0.05)

# Change metabolite names to those associated with Tyrosine metabolism --------------
annotated_sig_ecs_ee <- annotated_sig_ecs_ee_tyrosine %>% 
  tidylog::mutate(met_name_tyr_pw = case_when(
    met_name == "Ascorbate; Glucuronolactone; D-glucurono-3,6-lactone" ~ 
      "Ascorbate", 
    met_name == "Pyridoxine; Norepinephrine" ~ "Norepinephrine",
    met_name == "Dopaquinone; Leucodopachrome" ~ "Dopaquinone",
    met_name == "Formylanthranilic acid; Noradrenochrome" ~ 
      "Noradrenochrome",
    met_name == "Vanylglycol; Phosphorylcholine" ~ "Vanylglycol",
    met_name == "Adrenochrome; Hippuric acid" ~ "Hippuric acid", 
    met_name == "Homovanillic acid; 3-Methoxy-4-hydroxyphenylglycolaldehyde" ~
      "Homovanillic acid", 
    met_name == "3-Methoxytyramine; Epinine" ~ "3-Methoxytyramine", 
    met_name == "5-Hydroxytryptophol; 1,2-dehydrosalsolinol" ~ 
      "1,2-dehydrosalsolinol", 
    met_name == "Phenylacetylglutamine; Acetyl-N-formyl-5-methoxykynurenamine" ~ 
      "Phenylacetylglutamine", 
    met_name == "L-Glutamic acid; D-Glutamic acid; L-4-Hydroxyglutamate semialdehyde" ~
      "L-Glutamic acid",
    met_name == "L-Glutamic acid; D-Glutamic acid; DL-Glutamate; L-4-Hydroxyglutamate semialdehyde" ~ 
      "L-Glutamic acid",
    met_name == "sulfuric acid [4-(2-aminoethyl)phenyl] ester" ~ 
      "Tyramine-O-sulfate", 
    met_name == "4-Hydroxyphenylacetaldehyde; Phenylacetic acid" ~ 
      "4-Hydroxyphenylacetaldehyde", 
    met_name == "3-methyl pyruvic acid; Acetoacetic acid; Succinic acid semialdehyde; 2-Methyl-3-oxopropanoic acid; (S)-Methylmalonic acid semialdehyde" ~ "Acetoacetic acid",
    met_name == "Pyruvic acid; Malonic semialdehyde" ~ "Pyruvic acid",
    TRUE ~ met_name))

# Add tyrosine sub pathways
annotated_sig_ecs_ee <- annotated_sig_ecs_ee %>% 
  tidylog::mutate(
    tyr_subpath = case_when(
      met_name_tyr_pw %in% c("3,4-Dihydroxyphenylglycol","Metanephrine",
                             "Noradrenochrome", "Vanylglycol","Ascorbate",
                             "Adrenochrome","Norepinephrine",
                             "Norepinephrine sulfate", "L-Dopa", 
                             "1,2-dehydrosalsolinol","3-Methoxytyramine",
                             "3-O-methyldopa","Homovanillic acid",
                             "Homovanillin") ~
        "Catecholamine biosynthesis\nand degredation", 
      
      met_name_tyr_pw %in%  c("Phenylacetaldehyde",
                              "Phenylacetylglutamine", 
                              "Hippuric acid") ~ 
        "Phenylalanine metabolism",
      
      met_name_tyr_pw %in%  c("4-Hydroxyphenylacetaldehyde", 
                              "L-Glutamic acid", 
                              "Tyramine-O-sulfate", 
                              "Acetoacetic acid",
                              "Pyruvic acid") ~
        "Tyrosine metabolism and\ndegredation", 
      
      met_name_tyr_pw == "Dopaquinone" ~ "Melanin biosynthesis", 
      met_name_tyr_pw == "Thyroxine" ~ "Thyroid Hormone Biosynthesis"))


# Modify table based on Dougs feedback (comments are directly from Doug)
sig_ecs_reduced <- annotated_sig_ecs_ee %>% 
  tidylog::filter(
    #For query_mass == 167.0701851, this is homovanillin, not vanyglycl-
    # multiple adducts at the same retention time corresponding to this compound (~125 sec)
    !(met_name == "Vanylglycol; Phosphorylcholine" & query_mass == 167.0701851), # Removes 1 row
    
    # This is not homovanilian, since it does not match the other 3 adduct retention times:
    !(met_name == "Homovanillin" & retention_time  == "72.3"), # Removes 1 row
    
    #This is a very weird adduct, I would say this is not correct:
    matched_form != "M-HCOOH+H[1+]", # Removes 1 row
    
    # Duplicate mz/rt: 
    !(matched_form == "M-CO2+H[1+]" & query_mass == 152.070757656544),
    
    # Duplicate mz/rt: 
    !(matched_form == "M-H2O+H[1+]" & query_mass == 152.070757656544),
    
    # Duplicate mz/rt: 
    !(matched_form == "M-H2O+H[1+]" & query_mass == 180.0655298),

    !(matched_form == "M-H[-]" & name == "146.045966215623_339.881479830506")
    
  ) 


# Get Mass difference for the tyrosine metabolites 
sig_ecs_reduced2 <- sig_ecs_reduced %>% 
  group_by(met_name, name) %>% 
  tidylog::mutate(
    mass_diff_ls = str_split(mass_diff, "; "), 
    mass_diff_tyr_met = case_when(
      met_name == "Formylanthranilic acid; Noradrenochrome" ~ mass_diff_ls[[1]][2],
      met_name == "Vanylglycol; Phosphorylcholine" ~ mass_diff_ls[[1]][2],
      met_name == "Ascorbate; Glucuronolactone; D-glucurono-3,6-lactone" ~ 
        mass_diff_ls[[1]][1], 
      TRUE ~ mass_diff))  %>% 
  ungroup() %>% 
  select(met_name_tyr_pw,  tyr_subpath, mass_diff_tyr_met,  everything(), 
         -mass_diff_ls,  -casnum, -HRE_standard)   


# Reorder tyrosine sub path by number of cpds within path
sig_ecs_reduced3 <- sig_ecs_reduced2 %>% 
  mutate(tyr_subpath = fct_infreq(tyr_subpath), 
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
sol <- sig_ecs_reduced3 %>% 
  select(met_name_tyr_pw:mass_diff,empirical_compound, 
         contains("sol")) %>% 
  tidylog::filter(ucl_beta_sol_mixture <= 0 | 
                    lcl_beta_sol_mixture >= 0) %>% 
  group_by(met_name_tyr_pw) %>% 
  arrange(-abs(estimate_beta_sol_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name_tyr_pw = fct_reorder(met_name_tyr_pw, 
                                       estimate_beta_sol_mixture))


## Subset CHS (filter sig features) ----------------
chs <- sig_ecs_reduced3 %>% 
  select(met_name_tyr_pw:mass_diff,empirical_compound, 
         contains("chs")) %>% 
  tidylog::filter(ucl_beta_chs_mixture <= 0 | 
                    lcl_beta_chs_mixture >= 0) %>% 
  group_by(met_name_tyr_pw) %>% 
  arrange(-abs(estimate_beta_chs_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name_tyr_pw = fct_reorder(met_name_tyr_pw, 
                                       estimate_beta_chs_mixture))


## Subset Pooled (filter sig features) ----------------
pooled <- sig_ecs_reduced3 %>% 
  select(met_name_tyr_pw:mass_diff,empirical_compound, 
         contains("pooled")) %>% 
  tidylog::filter(ucl_beta_pooled_mixture <= 0 | 
                    lcl_beta_pooled_mixture >= 0) %>% 
  group_by(met_name_tyr_pw) %>% 
  arrange(-abs(estimate_beta_pooled_mixture)) %>% 
  tidylog::slice_head() %>% 
  ungroup() %>% 
  mutate(met_name_tyr_pw = fct_reorder(met_name_tyr_pw, 
                                       estimate_beta_pooled_mixture))


## Get data on overlap ---------------------
sol_chs_overlap <- tidylog::full_join(
  sol %>% 
    select(met_name_tyr_pw, 
           estimate_beta_sol_mixture), 
  chs %>% 
    select(met_name_tyr_pw, 
           estimate_beta_chs_mixture)) %>% 
  mutate(sig_dir = estimate_beta_sol_mixture*estimate_beta_chs_mixture, 
         sig_both = !is.na(sig_dir), 
         sig_with_direction = case_when(
           sig_both & sig_dir < 0 ~  
             "Opposite direction of association in SOLAR vs. CHS",
           sig_both & sig_dir > 0 ~ 
             "Same direction of association in SOLAR and CHS", 
           TRUE ~  "Only significant in one cohort") %>% 
           fct_relevel(
             "Same direction of association in SOLAR and CHS", 
             "Only significant in one cohort",
             "Opposite direction of association in SOLAR vs. CHS")) %>% 
  select(-estimate_beta_chs_mixture, -estimate_beta_sol_mixture)

# Join with solar/chs data
sol_tyr <- sol %>% 
  tidylog::left_join(sol_chs_overlap) 
# Check number of unique mz/rt 
length(unique(sol_tyr$name))/length(sol_tyr$name)

chs_tyr <- chs %>% 
  tidylog::left_join(sol_chs_overlap) %>% 
  mutate(qsig = if_else(q_value_chs_mixture <0.05, "*", ""),
         met_name_sig = str_c(met_name_tyr_pw, qsig))
# Check number of unique mz/rt
length(unique(chs_tyr$name))/length(chs_tyr$name)


pooled_tyr <- pooled %>% 
  tidylog::left_join(sol_chs_overlap)  %>% 
  mutate(qsig = if_else(q_value_pooled_mixture <0.05, "*", ""),
         met_name_sig = str_c(met_name_tyr_pw, qsig))
# Check number of unique mz/rt
length(unique(pooled_tyr$name))/length(pooled_tyr$name)


# Make List 
tyrosine_met <- list(solar = sol_tyr, chs = chs_tyr, pooled = pooled_tyr)


# Save data 
write_rds(tyrosine_met,
          fs::path(dir_results_mixtures,
                   "ind_met_effect_ests",
                   "aeromatic amino acid pooled effect estimates.rds"))

# Make Plots ----------------------------------------------------------
## Plot solar------------------------------
max(tyrosine_met$solar$ucl_beta_sol_mixture)
(sol_met_plot <- ggplot(tyrosine_met$solar, 
                        aes(x = met_name_tyr_pw, 
                            y = estimate_beta_sol_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_sol_mixture , 
                     ymax = ucl_beta_sol_mixture), 
                 width = 0, size = .75) + 
   geom_point(aes(shape = sig_with_direction, 
                  fill = sig_with_direction, 
                  size = sig_with_direction)) + 
   # geom_text(aes(label = qsig_sol, y = 2.1),
   #           size = 3.5, hjust = 0, vjust = .5) + 
    
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   facet_grid(tyr_subpath ~ .,  scales = "free_y", space = "free_y") +
   
   # Scales
   scale_fill_manual(values = c( "black", "grey50", "white"), name = NULL) +
   scale_shape_manual(values = c(21, 21, 22), name = NULL) +
   scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
   scale_y_continuous(limits = c(-1.75, 3.2)) +
   coord_flip(clip = "off") +    
   # ylab("PFAS Exposure Effect Estimate ψ (95% BCI)") +
   theme(axis.title.y = element_blank(), 
         axis.title.x = element_blank(), 
         strip.text.x = element_blank(), 
         panel.background = element_rect(fill="grey95"), 
         strip.background = element_rect(fill = "white"),
         
         legend.position = "none",
         legend.direction = "vertical",
         strip.text.y = element_text(angle = 0, hjust = 0),
         text = element_text(size = 15),
         panel.grid = element_blank(),
         axis.line.x = element_line(color = "black"),
         axis.line.y = element_line(color = "black"))) 

## Plot CHS ------------------------------
(chs_metplot <- ggplot(tyrosine_met$chs,
                       aes(x = met_name_tyr_pw, 
                           y = estimate_beta_chs_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_chs_mixture , 
                     ymax = ucl_beta_chs_mixture, 
                     linetype = sig_with_direction), 
                 width = 0, size = .75) + 
   geom_point(aes(shape = sig_with_direction, 
                  fill = sig_with_direction, 
                  size = sig_with_direction)) + 
   # geom_text(aes(label = qsig_chs, y = 2.1),
   #                      size = 3.5, hjust = 0, vjust = .5) + 
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   facet_grid(tyr_subpath~"1",  scales = "free_y", space = "free_y") +
   # Scales
   scale_fill_manual(values = c( "black", "grey50", "white"), name = NULL) +
   scale_shape_manual(values = c(21, 21, 22), name = NULL) +
   scale_linetype_manual(values = c(1,1,1), name = NULL) +
   scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
   scale_y_continuous(limits = c(-1.75, 3.2)) +
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
         strip.text.y = element_text(angle = 0, hjust = 0),
         text = element_text(size = 15),
         panel.grid = element_blank(),
         axis.line.x = element_line(color = "black"),
         axis.line.y = element_line(color = "black")))


# Combine Plots
(figure_2 <- plot_grid(NULL, sol_met_plot, 
                       NULL, chs_metplot, 
                       ncol = 1, align = "v", 
                       rel_heights = c(.05, 1,.05, .85),
                       labels = c("A) SOLAR","",
                                  "B) CHS", "")))

#Save 
ggsave(figure_2,
       filename = fs::path(dir_reports,
                           "Figure 2 R1 associations of PFAS and tyr metabolites.jpg"),
       width = 7.5, height = 10, bg= "white")


# rm(sol, chs, sol_chs_overlap, sol_tyr, chs_tyr)