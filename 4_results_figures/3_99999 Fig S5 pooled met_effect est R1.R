# Make coefficient plots of pooled data
rm(list = ls())
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))


# Read data
tyrosine_met <- read_rds(
  fs::path(
    dir_results_mixtures, 
    "ind_met_effect_ests", 
    "aeromatic amino acid pooled effect estimates.rds"))$pooled %>% 
  select(-met_name) %>% 
  rename(met_name = met_name_tyr_pw)


nonaro_met <- read_rds(
  fs::path(dir_results_mixtures, 
           "ind_met_effect_ests", 
           "nonaeromatic aa effect estimates.rds"))$pooled

fa_met <- read_rds(
  fs::path(dir_results_mixtures, 
           "ind_met_effect_ests", 
           "lipid effect estimates.rds"))$pooled

other_met <- read_rds(
  fs::path(dir_results_mixtures, 
           "ind_met_effect_ests", 
           "other_met_pathway effect estimates.rds"))$pooled

# Bind rows 
pooled_met_ee <- bind_rows(tyrosine_met, nonaro_met, fa_met, other_met)


# Percentage unique mzrt
length(unique(pooled_met_ee$name))/length((pooled_met_ee$name))

dup_mzrt(pooled_met_ee)


pooled_met_ee_fin <- pooled_met_ee %>% 
  mutate(sig_with_direction = as.character(sig_with_direction), 
         sig_with_direction = replace_na(
           sig_with_direction, 
           "Not significant in either cohort individually") %>% 
           fct_relevel("Same direction of association in SOLAR and CHS",
                       "Only significant in one cohort", 
                       "Not significant in either cohort individually"), 
         pathway = str_remove(pathway,
                              "; Urea cycle/amino group metabolism") %>% 
           str_replace("Prostaglandin formation from arachidonate", 
                       "Prostaglandin formation\nfrom arachidonate")) %>% 
  select(qsig_pooled, ucl_beta_pooled_mixture, everything(), -contains("tyr_"))
# qsig = if_else(q_value_pooled_mixture <0.05, "*", ""),
# met_name_sig = str_c(met_name, qsig))

colnames(pooled_met_ee_fin)
# Make Plots ----------------------------------------------------------
(pooled_metplot <- ggplot(pooled_met_ee_fin,
                          aes(x = met_name, 
                              y = estimate_beta_pooled_mixture)) + 
   geom_errorbar(aes(ymin = lcl_beta_pooled_mixture , 
                     ymax = ucl_beta_pooled_mixture, 
                     linetype = sig_with_direction), 
                 width = 0, size = .75) +
   geom_point(aes(shape = sig_with_direction, 
                  fill = sig_with_direction, 
                  size = sig_with_direction)) +  
   geom_text(aes(label = qsig_pooled, 
                 y = ucl_beta_pooled_mixture+.1),
             size = 3.5, hjust = 0, vjust = .5) + 
   geom_hline(yintercept = 0, linetype = 3, color = "grey50") +
   facet_grid(pathway~"1",  scales = "free_y", space = "free_y") +
   # Scales
   scale_fill_manual(values = c( "black", "grey50", "white"), name = NULL) +
   scale_shape_manual(values = c(21, 21, 22), name = NULL) +
   scale_linetype_manual(values = c(1,1,1), name = NULL) +
   scale_size_manual(values = c(2.5, 2.5, 2.5), name = NULL) +
   scale_y_continuous(limits = c(-1.75, 3.5)) +
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



#Save 
ggsave(pooled_metplot,
       filename = fs::path(
         dir_reports,
         "Figure S5 R1 pooled metabolite effect ests.jpg"),
       width = 8.5, height = 5)


# rm(sol, chs, sol_chs_overlap, sol_tyr, chs_tyr)
