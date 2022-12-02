# Upset Plots
library(ComplexUpset)
rm(list = ls())
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!functions.R"))

# read data from 2_5 
# in data v3, the change is that we selected metabolites associated with significant pathway changes in either cohort
upplot_data <- read_rds(
  here::here(
    dir_results_mixtures,
    # "Ecd effect estimates and pathways sol chs for upset.RDS")) # Original
    "Ecd effect estimates and pathways sol chs for upset_v3.RDS"))

# Update superpathway
upplot_data <- upplot_data %>%
  mutate(super_pathway_new = case_when(str_detect(pathway, "Alkaloid") ~ "Other", 
                                       str_detect(pathway, "Nitrogen") ~ "Other", 
                                       str_detect(pathway, "cytochrome P450") ~ "Other",
                                       str_detect(super_pathway, "Aromatic Amino Acids") ~ 
                                         "Aromatic amino acid\nmetabolism",
                                       str_detect(super_pathway, "Amino acid") ~ 
                                         "Non-aromatic amino\nacid metabolism",
                                       str_detect(super_pathway, 
                                                  "cofactors") ~ 
                                         "Metabolism of cofactors\nand vitamins", 
                                       TRUE ~ super_pathway) %>% 
           fct_relevel("Aromatic amino acid\nmetabolism", 
                       "Non-aromatic amino\nacid metabolism"))



# Data contains all 14 pathways:
length(unique(upplot_data$pathway))

# SOLAR --------------------------------------------
ecds_for_upset <- upplot_data %>% 
  group_by(met_name) %>% 
  mutate(empirical_compound_all = str_c(empirical_compound, collapse = ";"), 
         empirical_compound = empirical_compound[1]) %>% 
  ungroup() %>% 
  tidylog::filter(
    # From Doug: This is a weird adduct that probably is not correct
    matched_form != "M-HCOOH+H[1+]", 
    # This is not homovanilian, since it does not match the other 3 adduct retention times:
    !(met_name == "Homovanillin" & name == "185.080994057582_72.2749441021845")) 


# SOLAR -------------------------------------------------------------------
# Calculate max estimate and sig
ecds_for_upset_sol <- ecds_for_upset %>% 
  tidylog::filter(ucl_beta_sol_mixture <= 0 |
                    lcl_beta_sol_mixture >= 0) %>% 
  group_by(super_pathway_new, empirical_compound) %>%
  tidylog::summarise(max_beta = max(estimate_beta_sol_mixture), 
                     min_beta = min(estimate_beta_sol_mixture), 
                     max_est = if_else(abs(max_beta) > abs(min_beta), 
                                       max_beta, 
                                       min_beta), 
                     ec = str_c(empirical_compound, collapse = ";"), 
                     met_name = str_c(met_name, collapse = ";"), 
                     mzrt = unique(name) %>% str_c(collapse = "; ")) %>% 
  ungroup() %>% 
  group_by(met_name) %>% 
  arrange(-abs(max_est)) %>% 
  tidylog::slice_head() %>% 
  ungroup()


length(unique(ecds_for_upset_sol$met_name))/
  length(ecds_for_upset_sol$met_name)

length(unique(ecds_for_upset_sol$mzrt))/
  length(upplot_data_sol$mzrt)

# Pivot wider to get table of 0/1s
upplot_data_sol <- ecds_for_upset_sol %>% 
  select(-contains("beta")) %>%
  mutate(inpath = 1) %>% 
  tidylog::pivot_wider(values_from = inpath, 
              names_from = super_pathway_new,
              id_cols = c(empirical_compound, 
                          met_name,
                          ec,
                          mzrt,
                          "max_est")) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)), 
         `Direction of Association` = if_else(max_est > 0, 
                                              "Positive Association",
                                              "Negative Association")) %>% 
  as.data.frame() 

# Number of metabolites from the 14 total pathways
length(unique(upplot_data_sol$empirical_compound))/length(upplot_data_sol$empirical_compound)
length(unique(upplot_data_sol$mzrt))/length(upplot_data_sol$mzrt)


## Set intersection order for both graphs ------------------------------
intersect_order <- c( "Aromatic amino acid\nmetabolism",
                      "Lipid metabolism", 
                    "Non-aromatic amino\nacid metabolism",
                     "Metabolism of cofactors\nand vitamins", 
                     "Other") %>% rev() # sol only

table(ecds_for_upset_sol$super_pathway_new)
## UpSet graph showing the number of empirical cpds from each pathway --------------
(upset_sol <- upset(
  upplot_data_sol, 
  min_size = 0, 
  stripes="white",
  sort_sets=FALSE, 
  intersect = intersect_order,
  # Adjusting the main bargraph (set intersection sizes)
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        counts=FALSE,
        mapping=aes(fill=`Direction of Association`)) 
      + theme_cowplot() 
      + ylim(c(0,20))
      + theme(axis.text.x = element_blank(), 
              legend.position = "none",
              axis.title.x = element_blank(), 
              axis.title.y = element_text(vjust = -35)
      ))),
  # Adjusting the lower left bar graph (set sizes)
  set_sizes=(
    upset_set_size(
      geom=geom_bar(
        aes(fill=`Direction of Association`, x=group),
        width=0.8), 
      position='right') 
    + theme_cowplot() 
    + theme(axis.text.y = element_blank(), 
            axis.title.y = element_blank(), 
            legend.position = "none", 
            axis.line.y = element_blank(), 
            axis.ticks.y = element_blank())), 
  # Adjusting the intersection matrix
  matrix=(
    intersection_matrix(
      geom=geom_point(shape='circle',size=3),
      outline_color=list(active='black',inactive='white'))) ,
  themes = list(intersections_matrix=
                  theme_cowplot() + 
                  theme(axis.text.x = element_blank(), 
                        axis.ticks = element_blank(), 
                        axis.title.y = element_blank()))))


# CHS --------------------------------------------
ecds_for_upset_chs <- ecds_for_upset %>% 
  tidylog::filter(ucl_beta_chs_mixture <= 0 |
                    lcl_beta_chs_mixture >= 0) %>%
  group_by(super_pathway_new, empirical_compound) %>% 
  tidylog::summarise(max_beta = max(estimate_beta_chs_mixture), 
                     min_beta = min(estimate_beta_chs_mixture), 
                     max_est = if_else(abs(max_beta) > abs(min_beta), 
                                       max_beta, 
                                       min_beta)) %>% 
  ungroup()

# Number of metabolites from the 7 enriched pathways
length(unique(ecds_for_upset_chs$empirical_compound))


# Pivot wider to get table of 0/1s
upplot_data_chs <- ecds_for_upset_chs %>% 
  select(-contains("beta")) %>%
  mutate(inpath = 1) %>% 
  pivot_wider(values_from = inpath, 
              names_from = super_pathway_new,
              id_cols = c(empirical_compound, "max_est")) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)), 
         `Direction of Association` = if_else(max_est > 0, 
                                              "Positive Association",
                                              "Negative Association")) %>% 
  as.data.frame()


## UpSet graph showing the number of empirical cpds from each pathway --------------
(upset_chs <- upset(
  upplot_data_chs, 
  min_size = 1, 
  stripes="white",
  sort_sets=FALSE, 
  intersect = intersect_order,
  # Adjusting the main bargraph (set intersection sizes)
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        counts=FALSE,
        mapping=aes(fill=`Direction of Association`)) 
      + theme_cowplot() 
      + ylim(c(0,15))
      + theme(axis.text.x = element_blank(), 
              legend.position = "none",
              axis.title.x = element_blank(), 
              axis.title.y = element_text(vjust = -35)
      ))),
  # Adjusting the lower left bargraph (set sizes)
  set_sizes=(
    upset_set_size(
      geom=geom_bar(
        aes(fill=`Direction of Association`, x=group),
        width=0.8), 
      position='right') 
    + theme_cowplot() 
    + theme(axis.text.y = element_blank(), 
            axis.title.y = element_blank(), 
            legend.position = "none", 
            axis.line.y = element_blank(), 
            axis.ticks.y = element_blank())), 
  # Adjusting the intersection matrix
  matrix=(
    intersection_matrix(
      geom=geom_point(shape='circle',size=3),
      outline_color=list(active='black',inactive='white'))) ,
  themes = list(intersections_matrix=
                  theme_cowplot() + 
                  theme(axis.text.x = element_blank(), 
                        # axis.text.y = element_blank(),
                        axis.ticks = element_blank(), 
                        axis.title.y = element_blank()))))



# Combine figures ----------------------------------------
# Get legend ----------------------------------
(upset_legend <- ggplot(upplot_data_chs, 
                        aes(x = `Direction of Association`, 
                            y = empirical_compound, 
                            color = `Direction of Association`)) + 
   geom_point(shape = "square") +
   theme(legend.position = "bottom", 
         legend.justification="center", 
         legend.box.background = 
           element_rect(color = "black"),
         legend.text = element_text(size = 14),
         legend.box.margin = margin(4,4,4,4)))



#Panel a
upsetfig_a <- plot_grid(upset_sol,
                        upset_chs, 
                        rel_widths = c(1., 1), 
                        labels = c("A. SOLAR", 
                                   "B. CHS"))

(upsetfig_final <- plot_grid(upsetfig_a,
                             cowplot::get_legend(upset_legend), 
                             nrow = 2,
                             rel_heights = c(2.5, .2)))


ggsave(upsetfig_final,
       filename = fs::path(dir_reports, "Figure S2 Upset Plot SOL CHS R1.jpg"), 
       height = 6.5, width = 13)




# Upset plot pooled analysis ------------------
ecds_for_upset_chs <- ecds_for_upset %>% 
  tidylog::filter(ucl_beta_pooled_mixture <= 0 |
                    lcl_beta_pooled_mixture >= 0) %>%
  group_by(pathway, empirical_compound) %>% 
  tidylog::summarise(max_beta = max(estimate_beta_pooled_mixture), 
                     min_beta = min(estimate_beta_pooled_mixture), 
                     max_est = if_else(abs(max_beta) > abs(min_beta), 
                                       max_beta, 
                                       min_beta)) %>% 
  ungroup()

# Number of metabolites from the 7 enriched pathways
length(unique(ecds_for_upset_chs$empirical_compound))


# Pivot wider to get table of 0/1s
upplot_data_chs <- ecds_for_upset_chs %>% 
  select(-contains("beta")) %>%
  mutate(inpath = 1, 
         pathway = str_replace_all(pathway, "Metabolism", "metabolism")) %>% 
  pivot_wider(values_from = inpath, 
              names_from = pathway,
              id_cols = c(empirical_compound, "max_est")) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0)), 
         `Direction of Association` = if_else(max_est > 0, 
                                              "Positive Association",
                                              "Negative Association")) %>% 
  as.data.frame()


## UpSet graph showing the number of empirical cpds from each pathway --------------
(upset_chs <- upset(
  upplot_data_chs, 
  min_size = 1, 
  stripes="white",
  sort_sets=FALSE, 
  intersect = intersect_order,
  # Adjusting the main bargraph (set intersection sizes)
  base_annotations=list(
    'Intersection size'=(
      intersection_size(
        counts=FALSE,
        mapping=aes(fill=`Direction of Association`)) 
      + theme_cowplot() 
      + ylim(c(0,15))
      + theme(axis.text.x = element_blank(), 
              legend.position = "none",
              axis.title.x = element_blank(), 
              axis.title.y = element_text(vjust = -50)
      ))),
  # Adjusting the lower left bargraph (set sizes)
  set_sizes=(
    upset_set_size(
      geom=geom_bar(
        aes(fill=`Direction of Association`, x=group),
        width=0.8), 
      position='right') 
    + theme_cowplot() 
    + theme(axis.text.y = element_blank(), 
            axis.title.y = element_blank(), 
            legend.position = "none", 
            axis.line.y = element_blank(), 
            axis.ticks.y = element_blank())), 
  # Adjusting the intersection matrix
  matrix=(
    intersection_matrix(
      geom=geom_point(shape='circle',size=3),
      outline_color=list(active='black',inactive='white'))) ,
  themes = list(intersections_matrix=
                  theme_cowplot() + 
                  theme(axis.text.x = element_blank(), 
                        # axis.text.y = element_blank(),
                        axis.ticks = element_blank(), 
                        axis.title.y = element_blank()))))
