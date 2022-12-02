# Volcano Plots
library("jag2")
library("tidyverse")
source(here::here("0_project_setup", "!directories.r"))

# SOLAR Volcano Plots ----------------------------------------------
mwas_results <- read_csv(
  file = fs::path(dir_results_mixtures, 
                  "all_pfas_mixtures_results_hyper_g.csv")) %>% 
  as_tibble()

# Get reduced dataset
mwas_reduced <- mwas_results %>% 
  filter(p_value  < 0.99)

# Volcano Plot
solar_volcano_plot <- ggplot(mwas_reduced %>% filter(cohort == "SOL"),
                              aes(x = estimate_beta, y = -log10(p_value))) +
    geom_point(size = 1, alpha = 0.5) + 
    geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
    facet_wrap(~exposure, scales = "free") +
    ylab("-log10 p") +
    xlab("Effect Estimate") +
    ggtitle("SOLAR")

# Save and clean env.
ggsave(solar_volcano_plot,
       filename = fs::path(
         dir_reports,
         "Volcano plots",
         "SOLAR Mixtures analysis volcano plots p_original_hyper_g.jpg"), 
       width = 6, height = 5)
rm(solar_volcano_plot)


# CHS Volcano Plots ----------------------------------------------
chs_volcano_plot <- ggplot(mwas_reduced %>% filter(cohort == "CHS"),
                           aes(x = estimate_beta, y = -log10(p_value))) +
  geom_point(size = 1, alpha = 0.5) + 
  geom_vline(xintercept = 0, color = "grey20", linetype = 2) +
  facet_wrap(~exposure, scales = "free") +
  ylab("-log10 p") +
  xlab("Effect Estimate") +
  ggtitle("CHS")

# # Save results
ggsave(chs_volcano_plot,
       filename =  fs::path(
         dir_reports,
         "Volcano plots",
         "chs Mixtures analysis volcano plots_p_from_original_hyper_g.jpg"),
       width = 6.5, height = 5)

