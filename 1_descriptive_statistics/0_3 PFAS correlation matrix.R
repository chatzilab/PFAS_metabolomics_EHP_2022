# Correlation plots
library(correlation)
library(GGally)

source(here::here("0_project_setup", "!functions.R"))
source(here::here("0_project_setup", "!set_exposure_outcome_vars.R"))

# Function for calculating correlations for plot:
cor_func <- function(data, mapping, method, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor(x, y, method=method, use='complete.obs')
  
  ggally_text(
    label = formatC(corr, format = "f", 2), 
    mapping = aes(),
    xP = 0.5, yP = 0.5,
    color = 'black',
    ...
  )
}

# Get analysis data set, select important vars ------------------

# SOLAR
solar_exposure <- exposure_outcome$solar %>% 
  select(wave, all_of(exposures_continuous))%>% 
  rename_all(toupper) %>% 
  rename(PFHxS = PFHXS, 
         PFHpS = PFHPS)

# CHS
chs_exposure <- exposure_outcome$chs %>% 
  select(all_of(exposures_continuous))%>% 
  rename_all(toupper) %>% 
  rename(PFHxS = PFHXS, 
         PFHpS = PFHPS)

# Make Final Plots -----------------------------------------------------
# Solar 
sol_corplot <- ggpairs(solar_exposure, columns = 2:7,
                       upper = list(continuous = wrap(cor_func,
                                                      method = "spearman", 
                                                      size = 5)), 
                       lower = list(continuous = wrap("points",alpha = .5, size=.75))) + 
  xlab("Plasma Concentration (\u00b5g/L)")  +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12), 
        axis.text.y = element_text(size=10))

# CHS
chs_corplot <- ggpairs(chs_exposure, 
                       upper = list(continuous = wrap(cor_func,
                                                      method = "spearman", 
                                                      size = 5)), 
                       lower = list(continuous = wrap("points",alpha = .5, size=.75))) + 
  xlab("Plasma Concentration (\u00b5g/L)") +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, size=12),
        axis.text.y = element_text(size=10))



# Combine SOLAR and CHS figs ------------------
sol_g <- grid::grid.grabExpr(print(sol_corplot))
chs_g <- grid::grid.grabExpr(print(chs_corplot))

# Combine figures
fig_s1 <- plot_grid(NULL, NULL, 
                   sol_g, 
                   chs_g,
                   labels = c("A. SOLAR Cohort", 
                              "B. CHS Cohort", 
                              "", ""),
                   rel_widths = c(1,1),
                   rel_heights = c(.02,1),
                   hjust = 0, vjust = 1,
                   label_size = 14,
                   # align = "vh",
                   axis = "lr",
                   nrow = 2)


## Save Correlation Plot All PFAS -------------------------
ggsave(plot = fig_s1,
       filename = path(dir_reports, 
                       "Figure S1 SOLAR CHS Correlation Plot All PFAS.jpeg"), 
       width = 16, 
       height = 7)
  
