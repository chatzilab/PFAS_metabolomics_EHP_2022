# Set up Pooled data_for_mixtures_analysis on HPC

source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!set_exposure_outcome_vars.R"))
source(here::here("0_project_setup", "!load_data.R"))

# Get BHRMA Function
source(here::here("2_mixtures_MWAS",
                  "0_0_BHRMA.g_function.R"))
# Get mixtures components
source(here::here("0_project_setup", "!set_exposure_outcome_vars.R"))
mixtures_name = "all_pfas"
mixtures_components = exposures_continuous

# Modify data -------------------------------------------
# Change ftdata from list to dataframe
sol_metab_dat <- purrr::reduce(ftdata$solar, .f = left_join)
chs_metab_dat <- purrr::reduce(ftdata$chs, .f = left_join)
metab_dat <- bind_rows(sol_metab_dat, chs_metab_dat)

# Covariates to be included in models: 
# sex
# tanner (tanner)
# ses.num (edu_house)
# cohort_wave_num
# Note: age and BMI are excluded as covariates because they are extremly highly correlated with 
# tanner stage and cohort/wave

## Select Covariates and exposures from SOLAR and CHS --------------------
# SOLAR
solar_eo <- exposure_outcome$solar %>% 
  select(id, age, sex, edu_house, wave, tanner, 
         all_of(mixtures_components)) 

# CHS
chs_eo   <- exposure_outcome$chs   %>% 
  select(id, age, sex, edu_house, all_of(mixtures_components)) %>% 
  mutate(id = as.character(id), 
         wave = "CHS", 
         edu_house = str_replace(edu_house, 
                                 "Standard college or university graduation",
                                 "Completed college/university" ), 
         tanner = as.factor(5))

# Bind Rows
exp_out <- bind_rows(solar_eo, chs_eo)

# factors to numeric
library(fastDummies)
exp_out_num <- exp_out %>% 
  mutate(
    edu_house = fct_collapse(
      edu_house, 
      "College degree" = c("Graduate professional training (graduate degree)",
                           "Completed college/university")) %>% 
      str_remove(" \\(at least one year\\) or specialized training")) %>% 
  dummy_cols(select_columns = c("sex", "edu_house", "wave", "tanner"), 
             remove_most_frequent_dummy = TRUE, 
             remove_selected_columns = TRUE) %>% 
  janitor::clean_names()

# Get new covariate colnames
cov_num_colnames <- colnames(exp_out_num %>% 
                               select(#age, #bmi, 
                                 sex_female:tanner_4))


# Join exposures and metabolites 
full_data <- tidylog::left_join(exp_out_num, 
                                metab_dat, 
                                by = "id")


# Set up analysis datasets ------------------
# X: PFAS
X.obs = full_data[mixtures_components] %>%
  mutate(across(everything(),
                ~scale(log2(.))))

# Y: Metabolite
Y = full_data %>%
  dplyr::select(colnames(metab_dat)[2]:ncol(full_data)) %>% 
  scale(center = F, scale = F) 

# U: Covariates
# sex
# tanner (tanner)
# ses.num (edu_house)
# cohort_wave_num
U = full_data %>% 
  select(contains("sex"), contains("tanner"), contains("edu"), contains("wave"))

P = ncol(X.obs)
LOD = c(0.01,0.05, 0.01, 0.01, 0.01, 0.43)
profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P)
# exposure.Names = colnames(X.obs)


rm(list = setdiff(ls(), c("X.obs",
                          "Y",
                          "U",
                          "LOD",
                          "profiles", 
                          "ridge.BDL.model", 
                          "mixtures_name",
                          "BHRMA.g")))

save.image(file = fs::path(dirname(here::here()),
                           "0_Data", 
                           "data_for_mixtures_analysis", 
                           str_c("Pooled_mixtures_datasets_", 
                                 mixtures_name, 
                                 "_hyper_g_prior.Rdata")))

rm(list = ls())
