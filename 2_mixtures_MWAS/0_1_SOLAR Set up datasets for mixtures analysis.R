# Set up SOLAR data_for_mixtures_analysis on HPC

source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!set_exposure_outcome_vars.R"))
source(here::here("0_project_setup", "!load_data.R"))

# Solar -------------------------------------------

# Get BHRMA Function
source(here::here("2_mixtures_MWAS",
                  "0_0_BHRMA.g_function.R"))

# Change ftdata from list to dataframe
sol_metab_dat <- purrr::reduce(ftdata$solar, .f = left_join)

# Create numeric Vars
solar_eo <- exposure_outcome$solar %>% 
  mutate(id = as.character(id), 
         sex.num = ifelse(sex == "Female",1,0),
         ses.num = recode(ses, "[3,11]" = 1, 
                          "(11,15.5]" = 2, 
                          "missing" = 2.5,
                          "(15.5,22]" = 3,  
                          "(22,63.5]" = 4),
         wave.num = ifelse(wave == "first wave", 1 ,2))


# Join exposures and metabolites 
solar <- left_join(solar_eo, 
                   sol_metab_dat %>% mutate(id = as.character(id)), 
                   by = "id")

# PFAS
X.obs = solar[mixtures_components] %>%
  mutate(across(everything(), ~scale(log2(.))))

# exclude outcome, leave only predictors:
Y = solar %>%
  dplyr::select(colnames(sol_metab_dat)[2]:ncol(solar)) %>% 
  scale(center = F, scale = F) 

# Covariates
U = cbind.data.frame(age = as.numeric(solar$age), 
                     sex.num = as.numeric(solar$sex.num), 
                     bmi = as.numeric(solar$bmi), 
                     tanner = as.numeric(solar$tanner),
                     ses.num = as.numeric(solar$ses.num), 
                     wave.num = as.numeric(solar$wave.num))


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
                           str_c("SOLAR_mixtures_datasets_", 
                                 mixtures_name, 
                                 "_hyper_g_prior.Rdata")))

rm(list = ls())
