# Set up CHS data_for_mixtures_analysis on HPC
# Note: run this code from "0_set up datasets for mixtures analysis.R"
# rm(list = ls())
source(here::here("0_project_setup", "!directories.R"))
source(here::here("0_project_setup", "!set_exposure_outcome_vars.R"))
source(here::here("0_project_setup", "!load_data.R"))
source(here::here("2_mixtures_MWAS",
                  "0_0_BHRMA.g_function.R"))

# Get all metabolite data in single dataframe
chs_metab_dat <- purrr::reduce(ftdata$chs, .f = left_join)

# Create numeric vars
chs_eo <- exposure_outcome$chs %>% 
  mutate(id = as.character(id), 
         sex.num = ifelse(sex == "Female",1,0),
         ses.num = recode(ses, "[1,2]" = 1, "(2,4]" = 2, "(4,9]" = 3))

# Join exposures and metabolites 
chs <- left_join(chs_eo, 
                 chs_metab_dat %>% mutate(id = as.character(id)), 
                 by = "id")

# PFAS
X.obs = chs[mixtures_components] %>%
  mutate(across(everything(), ~scale(log2(.))))


Y = chs %>%
  dplyr::select(colnames(chs_metab_dat)[2]:ncol(chs)) %>% # exclude outcome, leave only predictors
  scale(center = F, scale = F) 

dim(Y)

# Covariates
U = cbind.data.frame(age = as.numeric(chs$age), 
                     sex.num = as.numeric(chs$sex.num), 
                     bmi = as.numeric(chs$bmi), 
                     ses.num = as.numeric(chs$ses.num))


P = ncol(X.obs)
LOD = c(0.01, 0.01, 0.05, 0.01, 0.01, 0.01, 0.43, 0.01)
profiles = c(-1,1)*matrix(.5, nrow=2, ncol=P)


rm(list = setdiff(ls(), 
                  c("X.obs",
                    "Y",
                    "U",
                    "LOD",
                    "profiles", 
                    "mixtures_name",
                    "ridge.BDL.model", 
                    "BHRMA.g")))

save.image(
  file = fs::path(dirname(here::here()),
                  "0_Data", 
                  "data_for_mixtures_analysis",  
                  str_c("CHS_mixtures_datasets_", 
                        mixtures_name, 
                        "_hyper_g_prior.Rdata")))

rm(list = ls())
