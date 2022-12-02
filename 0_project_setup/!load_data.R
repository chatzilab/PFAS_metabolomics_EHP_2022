library(janitor)
library(jag2)

# Load Metabolomics Feature Tables --------------------------------------
ftdata <- read_rds(fs::path(dir_data, 
                            "sol_chs_batch_cor_scaled_untargeted_fts.rds"))

# Obtain Feature metadata 
ft_metadata <- ftdata$solar %>% 
  modify(~data.frame(feature = colnames(.)[-1])) %>% 
  bind_rows(.id = "mode")

# Load Exposure Outcome Data from drive  ------------------------
sol <- read_rds(fs::path(dir_data,
                         "SOLAR exposure outcome data HRE PFAS.rds"))

chs <- read_rds(fs::path(dir_data,
                         "CHS MetaAir exposure outcome data HRE PFAS.rds"))


# Rename Datasets and remove emory pfas
solar_exposure_outcome <- sol$baseline %>% 
  select(-contains("emory")); rm(sol)

chs_exposure_outcome <- chs$baseline %>% 
  select(-contains("emory")); rm(chs)


# Calculate additional exposure variables: OC Chemicals --------------------------
solar_exposure_outcome <- solar_exposure_outcome  %>%
  mutate(across(contains("detect"), 
                ~if_else(str_detect(.,"non"), 0, 1)), 
         across(all_of(exposures_continuous), 
                ~log2(.),
                .names = "lg2_{col}")) %>% 
  mutate(pcb_num_detect = pcb_180_ngml_detect+pcb_153_ngml_detect+
           pcb_138_ngml_detect+pcb_118_ngml_detect, 
         pbde_num_detect = pbde_154_ngml_detect + 
           pbde_47_ngml_detect + 
           pbde_100_ngml_detect + 
           pbde_153_ngml_detect+
           pbde_85_ngml_detect, 
         ocs = dde_impute + 
           hexachlorobenzene_impute) %>% 
  filter(id %in% ftdata$solar$c18neg$id) 


#CHS 
chs_exposure_outcome <- chs_exposure_outcome %>% 
  mutate(across(contains("detect"), 
                ~if_else(str_detect(.,"non"), 0, 1)), 
         across(all_of(exposures_continuous), 
                ~log2(.),
                .names = "lg2_{col}")) %>% 
  mutate(pcb_num_detect = pcb_180_ngml_detect+pcb_153_ngml_detect+ pcb_138_ngml_detect+pcb_118_ngml_detect, 
         pbde_num_detect = pbde_154_ngml_detect + pbde_47_ngml_detect + 
           pbde_100_ngml_detect + pbde_153_ngml_detect+pbde_85_ngml_detect, 
         ocs = dde_impute + hexachlorobenzene_impute) 




#bind to list
exposure_outcome = list(solar = solar_exposure_outcome,
                        chs = chs_exposure_outcome)


# Modify exposure outcome: change naming of vars with below lod as NA
exposure_outcome <- exposure_outcome %>% 
  modify(~rename_with(.,
                      .cols = contains("conc_below_lod_na_"), 
                      ~str_c(., "_w_na") %>% 
                        str_remove_all("conc_below_lod_na_")))

# Clean Environment
rm(chs_exposure_outcome, solar_exposure_outcome)


# Load PFAS LOD Data 
lod <- read_csv(
  fs::path(dir_data, 
           "hre_pfas_lod.csv"))

