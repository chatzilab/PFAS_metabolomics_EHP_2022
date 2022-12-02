# Compound ID to Compound Names
library(janitor)

# read in hand created key -----------------------------------
file_location <- fs::path(
  dir_data_local,
  "Supporting Files", 
  "Cpd id to name keys",
  "Mummichog cpd id to cpd name key with db matches.xlsx")

# Modify annotated data to get common compound names ------------------------
# read annotated cmpd data
annotated_fts_from_mum <- read_rds(
  fs::path(dir_data_local,
           "Supporting Files", 
           "mummichog_pw_ec_feature_key_cpd_id_only.rds"))


# Read in metaboanalyst key 
nm_conv_metaboanalyst <- readxl::read_xlsx(file_location, 
                                           sheet = "metaboanalyst", 
                                           na = "NA") %>% 
  janitor::clean_names() %>% 
  rename(chem_id=query) %>% 
  filter(comment == 1) # Comment = 1 for compounds with a match

# Read in metanet key
nm_conv_metanet <- readxl::read_xlsx(file_location, 
                                     sheet = "MetaNetX") %>% 
  janitor::clean_names() %>% 
  rename(chem_id = number_query)

# Read in mbrole name conversions
nm_conv_mbrole <- readxl::read_xlsx(file_location, 
                                    sheet = "mbrole") %>% 
  janitor::clean_names() %>% 
  rename(chem_id=input) 


# Modify mbrole 
mbrole_w <- nm_conv_mbrole %>% 
  group_by(chem_id, output_source) %>% 
  summarise(input_source = paste(unique(input_source),collapse = "; "), 
            output = str_c(unique(output), collapse = "; ")) %>% 
  tidylog::pivot_wider(names_from = "output_source", 
                       values_from = c("output") , 
                       values_fn = list) %>% 
  janitor::clean_names() %>%
  unnest(cas:lipid_maps) %>% 
  tidylog::select(-ymdb)


# Combine all annotations
name_cpd_key <- tidylog::full_join(nm_conv_metaboanalyst, 
                                   mbrole_w, by = "chem_id", 
                                   suffix = c("_metab", "_mbrole")) %>%
  tidylog::full_join(nm_conv_metanet,  
                     by = "chem_id", 
                     suffix = c("", "_metanet"))

# Merge matching compounds across databases 
name_cpd_key2 <- name_cpd_key %>% 
  select(-input_source) %>%  
  mutate(hmdb = if_else(is.na(hmdb_metab), hmdb_mbrole, hmdb_metab)) %>% 
  select(chem_id, match, name, hmdb, everything(), -hmdb_mbrole, -hmdb_metab) %>% 
  mutate(name_2 = if_else(is.na(match), name, match)) %>% 
  select(chem_id, match, hmdb, everything(), -name) %>% 
  janitor::remove_empty(which = c("rows", "cols")) %>%
  rename(matched_compound = chem_id, 
         met_name = name_2)

colnames(name_cpd_key2)

# Join Annotated fts from mumichog with names of compounds 
final_annotated_cpds <- tidylog::left_join(annotated_fts_from_mum,
                                           name_cpd_key2)

# Change any remaining missing values to the compound id
final_annotated_cpds <- final_annotated_cpds %>% 
  mutate(met_name = if_else(is.na(met_name), 
                            matched_compound, 
                            met_name))

# Save Data
write_rds(final_annotated_cpds,
          fs::path(dir_data_local,
                   "Supporting Files", 
                   "mummichog_pw_ec_feature_key_with_cpd_names.rds"))