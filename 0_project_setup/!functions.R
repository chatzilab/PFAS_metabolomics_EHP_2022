# Functions

# Read Data from HPC -----------------------------------------------------
read_data_from_hpc <- function(file_path, n_col_in_df){
  
  # Read in data, without headers.
  df_for_colnames <- read.table(file_path,
                                sep = ",", 
                                na.strings = "",
                                as.is = TRUE, 
                                fill = TRUE,
                                header = FALSE)
  
  # Get temp col names. Col names are not indexed correctly though- they 
  # Need to be shifted to the right by 1, and we need to 
  row1 <- df_for_colnames[1,] %>% as.character(.)
  # Create column names for readining in new data
  col_names <- c("model_term",
                 row1[-length(row1)], 
                 paste0(row1[1:n_col_in_df], ".999999"))
  
  # Read in final data with headers
  suppressWarnings(
    df_original <- read.table(file_path,
                              sep = ",",col.names = col_names,
                              na.strings = "",
                              as.is = TRUE, fill = TRUE, header = TRUE) %>% 
      dplyr::select(-model_term))
  
  # Clean Col names
  df_clean_names <- df_original %>% 
    janitor::clean_names() %>% 
    rename_all(~str_replace(., "x2_5", "lcl") %>% 
                 str_replace(., "x97_5", "ucl") %>% 
                 str_replace(., "p_val", "p") %>% 
                 str_replace(., "var_names", "var"))
  
  # Get new column names in a dataframe
  new_colnames <- tibble(cnms = colnames(df_clean_names)) %>% 
    mutate(variable = str_split_fixed(cnms, "_", n = 2)[,1], 
           group = str_split_fixed(cnms, "_", n = 2)[,2] %>% 
             if_else(. == "", "0", .) %>% 
             as.numeric)
  
  # create a list of colnames by column group
  cnms_bygroup <- split(new_colnames, new_colnames$group)
  
  # create a list of sol result by column groups
  results_list <- cnms_bygroup %>% 
    modify(~dplyr::select(df_clean_names, all_of(.$cnms)))
  
  # rename all cols to match names, then cbind
  results_df <- results_list %>% 
    modify(~rename_all(., ~str_split_fixed(., "_", 2)[,1])) %>% 
    bind_rows() %>% 
    filter(!is.na(metabolite))
  
  # Separate beta and pips
  results_final <- results_df %>% 
    mutate(term = case_when(str_detect(var, "beta") ~ "beta", 
                            str_detect(var, "gamma") ~ "pip",
                            str_detect(var, "psi") ~ "beta", 
                            str_detect(var, "eta") ~ "eta", 
                            TRUE ~ var), 
           var = str_remove(var, ".beta") %>% 
             str_remove(".gamma") %>% 
             str_remove("eta.") %>%
             str_replace("psi", "mixture")) %>% 
    rename(exposure = var, 
           estimate = mean, 
           p_value = p) %>% 
    dplyr::select(exposure, term, everything())
  
  return(results_final)
}

# # Rename Compounds -----------------------
rename_pfas <- function(pfas_names, include_asterisk = FALSE, 
                        arrange_by_class = FALSE){
  x <- tibble(pfas = pfas_names)
  
  suppressWarnings(
    pfas2 <-  x %>%
      mutate(pfas2 = case_when(
        pfas == "pfhxs" ~ "PFHxS",
        pfas == "pfhps" ~ "PFHpS",
        pfas == "pfpes" ~ "PFPeS",
        pfas == "pfhpa" ~ "PFHpA",
        pfas == "nmefosaab" ~ "N-MeFOSAA-b†", 
        pfas == "pfuda" ~ "PFUnDA†",
        pfas == "pfds" ~ "PFDS†",
        pfas == "netfosaa" ~ "N-EtFOSAA†",
        pfas == "pfns" ~ "PFNS†",
        pfas == "pfbs" ~ "PFBS†",
        pfas == "x82fts" ~ "8:2 FTS†", 
        pfas == "pfhxa" ~ "PFHxA†", 
        pfas == "pfdoa" ~ "PFDoDA†",
        pfas == "Mixture effect" ~ "Mixture effect",
        TRUE ~ toupper(pfas)) %>% 
          as.factor() %>% 
          fct_relevel(., 
                      "PFOS", "PFOA", "PFHxS", "PFNA", "PFHpS","PFDA", "PFPeS", 
                      "PFHpA","N-MeFOSAA-b†","N-EtFOSAA†","PFDS†","PFBS†", 
                      "8:2 FTS†", "PFDoDA†", "PFUnDA†","PFNS†","PFHxA†",
                      "Mixture effect")) 
    )
    
    if(include_asterisk == TRUE){ 
      suppressWarnings(
        
       pfas2 <-  x %>%
        mutate(pfas2 = case_when(
          pfas == "pfhxs" ~ "PFHxS",
          pfas == "pfhps" ~ "PFHpS",
          pfas == "pfpes" ~ "PFPeS",
          pfas == "pfhpa" ~ "PFHpA",
          pfas == "nmefosaab" ~ "N-MeFOSAA-b*", 
          pfas == "pfuda" ~ "PFUnDA*",
          pfas == "pfds" ~ "PFDS*",
          pfas == "netfosaa" ~ "N-EtFOSAA*",
          pfas == "pfns" ~ "PFNS*",
          pfas == "pfbs" ~ "PFBS*",
          pfas == "x82fts" ~ "8:2 FTS*", 
          pfas == "pfhxa" ~ "PFHxA*", 
          pfas == "pfdoa" ~ "PFDoDA*", 
          pfas == "Mixture effect" ~ "Mixture effect",
          TRUE ~ toupper(pfas)) %>% 
            as.factor() %>% 
            fct_relevel(., 
                        "PFOS", "PFOA", "PFHxS", "PFNA", "PFHpS","PFDA", "PFPeS",
                        "PFHpA","N-MeFOSAA-b*","N-EtFOSAA*","PFDS*","PFBS*", 
                        "8:2 FTS*", "PFDoDA*", "PFUnDA*","PFNS*","PFHxA*", 
                        "Mixture effect")) )
    }
    
    if(arrange_by_class == TRUE){ 
      suppressWarnings(
        
      pfas2 <-  pfas2 %>% 
        # left_join(lod, by = "pfas") %>% 
        mutate(pfas2 = fct_relevel(pfas2, 
                                   "PFOS", 
                                   "PFHxS",
                                   "PFHpS",
                                   "PFPeS",
                                   "PFOA", 
                                   "PFNA", 
                                   "PFDA", 
                                   "PFHpA")))
    }
    
    return(pfas2$pfas2)
}


# Get duplicate mz_rt
dup_mzrt <-  function(x){
  x2 <- x %>% 
    group_by(name) %>% 
    summarise(n_dup = length(name), 
              matched_form = str_c(unique(matched_form), collapse = ";"), 
              met_name = str_c(unique(met_name), collapse = ";"), ) %>% 
    ungroup() %>% 
    filter(n_dup > 1)
  
  return(x2)
}
