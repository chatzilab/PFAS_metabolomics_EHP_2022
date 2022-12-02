# directory fs::paths for file architecture
# ignore file on github
library(fs)

# home directory for project
dir_home <- dirname(here::here())

# data folder
dir_data <- fs::path("G:", 
                     "My Drive", 
                     "SOL CHS PFAS Metabolomics", 
                     "0_Data_mirror_do_not_edit")

# results results folder
dir_results <- fs::path(dir_home, "2_Results")

# reports folder
dir_reports <- fs::path(dir_home, "3_Reports")

# Local Data folder
dir_data_local <- fs::path(dir_home, "0_Data")

# Mixtures Data folder
dir_data_mixtures <- fs::path(dir_data_local, "data_for_mixtures_analysis")

# Mixtures results
dir_results_mixtures <- fs::path(dir_results, "PFAS_mixtures")

# Mummichog results
dir_results_mum_mixtures <- fs::path(dir_results_mixtures, "mummichog")
