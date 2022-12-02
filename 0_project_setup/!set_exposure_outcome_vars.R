# Set exposure outcome vars

# List PFAS with > 25% above LOD
exposures_continuous <- c("pfos", "pfhxs","pfhps",
                          "pfoa","pfna","pfda")

# Log transformed or categorical PFAS
exposures_for_analysis <- c("lg2_pfda",
                            "lg2_pfhps",
                            "lg2_pfhxs",
                            "lg2_pfna",
                            "lg2_pfoa",
                            "lg2_pfos") 

# LC-MS Modes
modes = c("c18pos","c18neg", "hilicpos", "hilicneg")

# Cohort Names
cohort = c("solar", "chs")