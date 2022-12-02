# script for metabolomics mixtures analysis
# Jesse Goodrich 111221

# Load Libraries 
library(rjags)
library(R2jags)
library(pbdMPI)

init()

#Set working directory
setwd("/project/dconti_624/Users/jagoodri/chs")

# Load Data
load("chs_mixtures_datasets_w_09.Rdata")

# Get number of metabolites
n_met <- ncol(Y)

# Model
model <- function(i){
  output <- BHRMA.g(X=X.obs,
                    Y=Y[,i],
                    U=U,
                    LOD=LOD,
                    profiles=profiles)
  
  # Name Metabolite
  output$metabolite = colnames(Y)[i]
  #Return Output
  return(output)
}

# Run model
coefs <- pbdLapply(1:n_met, model, pbd.mode = "spmd")

# Save results
comm.write.csv(coefs,  file = "results/chs_pfas_mixtures_mwas_w_09.csv")

# message(paste("SUCCESS from rank", comm.rank()))

finalize()

slurmR::Slurm_clean(coefs)