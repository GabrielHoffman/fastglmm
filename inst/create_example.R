
# 
# Sept 16, 2025

# Run R
/sc/arion/projects/CommonMind/leed62/conda/envs/dreamlet430/bin/R --vanilla
# Set path
.libPaths("/sc/arion/projects/CommonMind/leed62/conda/envs/dreamlet430/lib/R/library")

library(SingleCellExperiment)
library(tidyverse)
library(zellkonverter)
library(matrixStats)

folder = "/sc/arion/projects/psychAD/NPS-AD/public_release_0/" 
file = paste0(folder, "PsychAD_r0_Dec_28_2022.h5ad")
sce = readH5AD(file, use_hdf5=TRUE, verbose=TRUE, reader="R")
assayNames(sce)[1] = "counts"
sce$Dx = factor(sce$Dx_AD, c('Control','AD'))

libSize = colSums2(counts(sce))


info = with(colData(sce), 
        tibble(CellID = colnames(sce),
          SubID, 
          Sex, 
          CellType = subclass, 
          Dx, 
          Age, 
          libSize)) %>%
  filter(CellType == 'Micro_PVM') %>%
  select(-CellType)

info$PTPRG = counts(sce)["PTPRG",info$CellID]

PsychAD = info %>%
  select(-CellID)

save(PsychAD, file="~/PsychAD.RData")



