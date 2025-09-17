
# 
# Sept 16, 2025

# Run R
/sc/arion/projects/CommonMind/leed62/conda/envs/dreamlet430/bin/R --vanilla
# Set path
.libPaths("/sc/arion/projects/CommonMind/leed62/conda/envs/dreamlet430/lib/R/library")

library(SingleCellExperiment)
library(zellkonverter)

folder = "/sc/arion/projects/psychAD/NPS-AD/public_release_0/" 
file = paste0(folder, "PsychAD_r0_Dec_28_2022.h5ad")
sce = readH5AD(file, use_hdf5=TRUE, verbose=TRUE)
assayNames(sce)[1] = "counts"
sce$Dx = factor(sce$Dx_AD, c('Control','AD'))


