
#!/usr/bin/env Rscript

# Define arguments for Rscript
library(getopt)
library(Seurat)
library(STACAS)
library(future)
library(tidyverse)

spec = matrix(c(
  'runtype', 'l', 1, "character",
  'cores', 'c', 1, "integer",
  'id', 'a', 1, "character"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

# Set paths and load data
{
  if(length(commandArgs(trailingOnly = TRUE)) == 0){
    cat('No command line arguments provided, paths are set for running interactively in Rstudio server\n')
    
    plot_path = "./output/seurat/01_preprocessing/plots/"
    rds_path = "./output/seurat/01_preprocessing/rds_files/"
    data_path = "./output/seurat/01_preprocessing/cellranger/count/filtered_feature_bc_matrix"
    ncores = 4
    
  } else if (opt$runtype == "nextflow"){
    cat('pipeline running through Nextflow\n')
    
    plot_path = "./plots/"
    rds_path = "./rds_files/"
    data_path = "./input/rds_files/"
    ncores = opt$cores
    
    # Multi-core when running from command line
    plan("multicore", workers = ncores)
    options(future.globals.maxSize = 64* 1024^3)

  } else {
    stop("--runtype must be set to 'nextflow'")
  }
  
  cat(paste0("script ran with ", ncores, " cores\n"))
  dir.create(plot_path, recursive = T)
  dir.create(rds_path, recursive = T)
}

seurat_split <- readRDS(list.files(data_path, full.names = TRUE))

# Log normalize data and find variable features
seurat_split <- lapply(seurat_split, function(x) {
  NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Get array of all genes across all datasets in order to integrate using all features
all_features <- lapply(seurat_split, row.names) %>% Reduce(intersect, .)
# Find anchors used for integration
integration_data <- Run.STACAS(seurat_split, dims = 1:20, anchor.features = 500)
# Integrate data
integration_data <- IntegrateData(anchorset = integration_data, dims = 1:20, features.to.integrate = all_features)


# unmodified data still resides in the 'RNA' assay
# Scale both RNA and integrated assays. RNA used for heatmaps and DEA, integration for clustering
DefaultAssay(integration_data) <- "RNA"
integration_data <- ScaleData(integration_data, features = rownames(integration_data), vars.to.regress = "percent.mt")

DefaultAssay(integration_data) <- "integrated"
integration_data <- ScaleData(integration_data, features = rownames(integration_data), vars.to.regress = "percent.mt")

# Save RDS after integration
saveRDS(integration_data, paste0(rds_path, "integration_data.RDS"), compress = FALSE)
