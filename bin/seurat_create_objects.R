#!/usr/bin/env Rscript

library(getopt)
library(stringr)
library(Seurat)

# Init
system('type R')

spec = matrix(c(
  'cores', 'c', 1, "integer",
  'id', 'a', 1, "character",
  'output', 'b', 1, "character",
  'folder', 'd', 1, "character",
  'mincells', 'e', 1, "integer",
  'minfeatures', 'f', 1, "integer"
), byrow=TRUE, ncol=4)
opt = getopt(spec)

output_path <- paste0(".", "/", opt$output, ".rds")

# Load 10x object from folder
data.10x <- Read10X(opt$folder)
data <- CreateSeuratObject(counts = data.10x, project = opt$id, min.cells = opt$mincells, min.features = opt$minfeatures)

saveRDS(data, file = output_path)
