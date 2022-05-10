library(Seurat)
library(patchwork)
library(stringr)

#library(future) #needed for multicore
#plan("multiprocess", workers = 18)
#future.seed=TRUE
#https://vertesy.github.io/Seurat.multicore/
#use htop to monitor

T = TRUE
F = FALSE

single_dataset <- function()
{
    filename <- '/hdd1/clarisse2/R/Skin dataset_SB/GSE130973_filtered'
    data.data <- Read10X(filename)
    data <- CreateSeuratObject(counts = data.data, project = "skin", min.cells = 3, min.features = 200)
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    data <- NormalizeData(data)
    data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
    data <- RunPCA(data, features = VariableFeatures(object = data))
    data <- FindNeighbors(data, dims = 1:10)
    data <- FindClusters(data, resolution = 0.5)
    data <- RunUMAP(data, dims = 1:10)
    print(DimPlot(data, reduction = "umap"))
}

load_dataset <- function(filename,ident,type='10x')
{
    if(type=='10x') data.data <- Read10X(filename)
    else if(type == 'tabib') {
        data.data <- read.csv(filename)
        rownames(data.data) <- data.data[,1]
        data.data <- data.data[,-1]
    }
    else if(type == 'newcastle') {
        data.data <- read.csv(filename)
        rownames(data.data) <- data.data[,1]
        data.data <- data.data[,-1]
        data.data <- t(data.data)
    }
    data <- CreateSeuratObject(counts = data.data, project = ident, min.cells = 3, min.features = 200)
    if(type=='tabib') data$orig.ident <- ident
    data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
    data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    return(data)
}

if(F) {
    body_abdomen1 = load_dataset('/hdd1/clarisse2/R/cellranger/WS_SKN_KCL9369535/filtered_feature_bc_matrix','body_abdomen1')
    body_solebordo = load_dataset('/hdd1/clarisse2/R/Skin dataset_SB/GSE130973_filtered','body_solebordo')
    body_tabib = load_dataset('/hdd1/clarisse2/R/Skin dataset_Tabib/Tabib data/Control_Human_Skin_6People_Fibroblasts.csv','body_tabib','tabib')
    face_cheek1 = load_dataset('/hdd1/clarisse2/R/cellranger/WS_SKN_KCL9369531/filtered_feature_bc_matrix','face_cheek1')
    face_cheek2 = load_dataset('/hdd1/clarisse2/R/cellranger/WS_SKN_KCL9369534/filtered_feature_bc_matrix','face_cheek2')
    face_cheek3a = load_dataset('/hdd1/clarisse2/R/cellranger 2nd round/counts/cellranger310_count_36610_WS_SKN_KCL9369627_GRCh38-3_0_0/filtered_feature_bc_matrix','face_cheek3a')
    face_cheek3b = load_dataset('/hdd1/clarisse2/R/cellranger 2nd round/counts/cellranger310_count_36610_WS_SKN_KCL9369632_GRCh38-3_0_0/filtered_feature_bc_matrix','face_cheek3b')
    face_cheek4 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525739/filtered_feature_bc_matrix','face_cheek4')
    face_ear1 = load_dataset('/hdd1/clarisse2/R/cellranger/WS_SKN_KCL9369529/filtered_feature_bc_matrix','face_ear1')
    face_ear2 = load_dataset('/hdd1/clarisse2/R/cellranger 2nd round/counts/cellranger310_count_36610_WS_SKN_KCL9369625_GRCh38-3_0_0/filtered_feature_bc_matrix','face_ear2')
    face_forehead1 = load_dataset('/hdd1/clarisse2/R/cellranger/WS_SKN_KCL9369532/filtered_feature_bc_matrix','face_forehead1')
    face_forehead2 = load_dataset('/hdd1/clarisse2/R/cellranger/WS_SKN_KCL9369533/filtered_feature_bc_matrix','face_forehead2')
    face_forehead3 = load_dataset('/hdd1/clarisse2/R/cellranger 2nd round/counts/cellranger310_count_36610_WS_SKN_KCL9369629_GRCh38-3_0_0/filtered_feature_bc_matrix','face_forehead3')
    face_forehead4 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525738/filtered_feature_bc_matrix','face_forehead4')
    face_forehead5 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525741/filtered_feature_bc_matrix','face_forehead5')
    face_nose1 = load_dataset('/hdd1/clarisse2/R/cellranger/WS_SKN_KCL9369530/filtered_feature_bc_matrix','face_nose1')
    face_temple1 = load_dataset('/hdd1/clarisse2/R/cellranger 2nd round/counts/cellranger310_count_36610_WS_SKN_KCL9369628_GRCh38-3_0_0/filtered_feature_bc_matrix','face_temple1')
    face_temple2 = load_dataset('/hdd1/clarisse2/R/cellranger 2nd round/counts/cellranger310_count_36610_WS_SKN_KCL9369631_GRCh38-3_0_0/filtered_feature_bc_matrix','face_temple2')
   
bcc_cheek1 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525740/filtered_feature_bc_matrix','bcc_cheek1')
    bcc_cheek2 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525747/filtered_feature_bc_matrix','bcc_cheek2')
    bcc_ear1 = load_dataset('/hdd1/clarisse2/R/cellranger 2nd round/counts/cellranger310_count_36610_WS_SKN_KCL9369626_GRCh38-3_0_0/filtered_feature_bc_matrix','bcc_ear1')
    bcc_ear2 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525746/filtered_feature_bc_matrix','bcc_ear2')
    bcc_forehead1 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525742/filtered_feature_bc_matrix','bcc_forehead1')
    bcc_nose1 = load_dataset('/hdd1/clarisse2/R/cellranger 2nd round/counts/cellranger310_count_36610_WS_SKN_KCL9369630_GRCh38-3_0_0/filtered_feature_bc_matrix','bcc_nose1')
    bcc_nose2 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525744/filtered_feature_bc_matrix','bcc_nose2')
    bcc_temple1 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525743/filtered_feature_bc_matrix','bcc_temple1')
    bcc_temple2 = load_dataset('/hdd1/clarisse2/R/cell ranger 3rd round/cellranger.tar/cellranger/WS_SKN_KCL10525745/filtered_feature_bc_matrix','bcc_temple2')
}

standard_integration <- function()
{
    all_data.list <- c(body_abdomen1,body_tabib,body_solebordo)
    #all_data.list <- c(body_abdomen1,body_solebordo,body_tabib,face_cheek1,face_cheek2,face_cheek3a,face_cheek3b,face_cheek4,face_ear1,face_ear2,face_forehead1,face_forehead2,face_forehead3,face_forehead4,face_forehead5,face_nose1,face_temple1,face_temple2)
    all_data.list <- lapply(X = all_data.list, FUN = function(x) {
    x <- NormalizeData(x)
x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = all_data.list)
anchors <- FindIntegrationAnchors(object.list = all_data.list, anchor.features = features)
data_combined <- IntegrateData(anchorset = anchors)
    return(data_combined)
}

fast_integration <- function(all_data.list,reference)
{
    #all_data.list <- c(body_solebordo,face_cheek1,body_abdomen1,body_tabib,face_cheek2,face_cheek3a,face_cheek3b,face_cheek4,face_ear1,face_ear2,face_forehead1,face_forehead2,face_forehead3,face_forehead4,face_forehead5,face_nose1,face_temple1,face_temple2,bcc_cheek1,bcc_cheek2,bcc_ear1,bcc_ear2,bcc_forehead1,bcc_nose1,bcc_nose1,bcc_temple1,bcc_temple2)
    #reference = c(1,2)
all_data.list <- lapply(X = all_data.list, FUN = function(x) {
x <- NormalizeData(x, verbose = FALSE)
x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = all_data.list)
all_data.list <- lapply(X = all_data.list, FUN = function(x) {
x <- ScaleData(x, features = features, verbose = FALSE)
x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = all_data.list, reference = reference, reduction = "rpca",
dims = 1:50)
all_data.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
all_data.integrated <- ScaleData(all_data.integrated, verbose = FALSE)
all_data.integrated <- RunPCA(all_data.integrated, verbose = FALSE)
all_data.integrated <- RunUMAP(all_data.integrated, dims = 1:50)
DimPlot(all_data.integrated, group.by = "orig.ident")
    return(all_data.integrated)
}
plot <- function()
{
    if(T) {
        DefaultAssay(data_combined) <- "integrated"
        data_combined <- ScaleData(data_combined, verbose = FALSE)
        data_combined <- RunPCA(data_combined, npcs = 30, verbose = FALSE)
        data_combined <- RunUMAP(data_combined, reduction = "pca", dims = 1:30)
        data_combined <- FindNeighbors(data_combined, reduction = "pca", dims = 1:30)
        data_combined <- FindClusters(data_combined, resolution = 0.5)
    }
    if(T) {
        p1 <- DimPlot(data_combined, reduction = "umap", group.by = "orig.ident")
        print(p1)
    }
}
plot_markers <- function(data)
{
#print(FeaturePlot(data, features = c('LUM','CD3D','CLDN5','PROX1','PMEL','KRT14','TAGLN','CD68','TPSAB1','CCL5','CD79A','CLEC9A')))
print(FeaturePlot(data, features = c('DCN','LUM','APOE','CCL2','CXCL2','PTGDS','ASPN','WIF1','POSTN','WISP2','MFAP5','COL6A5','COL23A1','APCDD1','COL1A1')))
#print(FeaturePlot(data, features = c('WIF1','CD133','PROM1','SOX2')))
#print(FeaturePlot(data, features = c('ACTB')))
#print(FeaturePlot(data, features = c('WIF1','CCL20','CXCL1','COL1A1','CCL21','CXCL5','FABP4','PLA2G2A','CXCL3','CCL2','HBA2')))
#print(FeaturePlot(data, features = c('WIF1','VCAN','MZF1','MEF2C','ACTA2','GREM2','SM22','PRDM1','LEF1')))
#plot()
#print(FeaturePlot(data, features = c('LUM','WIF1')))
#print(FeaturePlot(data, features = c('CCL19','CXCL1','CCL21','PMEL','RNASE1','SELE')))
}
plot_sites <- function()
{
#print(DimPlot(data_combined, group.by = "orig.ident") )
#print(DimPlot(data_combined, reduction = "umap", group.by = "orig.ident"))
#print(DimPlot(data_combined, cols=c(), reduction = "umap", group.by = "orig.ident"))
#body <- subset(x = data_combined, subset = orig.ident %in% c("body_tabib","body_solebordo","body_abdomen1"))
#face <- subset(x = data_combined, subset = orig.ident %in% c('face_cheek1','face_cheek2','face_cheek3a','face_cheek3b','face_cheek4','face_ear1','face_ear2','face_forehead1','face_forehead2','face_forehead3','face_forehead4','face_forehead5','face_nose1','face_temple1','face_temple2'))
#print(DimPlot(body, reduction = "umap", group.by = "orig.ident", cols = c('red','red','red')))
#print(DimPlot(face, reduction = "umap", group.by = "orig.ident", cols = c('red','red','red','red','red','green','green','blue','blue','blue','blue','blue','yellow','orange','orange')))
#print(DimPlot(body, reduction = "umap", group.by = "orig.ident", cols = c('blue','blue','blue')))
#print(DimPlot(face, reduction = "umap", group.by = "orig.ident", cols = c('red','red','red','red','red','red','red','red','red','red','red','red','red','red','red')))
}
load <- function()
{
    non_nc_normal <- readRDS('/data/singlecell/ourdata/all_normal.rds')
    return(data_combined)
}
save <- function()
{
    saveRDS(data_combined,'/data/singlecell/bcc_and_normal.rds')
    saveRDS(integrated_fibroblasts,'/data/singlecell/integrated_fibroblasts_ncbody_base.rds')
#saveRDS(fibroblasts,'/data/singlecell/ourdata/fibroblasts.rds')
}
#data_combined <- standard_integration()
#data_combined <- fast_integration()
plot_clusters <- function (data)
{
    data <- FindNeighbors(data,dims=1:20)
    data <- FindClusters(data, resolution = 1.0)
    #v = rep('red',27)
    #v[5] = 'black'
    print(DimPlot(data, reduction = "umap"))
#markers <- FindMarkers(object = data, ident.1 = 5)
    return(data)
}
extract_cells <- function(data,ids)
{
    fibroblasts = c(1,8,10,12,15,23,30,33,34)
    vasc_endothelial = c(13)
    t_cells = c(0,2,3,9,11,24)
    pericytes = c(6,16)
    keratinocytes = c(18,19)
   
    selected = ids 
    #data_combined <- FindClusters(data_combined, resolution = 0.4)
    v = rep('red',39)
    for(x in selected) {
        v[x+1] = 'black'
    }
    #print(DimPlot(data, reduction = "umap", raster=FALSE, cols=v))
    extracted <- subset(data,idents=selected)
    DefaultAssay(extracted) <- "RNA"
    return(extracted)
}
plot_all_clusters <- function(data)
{
    n=length(unique(Idents(data)))
    for(i in 0:n+1) {
        v = rep('red',n)
        v[i] = 'black'
        print(i)
        print(DimPlot(data,  reduction = "umap",cols=v))
    }
}
make_new_idents <- function(data,prefix)
{
    v = c()
    for(col in colnames(data)) {
        s <- unlist(str_split(col,'-') )[3]
        #s <- paste('data_',s,sep='')
        s <- paste(prefix,s,sep='_')
        print(s)
        v <- c(v,s)
    }
    data$orig.ident <- as.factor(v)
    return(data)
}
#Reculster for example fibroblasts or T cells that have been extracted
recluster_celltype <- function(data)
{
    DefaultAssay(data) <- "integrated"
    #data <- NormalizeData(data)
    #data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(data)
    data <- ScaleData(data, features = all.genes)
    data <- RunPCA(data, features = VariableFeatures(object = data))
    data <- FindNeighbors(data, dims = 1:10)
    data <- FindClusters(data, resolution = 1.0)
    data <- RunUMAP(data, dims = 1:10)
    #print(DimPlot(data, reduction = "umap",group.by="orig.ident"))
    #print(DimPlot(data, reduction = "umap"))
    return(data)
}
#For each dataset
filter_ncells <- function(in_list,out_list,cutoff)
{
    for(x in in_list) {
        if(ncol(x$RNA)>cutoff) {
            out_list <- c(out_list,x)
        }
    }
    return(out_list)
}
#nc_fibroblasts <- readRDS('/data/singlecell/newcastle/nc_fibroblasts.rds')
#fibroblasts <- readRDS('/data/singlecell/fibroblasts.rds')
combine_face_with_newcastle <- function(face_cells,nc_cells)
{
    a <- SplitObject(face_cells,split.by='orig.ident')
    b <- SplitObject(nc_cells,split.by='orig.ident')
    combined_cells <- c()
    combined_cells <- filter_ncells(a,combined_cells,300)
    combined_cells <- filter_ncells(b,combined_cells,300)
    #integrated_cells <- fast_integration(combined_cells,integration_reference)
    #return(integrated_cells)
    return(combined_cells)
}
plot_fibroblasts_face_body <- function(data)
{
    cols <- rep('#0cB702',39)
    cols[1:2] <- 'red'
    cols[3:15] <- '#00a9ff'
    print(DimPlot(data, group.by = "orig.ident", cols =cols) )
}
plot_orig <- function(data)
{
    print(DimPlot(data, group.by = "orig.ident", raster=FALSE))
}
plot_face_vs_body <- function(data)
{
    cols <- rep('red',19)
    cols[1:3] <- 'black'
    print(DimPlot(data, group.by = "orig.ident", cols=cols) )
}
plot_body_site <- function(data)
{
    cols <- rep('#f8766d',18)
    cols[1:4] <- '#e68613' #body
    cols[4:9] <- '#0cb702' #cheek
    cols[9:11] <- '#c77cff' #ear
    cols[11:16] <- '#8494ff' #forehead
    cols[16:17] <- '#ff68a1' #nose
    cols[17:19] <- '#ed68ed' #temple
    print(DimPlot(data, group.by = "orig.ident", cols =cols) )
}
plot_pericytes_face_body <- function(pericytes)
{
    cols <- rep('#0cB702',28)
    cols[1:2] <- 'red'
    cols[10:28] <- '#00a9ff'
    print(DimPlot(pericytes, group.by = "orig.ident", cols =cols) )
}
plot_t_cells_face_body <- function(pericytes)
{
    cols <- rep('#0cB702',46)
    cols[1:2] <- 'red'
    cols[18:47] <- '#00a9ff'
    print(DimPlot(pericytes, group.by = "orig.ident", cols =cols,raster=FALSE) )
}
#plot_markers(integrated_fibroblasts)
#plot_face_vs_body()
#print(DimPlot(fibroblasts, reduction = "umap", group.by="orig.ident"))
#integrated_fibroblasts <- fast_integration(combined_fibroblasts,c(29,33))
#integrated_fibroblasts <- plot_clusters(integrated_fibroblasts)
#plot_all_clusters(non_nc_normal)
#print(DimPlot(non_nc_normal,  reduction = "umap", group.by="orig.ident"))
id_fibroblasts = c(1,8,10,12,15,23,30,33,34)
id_vasc_endothelial = c(13)
id_t_cells = c(0,2,3,9,11,24)
id_pericytes = c(6,16)
id_keratinocytes = c(18,19)
repeat_umap_subpop <- function()
{
    #
    #non_nc_normal <- readRDS('/data/singlecell/ourdata/all_normal.rds')
    #fibroblasts <- extract_cells(non_nc_normal,id_fibroblasts)
    #t_cells <- extract_cells(non_nc_normal,id_t_cells)
    #pericytes <- extract_cells(non_nc_normal,id_pericytes)
    #vasc_endothelial <- extract_cells(non_nc_normal,id_vasc_endothelial)
    #keratinocytes <- extract_cells(non_nc_normal,id_keratinocytes)
    fibroblasts <- recluster_celltype(fibroblasts)
    pericytes <- recluster_celltype(pericytes)
    t_cells <- recluster_celltype(t_cells)
    vasc_endothelial <- recluster_celltype(vasc_endothelial)
    keratinocytes <- recluster_celltype(keratinocytes)
}
#plot_face_vs_body(t_cells)
#plot_face_vs_body(t_cells)
#plot_face_vs_body(data)
#print(DimPlot(data, group.by = "orig.ident"))
#data <- read.csv('/data/singlecell/newcastle/nc_pericytes.csv')
#pericytes_nc <- make_new_idents(pericytes_nc,'nc_pericytes')
#combined_pericytes <- combine_face_with_newcastle(pericytes,nc_pericytes)
#integrated_pericytes <- fast_integration(combined_pericytes,c(1,6))
#plot_pericytes_face_body(integrated_pericytes)
#print(DimPlot(integrated_pericytes, group.by = "orig.ident") )
import_t_cells_nc <- function()
{
    newcastle_t_cells1 <- load_dataset('/data/singlecell/newcastle/nc_t_cells1.csv','nc_t_cells','newcastle')
    newcastle_t_cells2 <- load_dataset('/data/singlecell/newcastle/nc_t_cells2.csv','nc_t_cells','newcastle')
    newcastle_t_cells3 <- load_dataset('/data/singlecell/newcastle/nc_t_cells3.csv','nc_t_cells','newcastle')
    newcastle_t_cells4 <- load_dataset('/data/singlecell/newcastle/nc_t_cells4.csv','nc_t_cells','newcastle')
    newcastle_t_cells5 <- load_dataset('/data/singlecell/newcastle/nc_t_cells5.csv','nc_t_cells','newcastle')
    newcastle_t_cells6 <- load_dataset('/data/singlecell/newcastle/nc_t_cells6.csv','nc_t_cells','newcastle')
    newcastle_t_cells7 <- load_dataset('/data/singlecell/newcastle/nc_t_cells7.csv','nc_t_cells','newcastle')
    newcastle_t_cells8 <- load_dataset('/data/singlecell/newcastle/nc_t_cells8.csv','nc_t_cells','newcastle')
    newcastle_t_cells9 <- load_dataset('/data/singlecell/newcastle/nc_t_cells9.csv','nc_t_cells','newcastle')
    t_cells_nc <- merge(newcastle_t_cells1,c(newcastle_t_cells2,newcastle_t_cells3,newcastle_t_cells4,newcastle_t_cells5,newcastle_t_cells6,newcastle_t_cells7,newcastle_t_cells8,newcastle_t_cells9))
    t_cells_nc <- make_new_idents(t_cells_nc,'nc_t_cells')
    saveRDS(t_cells_nc,'/data/singlecell/newcastle/t_cells.rds')
}
#combined_t_cells <- combine_face_with_newcastle(t_cells,t_cells_nc)
#integrated_t_cells <- fast_integration(combined_t_cells,c(10,35))
#plot_orig(integrated_t_cells)
#plot_t_cells_face_body(integrated_t_cells)
#saveRDS(integrated_t_cells,'/data/singlecell/integrated_t_cells.rds')
#print(DimPlot(integrated_t_cells, group.by = "orig.ident", cols =cols,raster=FALSE) )
export_integrated_data <- function(integrated_data,filename_prefix)
{
    integrated <- GetAssayData(object = integrated_data, assay = "integrated", slot = "data")- GetAssayData(object = integrated_data, assay = "integrated", slot = "data")
    #integrated <- GetAssayData(object = integrated_data, assay = "integrated", slot = "data")-
    #file.remove(paste(filename_prefix,'_idents.txt',sep=''))
    lapply(integrated_data$orig.ident, write, paste(filename_prefix,'_idents.txt',sep=''), append=TRUE)
    write.csv(integrated,paste(filename_prefix,'.csv',sep=''))
    return(integrated)
}
#integrated_t_cells <- readRDS('/data/singlecell/integrated_t_cells.rds')
#integrated <- export_integrated_data(integrated_t_cells,'/data/singlecell/integrated_t_cells')
#write.csv(integrated,'/data/singlecell/integrated_t_cells.csv',col.names=colnames(integrated))
#write.table(integrated,'/data/singlecell/integrated_t_cells.csv',col.names=colnames(integrated),sep=",")
 