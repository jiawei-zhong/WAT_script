options(stringsAsFactors = FALSE)
options(scipen = 100)

# library(SPOTlight,lib.loc = "~/anaconda3/envs/R/lib/R/library/")
library(Seurat,lib.loc = "~/anaconda3/envs/R4/lib/R/library/")
# library(dplyr)
library(SeuratData)
library(data.table)
library(Matrix)
# library(future,lib.loc = "~/anaconda3/envs/R/lib/R/library/")

#plan("multiprocess", workers = 12)
#options(future.globals.maxSize = 10 * 1024^3)

st_s49 <- Load10X_Spatial(data.dir = "/home/zjw/zjw/20220419ki_adipose/s49",filename = "filtered_feature_bc_matrix.h5")
st_s49$row <- st_s49@images$slice1@coordinates$row
st_s49$col <- st_s49@images$slice1@coordinates$col
st_s49$imagerow <- st_s49@images$slice1@coordinates$imagerow
st_s49$imagecol <- st_s49@images$slice1@coordinates$imagecol
st_s49 <- SCTransform(st_s49, assay = "Spatial", verbose = FALSE)
st_s49 <- RunPCA(st_s49, verbose = FALSE)
st_s49 <- RunUMAP(st_s49, dims = 1:30, verbose = FALSE)
st_s49 <- FindNeighbors(st_s49, dims = 1:30, verbose = FALSE)
st_s49 <- FindClusters(st_s49, verbose = FALSE)

st_s50 <- Load10X_Spatial(data.dir = "/home/zjw/zjw/20220419ki_adipose/s50",filename = "filtered_feature_bc_matrix.h5")
st_s50$row <- st_s50@images$slice1@coordinates$row
st_s50$col <- st_s50@images$slice1@coordinates$col
st_s50$imagerow <- st_s50@images$slice1@coordinates$imagerow
st_s50$imagecol <- st_s50@images$slice1@coordinates$imagecol
st_s50 <- SCTransform(st_s50, assay = "Spatial", verbose = FALSE)
st_s50 <- RunPCA(st_s50, verbose = FALSE)
st_s50 <- RunUMAP(st_s50, dims = 1:30, verbose = FALSE)
st_s50 <- FindNeighbors(st_s50, dims = 1:30, verbose = FALSE)
st_s50 <- FindClusters(st_s50, verbose = FALSE)

st_s51 <- Load10X_Spatial(data.dir = "/home/zjw/zjw/20220419ki_adipose/s51",filename = "filtered_feature_bc_matrix.h5")
st_s51$row <- st_s51@images$slice1@coordinates$row
st_s51$col <- st_s51@images$slice1@coordinates$col
st_s51$imagerow <- st_s51@images$slice1@coordinates$imagerow
st_s51$imagecol <- st_s51@images$slice1@coordinates$imagecol
st_s51 <- SCTransform(st_s51, assay = "Spatial", verbose = FALSE)
st_s51 <- RunPCA(st_s51, verbose = FALSE)
st_s51 <- RunUMAP(st_s51, dims = 1:30, verbose = FALSE)
st_s51 <- FindNeighbors(st_s51, dims = 1:30, verbose = FALSE)
st_s51 <- FindClusters(st_s51, verbose = FALSE)

spotlight_s49 <- readRDS("/home/zjw/zjw/20220419ki_adipose/SPOTlight_s49.rds")
spotlight_s50 <- readRDS("/home/zjw/zjw/20220419ki_adipose/SPOTlight_s50.rds")
spotlight_s51 <- readRDS("/home/zjw/zjw/20220419ki_adipose/SPOTlight_s51.rds")

spotlight_s49 <- spotlight_s49[[2]]
spotlight_s49 <- spotlight_s49[,grep("res_ss",colnames(spotlight_s49),invert = T)]
colnames(spotlight_s49) <- paste0("SPOTlight_",colnames(spotlight_s49))
cell_types_all <- colnames(spotlight_s49)[which(colnames(spotlight_s49) != "res_ss")]
st_s49@meta.data <- cbind(st_s49@meta.data, spotlight_s49)

spotlight_s50 <- spotlight_s50[[2]]
spotlight_s50 <- spotlight_s50[,grep("res_ss",colnames(spotlight_s50),invert = T)]
colnames(spotlight_s50) <- paste0("SPOTlight_",colnames(spotlight_s50))
cell_types_all <- colnames(spotlight_s50)[which(colnames(spotlight_s50) != "res_ss")]
st_s50@meta.data <- cbind(st_s50@meta.data, spotlight_s50)

spotlight_s51 <- spotlight_s51[[2]]
spotlight_s51 <- spotlight_s51[,grep("res_ss",colnames(spotlight_s51),invert = T)]
colnames(spotlight_s51) <- paste0("SPOTlight_",colnames(spotlight_s51))
cell_types_all <- colnames(spotlight_s51)[which(colnames(spotlight_s51) != "res_ss")]
st_s51@meta.data <- cbind(st_s51@meta.data, spotlight_s51)


stereoscope_s49 <- read.table("/home/zjw/zjw/20220419ki_adipose/stereoscope_out/s49/st_s49_MTX/W.2022-04-24121117.973429.tsv",sep = "\t",row.names = 1,header = 1)
stereoscope_s50 <- read.table("/home/zjw/zjw/20220419ki_adipose/stereoscope_out/s50/st_s50_MTX/W.2022-04-24121117.860783.tsv",sep = "\t",row.names = 1,header = 1)
stereoscope_s51 <- read.table("/home/zjw/zjw/20220419ki_adipose/stereoscope_out/s51/st_s51_MTX/W.2022-04-24121117.850861.tsv",sep = "\t",row.names = 1,header = 1)

colnames(stereoscope_s49) <-paste0("stereoscope_",colnames(stereoscope_s49))
colnames(stereoscope_s50) <-paste0("stereoscope_",colnames(stereoscope_s50))
colnames(stereoscope_s51) <-paste0("stereoscope_",colnames(stereoscope_s51))

st_s49@meta.data <- cbind(st_s49@meta.data, stereoscope_s49)
st_s50@meta.data <- cbind(st_s50@meta.data, stereoscope_s50)
st_s51@meta.data <- cbind(st_s51@meta.data, stereoscope_s51)
