## Seurat analysis

## 2019/4/19

#load packages

library(Seurat,lib.loc = "~/anaconda3/envs/R/lib/R/library/")
library(cowplot)
library(dplyr)

#requre(Seurat)

library(Matrix)
library(future)
library(ggplot2)
#read file and return seurat object

plan("multiprocess", workers = 30)
#cca analysis in dims from 20 to 50
options(future.globals.maxSize = 30 * 1024^3)

CreateSeurat <- function(data = data, project.name = project.name, min.cells = 0,
                         min.features = 0, cell.ids = cell.ids, conditions = conditions,
                         var = "stim")
{
  object <- CreateSeuratObject(counts = data, project = project.name, min.cells = min.cells,
                               min.features = min.features)
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object <- subset(object, subset = nFeature_RNA > 0 & percent.mt < 100)
  object <- RenameCells(object, add.cell.id = cell.ids)
  object$var <- conditions
  object$multi <- cell.ids
  object <- NormalizeData(object = object, verbose = FALSE)
  object <- FindVariableFeatures(object = object, selection.method = "vst",
                                 nfeatures = 2000)
  object <- CellCycleScoring(object, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)
  return(object)
}


ReNameIdents <- function(object = object, condition = "var"){
  Idents(object) <- condition
  return(object)
}

plotl <- function(object = object, ctrl = ctrl, stim = stim, title = title){
  p1 <- ggplot(object, aes(ctrl, stim)) + geom_point() + ggtitle(title)
  return(p1)
}

ReNameIdents <- function(object = object){
  object$celltype.stim <- paste(Idents(object), object$var, sep = "_")
  object$celltype <- Idents(object)
  Idents(object) <- "celltype.stim"
  return(object)
}

ReDefaultAssay <- function(object, assays = "RNA"){
  DefaultAssay(object) <- assays
  return(object)
}


###########################################################
#        text file containing information of samples      #
#   file_name project_name cell_ids condition resolution  #
###########################################################

sample <- read.csv("work.csv", stringsAsFactors = F, header = T)
#
#create a list for integrate
object.list = list()

for (i in 1:nrow(sample)){
  file_name <-  sample[i,"file_name"]
  project_name <-  sample[i,"project_name"]
  cell_ids <- sample[i,"cell_ids"]
  condition <- sample[i,"condition"]
  resolution <- sample[i,"resolution"]
  species <- sample[i,"species"]

  file <- sample[i,"name"]
  name <- paste0("sample_", i)
  print(file_name)
  if (!dir.exists(file_name)) {
    #only support RData file type with dgTMatrix class
#    data <- load(file_name)
#    load(file_name)
#    data <- get(data)
    data <- data.table::fread(file_name, data.table=F)
    rownames(data) <- data[,1]
    data <- data[,-1]
  }
  if (dir.exists(file_name)) {
    data <- Read10X(file_name)
  }
  #doublet.file <- paste0(cell_ids,".txt")
  #doublet.bool <- as.logical(read.delim(doublet.file)[,1])
  #data <- data[,!doublet.bool]
  assign(name, CreateSeurat(data = data, project.name = project_name,
                            cell.ids = cell_ids, conditions = condition))

  tmp.list <- list(get(name))
  object.list <- c(object.list, tmp.list)
  rm(tmp.list)
}

rm(data)
rm(list = grep(x = ls(), pattern = "sample_", value = T))

if (!dir.exists("results")) dir.create("results")
setwd("results")

#plan("multiprocess", workers = 20)
#cca analysis in dims from 20 to 50
#options(future.globals.maxSize = 10 * 1024^3)


reduction_method <- "cca" # cca or rpca, rlsi is not available now
for (i in c(40)){
  gc()
  #plan("multiprocess", workers = 40)


  if (reduction_method == "rpca") {
    # zjw f
    # normalize and identify variable features for each dataset independently
    object.list <- lapply(X = object.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })

    # select features that are repeatedly variable across datasets for integration run PCA on each
    # dataset using these features
    features <- SelectIntegrationFeatures(object.list = object.list)
    object.list <- lapply(X = object.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
    # zjw l
  }



  immune.anchors <- FindIntegrationAnchors(object.list = object.list,
                                           dims = 1:i, verbose = T,reduction = reduction_method)
  combined_name <- paste0("combined_",reduction_method,"_", i)
  assign(combined_name, IntegrateData(anchorset = immune.anchors, dims = 1:i))
  assign(DefaultAssay(get(combined_name)), "integrated")
#  assign(combined_name, CellCycleScoring(combined_cca_40, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes))
#  assign(combined_name, ScaleData(object = get(combined_name), verbose = FALSE, vars.to.regress=c("percent.mt", "S.Score", "G2M.Score")))
  assign(combined_name, ScaleData(object = get(combined_name), verbose = FALSE, vars.to.regress=c("percent.mt", "S.Score", "G2M.Score")))
  #assign(combined_name, SCTransform(get(combined_name), vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")))
  assign(combined_name, RunPCA(object = get(combined_name), npcs = 100, verbose = FALSE))
  gc()
  # plan("multiprocess", workers = 10)
  assign(combined_name, JackStraw(object = get(combined_name), num.replicate = 100,
                                  dims = 100))
  assign(combined_name, ScoreJackStraw(object = get(combined_name), dims = 1:100))
 #plan("multiprocess", workers = 20)
  p1 <- JackStrawPlot(object = get(combined_name), dims = 1:75)
  PC_Score <- table(p1[["data"]][["PC.Score"]]) %>% names() %>% as.data.frame() %>%
    tidyr::separate(col = 1, into = c("1", "2"), sep = ":")

  PC_choose <- which(as.numeric(PC_Score[,2]) > 0.00001)[1] - 1
  #PC_choose <- 75
  print(PC_choose)
  # add if the pc value is NA, choose a default value 30;
  if (is.na(PC_choose)){
     PC_choose <- 30
     print("Change PC choosed, now value is 30")
  }
  assign(combined_name, RunUMAP(object = get(combined_name), reduction = "pca", dims = 1:PC_choose))
  assign(combined_name, FindNeighbors(object = get(combined_name), reduction = "pca", dims = 1:PC_choose))
  assign(combined_name, FindClusters(get(combined_name), resolution = resolution))
  saveRDS(get(combined_name), file=paste0(combined_name,".RDS"))
  p1 <- DimPlot(object = get(combined_name), reduction = "umap", group.by = "multi")
  p2 <- DimPlot(object = get(combined_name), reduction = "umap", label = TRUE)
  gc()
  cluster_file_name <- paste0(file,"_", i, "_cell_cluster.png")
  #saveRDS(get(combined_name), file="all.RDS")
  #pdf(file = cluster_file_name, width = 9, height = 4)
  p <- plot_grid(p1, p2)
  ggplot2::ggsave(p, filename = cluster_file_name, width = 16, height = 8, limitsize = FALSE)


#  assign(DefaultAssay(get(combined_name)), "RNA")
  ReDefaultAssay(get(combined_name))
#   plan("multiprocess", workers = 8)
  assign(combined_name, ScaleData(object = get(combined_name), verbose = FALSE, vars.to.regress="percent.mt"))
#   plan("multiprocess", workers = 24)
  markers <- FindAllMarkers(object = get(combined_name), min.pct = 0.25, logfc.threshold = 0.25)

  csv_file_name <- paste0(file,"_", i, "_cell_marker.csv")
  write.csv(markers, file = csv_file_name)

  #plot heatmap
  assign(combined_name, ScaleData(get(combined_name)))
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  heatmap_file_name <- paste0(file,"_", i, "_cell_heatmap.png")
  # pdf(file = cluster_file_name, width = 9, height = 4)
  p <- DoHeatmap(object = get(combined_name), features = top10$gene) + NoLegend()

  ggplot2::ggsave(p, filename = heatmap_file_name, height = 0.1 * length(top10$gene), width = 20, limitsize = FALSE)

  #Identify differential expressed genes across conditions
  #tmp.object <- ReNameIdents(object = get(combined_name))
  #DEG <- data.frame()
  #library(clusterProfiler)
  #for (l in levels(Idents(get(combined_name)))){
  #   plan("multiprocess", workers = 20)
  #   tmp.DEG <- try(FindMarkers(tmp.object, ident.1 = paste0(l,"_STIM"),
  #                              ident.2 = paste0(l,"_CTRL"), verbose = FALSE), silent = F)
  #   #DEG <- FindMarkers(tmp.object, ident.1 = paste0(l,"_STIM"), ident.2 = paste0(l,"_CTRL"), verbose = FALSE)
  #   if (class(tmp.DEG) == "data.frame"){
  #     tmp.DEG$cluster <- l
  #     DEG <- rbind(DEG, tmp.DEG)
  #
  #     if (species == "human"){
  #       library(org.Hs.eg.db)
  #       dir.create(paste0(file, "_GO_annotation_of_", i))
  #       setwd(paste0(file, "_GO_annotation_of_", i))
  #       ego_ALL <- enrichGO(gene = rownames(tmp.DEG),
                             #universe = names(geneList), #背景基因集
  #                           OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
  #                           keyType = 'SYMBOL',
  #                           ont = "ALL", #也可以是 CC  BP  MF中的一种
  #                           pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
  #                           pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
  #                           qvalueCutoff = 0.05)
  #                          #readable = TRUE)
  #       GO_file_name <- paste("GO_of", i, l, ".png",sep = "_")
  #       p = barplot(ego_ALL)
  #       ggplot2::ggsave(p, file = GO_file_name)
  #       setwd("../")
  #    }
  #  }

  #name <- paste0("DEG_", i, "_", l, ".csv")
  #write.csv(DEG, file = name)

#    cluster <- paste0(file,"_", l, "_cells.png")
#    assign(cluster, subset(x = get(combined_name), idents = l))
#    assign(cluster, ReNameIdents(object = get(cluster)))
#    avg.t.cells <- log1p(x = AverageExpression(object = get(cluster), verbose = FALSE)$RNA)
#    avg.t.cells$gene <- rownames(x = avg.t.cells)
#
#    p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle(cluster)
#    avg.t.cells$FC <- (avg.t.cells[,sample[1,"condition"]] + 0.1) / (avg.t.cells[,sample[2,"condition"]] + 0.1)
#    top_n <- avg.t.cells %>% top_n(10, wt = FC)
#    avg.t.cells$FC <- (avg.t.cells[,sample[2,"condition"]] + 0.1) / (avg.t.cells[,sample[1,"condition"]] + 0.1)
#    btm_n <- avg.t.cells %>% top_n(10, wt = FC)
#    top_n <- rbind(top_n, btm_n)
#
#    p1 <- LabelPoints(plot = p1, points = top_n$gene, repel = TRUE)
#    filename <- paste0(file, i, "_cluster_", l, "_cell_type_gene.png")
#    ggplot2::ggsave(p1, filename = filename)
#}
  filenames <- paste0(file, "_", i, ".RData")
  save(list = c(combined_name, "markers"), file = filenames)

  #plot cell proportion
  data.to.plot <- get(combined_name)@meta.data %>%
    group_by(var, seurat_clusters) %>% count()
  cell_counts <- get(combined_name)@meta.data %>% group_by(var) %>% count()
  data.to.plot <- dplyr::left_join(data.to.plot, cell_counts, "var")
  data.to.plot$freq <- data.to.plot$n.x/data.to.plot$n.y
  p <- ggplot2::ggplot(data.to.plot) +
    geom_bar(aes(x = seurat_clusters, y = freq, fill = var),
             stat = "identity", position = "dodge")  # + scale_fill_manual(values=c("blue", "red"))
  cell_proportion <- paste0(file,"_", i,"_cell_proportion.png")
  ggplot2::ggsave(p, file = cell_proportion, limitsize = FALSE)
  ls()
  rm(combined_name)
  gc()
}
