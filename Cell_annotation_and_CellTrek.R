#Untitled1

options(stringsAsFactors = FALSE)
options(scipen = 100)

.libPaths()

library(MASS)
library(Seurat,lib.loc = "~/anaconda3/envs/R/lib/R/library/")
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(VennDiagram)
library(Rmisc)
library(ggpubr)
library(BuenColors)
library(pheatmap)
library(RColorBrewer)
library(CellChat)
library(clusterProfiler)

# InstallData("stxBrain")

scRNAseq <- readRDS("/home/zjw/zjw/20220419ki_adipose/sc.rds")
snRNAseq <- readRDS("/home/zjw/zjw/20220419ki_adipose/sn.rds")

scRNAseq_meta <- scRNAseq@meta.data
snRNAseq_meta <- snRNAseq@meta.data

"B cell" "T cell" "Vasculature" "Macrophage" "Endothelial" "Adipocyte" "Fibroblast" "Fibro/adipogenic progenitor"

# scRNAseq_meta$cluster_anno[scRNAseq_meta$cluster_anno=="FAP"] <- "Fibro/adipogenic progenitor"
scRNAseq_meta$cluster_anno[scRNAseq_meta$cluster_anno=="Bcell"] <- "B cell"
scRNAseq_meta$cluster_anno[scRNAseq_meta$cluster_anno=="T"] <- "T cell"

snRNAseq_meta$cluster_anno[snRNAseq_meta$cluster_anno=="Fibroblasts"] <- "Fibroblast"
snRNAseq_meta$cluster_anno[snRNAseq_meta$cluster_anno=="Adipocytes"] <- "Adipocyte"




temp <- rbind(data.frame(cluster_anno=scRNAseq_meta$cluster_anno),
              data.frame(cluster_anno=snRNAseq_meta$cluster_anno))


scsncomb <- readRDS("/home/zjw/zjw/20220419ki_adipose/results/combined_rpca_40.RDS")

scsncomb_meta <- scsncomb@meta.data

DimPlot(scsncomb,split.by = "var")

scsncomb_meta$barcode <- rownames(scsncomb_meta)

# scsncomb_meta <- merge(scsncomb_meta,temp,by="barcode")
scsncomb_meta <- cbind(scsncomb_meta,temp)


scsncomb@meta.data <- scsncomb_meta

Idents(scsncomb) <- "cluster_anno"

scsncomb@active.ident <- factor(scsncomb@active.ident,levels = c("Adipocyte","Fibroblast","FAP","Vasculature","Endothelial","Macrophage","T cell","B cell"))

DimPlot(scsncomb,label = T)+
  scale_color_manual(values=jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","CD8","GMP-B")))

DimPlot(scsncomb,split.by = "var",group.by = "var")


scsncomb@active.ident <- factor(scsncomb@active.ident,levels = c("B cell","T cell","Macrophage","Endothelial","Vasculature","FAP","Fibroblast","Adipocyte"))


DotPlot(scsncomb,features = c("ADIPOQ","LAMA2","PDGFRA","STEAP4","JAM2","C1QA","CD3D","CD79A"),assay = "RNA") +
  theme(axis.text.x = element_text(angle = 315, hjust = 0,vjust = 1)) +
  labs(x="",y="")

scsncomb@active.ident <- factor(scsncomb@active.ident,levels = c("Adipocyte","Fibroblast","FAP","Vasculature","Endothelial","Macrophage","T cell","B cell"))



Idents(scRNAseq) <- "cluster_anno"
DimPlot(scRNAseq,label = T)


Idents(snRNAseq) <- "cluster_anno"
DimPlot(snRNAseq,label = T)




mk <- FindAllMarkers(object = scsncomb,only.pos = T,assay = "RNA")
mk$cluster <- as.character(mk$cluster)

ave_exp <- AverageExpression(scsncomb,assays = "RNA")[[1]]

gene_list <- c()
for (i in unique(mk$cluster)) {
  temp <- mk[mk$cluster==i,]
  temp <- temp[order(temp$avg_logFC,decreasing = T),]
  gene_list <- c(gene_list,temp[1:100,]$gene)
}


df <- data.frame()

for (i in gene_list) {
  df <- rbind(df,
              ave_exp[rownames(ave_exp)==i,])
}
colnames(df) <- colnames(ave_exp)





# Fig. 1c
pheatmap(df,
         scale = "row",
         cluster_cols = F,
         cluster_rows = F,
         #col=c(colorRampPalette(as.character(jdb_palette("solar_flare")))(250)),
         col=rev(colorRampPalette(brewer.pal(11, "RdBu"))(200)),
         show_rownames = F,
         treeheight_row = F,
         treeheight_col = 100,
         border = F,angle_col = 315)


all_term <- data.frame()
for (i in unique(mk$cluster)) {
  print(i)
  temp <- mk[mk$cluster==i,]
  temp <- temp[order(temp$avg_logFC,decreasing = T),]
  
  temp <- temp$gene[1:100]
  go_onto <- enrichGO(temp, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
  #go_onto <- simplify(go_onto)
  a <- data.frame(go_onto)
  all_term <- rbind(all_term,
                    data.frame(a,cluster=i))
}

all_term <- all_term[!all_term$Description %in% all_term$Description[duplicated(all_term$Description)],]



go_term <- data.frame()
for (i in unique(all_term$cluster)) {
  temp <- all_term[all_term$cluster==i,]
  temp <- temp[order(temp$p.adjust,decreasing = F),]
  go_term <- rbind(go_term,
                   temp[1:3,])
}

go_term$p.adjust <- (-log10(go_term$p.adjust))
go_term$cluster <- factor(go_term$cluster,levels = c("Adipocyte","Fibroblast","FAP","Vasculature","Endothelial","Macrophage","T cell","B cell"))
go_term <- go_term[rev(seq(1:nrow(go_term))),]



# go_term <- all_term[all_term$Description %in% c("neutrophil chemotaxis",
#                                                 "leukocyte chemotaxis",
#                                                 "MHC class II protein complex assembly",
#                                                 "extracellular matrix organization",
#                                                 "extracellular structure organization",
#                                                 "connective tissue development",
#                                                 "vascular process in circulatory system",
#                                                 "muscle contraction",
#                                                 "cell-substrate adhesion",
#                                                 "positive regulation of lymphocyte activation",
#                                                 "B cell mediated immunity",
#                                                 "lymphocyte mediated immunity"
# ),]

#有一个item富集了两次
go_term$Description[duplicated(go_term$Description)] <- paste0(go_term$Description[duplicated(go_term$Description)],"--",go_term$cluster[duplicated(go_term$Description)])


ggbarplot(go_term,x = "Description",y="p.adjust",fill = "cluster",color = "cluster") +
  theme(axis.text.x = element_text(face = "plain",size = 12,angle = 270, hjust = 0,vjust = 0.5),
        axis.text.y = element_text(face = "plain",size = 12,angle = 270, hjust = 0.5,vjust = 0.5),
        axis.title.y = element_text(angle = 270),
        legend.position="none") +
  labs(x="",y="-log10 adjusted p-value") +
  scale_fill_manual(values=jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","CD8","GMP-B")))+
  scale_color_manual(values=jdb_color_map(c("B","CD4","GMP","MPP","mono","Ery","CD8","GMP-B")))






## integrate sc(sn)RNA-seq and ST


st_s49 <- Load10X_Spatial(data.dir = "/home/zjw/zjw/20220419ki_adipose/s49",filename = "filtered_feature_bc_matrix.h5")
# st_s49 <- readRDS("/home/zjw/zjw/20220419ki_adipose/st_s49x.rds")
st_s49$row <- st_s49@images$slice1@coordinates$row
st_s49$col <- st_s49@images$slice1@coordinates$col
st_s49$imagerow <- st_s49@images$slice1@coordinates$imagerow
st_s49$imagecol <- st_s49@images$slice1@coordinates$imagecol
dim(st_s49)
# st_s49 <- subset(st_s49, subset = row > 52)
plot1 <- VlnPlot(st_s49, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(st_s49, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

st_s49 <- NormalizeData(st_s49)

st_s49 <- ScaleData(st_s49)

st_s49 <- SCTransform(st_s49, assay = "Spatial", verbose = FALSE, variable.features.n = 1500)
st_s49 <- RunPCA(st_s49, assay = "SCT", verbose = FALSE)
st_s49 <- FindNeighbors(st_s49, reduction = "pca", dims = 1:30)
st_s49 <- FindClusters(st_s49, verbose = FALSE, resolution = 1) #0.1 海马区 #(0.195 5haima)
st_s49 <- RunUMAP(st_s49, reduction = "pca", dims = 1:30)


SpatialDimPlot(st_s49, label = TRUE, label.size = 3)
SpatialDimPlot(st_s49, label = TRUE, label.size = 3,cols = jdb_palette("lawhoops")) + NoLegend()



## We first co-embed ST and scRNA-seq datasets using traint
WAT_traint <- CellTrek::traint(st_data=st_s49, sc_data=scsncomb, sc_assay='RNA', cell_names='cluster_anno')


## We can check the co-embedding result to see if there is overlap between these two data modalities
DimPlot(WAT_traint, group.by = "type")
DimPlot(WAT_traint, group.by = "type",split.by = "type")

## After coembedding, we can chart single cells to their spatial locations. Here, we use the non-linear interpolation (intp = T, intp_lin=F) approach to augment the ST spots.

WAT_celltrek <- CellTrek::celltrek(st_sc_int=WAT_traint, int_assay='traint', sc_data=scsncomb, sc_assay = 'RNA',
                                   reduction='pca', intp=T, intp_pnt=5000, intp_lin=F, nPCs=30, ntree=1000,
                                   dist_thresh=0.55, top_spot=5, spot_n=5, repel_r=20, repel_iter=20, keep_model=T)$celltrek

WAT_celltrek$cluster_anno <- factor(WAT_celltrek$cluster_anno, levels=sort(unique(WAT_celltrek$cluster_anno)))


## After cell charting, we can interactively visualize the CellTrek result using celltrek_vis

CellTrek::celltrek_vis(WAT_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, cluster_anno:id_new),
                       WAT_celltrek@images$slice1@image, WAT_celltrek@images$slice1@scale.factors$lowres)



SpatialFeaturePlot(st_s49,features = c("CD3E","CD96","IL7R"),alpha = c(0.1,1),pt.size.factor = 3) ## T cell marker

SpatialFeaturePlot(st_s49,features = c("C1QA","CD19","CD14"),alpha = c(0.1,1),pt.size.factor = 3) ## T cell marker


st_s49 <- AddModuleScore(st_s49, features=list(setdiff(mk$gene[mk$cluster=="Adipocyte"],mk$gene[!mk$cluster=="Adipocyte"])),
                     name='Adipocyte_marker', nbin=10, ctrl=50, seed=42)
SpatialFeaturePlot(st_s49, "Adipocyte_marker1",alpha = c(0.1,1))


st_s49 <- AddModuleScore(st_s49, features=list(setdiff(mk$gene[mk$cluster=="Macrophage"],mk$gene[!mk$cluster=="Macrophage"])),
                     name='Macrophage_marker', nbin=10, ctrl=50, seed=42)
SpatialFeaturePlot(st_s49, "Macrophage_marker1",alpha = c(0.1,1))

st_s49 <- AddModuleScore(st_s49, features=list(setdiff(mk$gene[mk$cluster=="T cell"],mk$gene[!mk$cluster=="T cell"])),
                     name='T_cell_marker_gene_set_score', nbin=10, ctrl=50, seed=42)
SpatialFeaturePlot(st_s49, "T_cell_marker_gene_set_score1",alpha = c(0.3,1))

## Based on the CellTrek result, we can summarize the colocalization patterns between different cell types using SColoc module. Here, we are using glutamatergic neuron cell types as an example (it is recommended to remove some cell types with very few cells, e.g., n<20). We first subset the glutamatergic neuron cell types from our charting result.

glut_cell <- scsncomb_meta$cluster_anno %>% unique()
names(glut_cell) <- make.names(glut_cell)
WAT_celltrek_glut <- subset(WAT_celltrek, subset=cluster_anno %in% glut_cell)
WAT_celltrek_glut$cluster_anno <- factor(WAT_celltrek_glut$cluster_anno, levels=glut_cell)


## Then we can use scoloc module to perform colocalization analysis.

WAT_sgraph_KL <- CellTrek::scoloc(WAT_celltrek_glut, col_cell='cluster_anno', use_method='KL', eps=1e-50)


## We extract the minimum spanning tree (MST) result from the graph

WAT_sgraph_KL_mst_cons <- WAT_sgraph_KL$mst_cons
rownames(WAT_sgraph_KL_mst_cons) <- colnames(WAT_sgraph_KL_mst_cons) <- glut_cell[colnames(WAT_sgraph_KL_mst_cons)]


## We then extract the metadata (including cell types and their frequencies)

WAT_cell_class <- WAT_celltrek@meta.data %>% dplyr::select(id=cluster_anno) %>% unique
WAT_celltrek_count <- data.frame(freq = table(WAT_celltrek$cluster_anno))
WAT_cell_class_new <- merge(WAT_cell_class, WAT_celltrek_count, by.x ="id", by.y = "freq.Var1")

CellTrek::scoloc_vis(WAT_sgraph_KL_mst_cons, meta_data=WAT_cell_class)



## Spatial-weighted gene co-expression analysis within the cell type of interest

WAT_celltrek_Adipocyte <- subset(WAT_celltrek, subset=cluster_anno=='T cell')
# WAT_celltrek_Adipocyte@assays$RNA@scale.data <- matrix(NA, 1, 1)
# WAT_celltrek_Adipocyte$cluster <- gsub('Adipocyte IT VISp ', '', WAT_celltrek_Adipocyte$cluster)
DimPlot(WAT_celltrek_Adipocyte)


## We select top 2000 variable genes (exclude mitochondrial, ribosomal and high-zero genes)

WAT_celltrek_Adipocyte <- FindVariableFeatures(WAT_celltrek_Adipocyte)
vst_df <- WAT_celltrek_Adipocyte@assays$RNA@meta.features %>% data.frame %>% mutate(id=rownames(.))
nz_test <- apply(as.matrix(WAT_celltrek_Adipocyte[['RNA']]@data), 1, function(x) mean(x!=0)*100)
hz_gene <- names(nz_test)[nz_test<20]
mt_gene <- grep('^MT-', rownames(WAT_celltrek_Adipocyte), value=T)
rp_gene <- grep('^RP1|^RPS', rownames(WAT_celltrek_Adipocyte), value=T)
vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) %>% arrange(., -vst.variance.standardized)
feature_temp <- vst_df$id[1:2000]


## We use scoexp to do the spatial-weighted gene co-expression analysis.

WAT_celltrek_Adipocyte_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=WAT_celltrek_Adipocyte, assay='RNA', approach='cc', gene_select = feature_temp, sigm=140, avg_cor_min=.4, zero_cutoff=3, min_gen=40, max_gen=400)


## We can visualize the co-expression modules using heatmap.

WAT_celltrek_Adipocyte_k <- rbind(data.frame(gene=c(WAT_celltrek_Adipocyte_scoexp_res_cc$gs[[1]]), G='co-exp gene set 1'),
                                  data.frame(gene=c(WAT_celltrek_Adipocyte_scoexp_res_cc$gs[[2]]), G='co-exp gene set 2')) %>%
  magrittr::set_rownames(.$gene) %>% dplyr::select(-1)
pheatmap::pheatmap(WAT_celltrek_Adipocyte_scoexp_res_cc$wcor[rownames(WAT_celltrek_Adipocyte_k), rownames(WAT_celltrek_Adipocyte_k)],
                   clustering_method='ward.D2', annotation_row=WAT_celltrek_Adipocyte_k, show_rownames=F, show_colnames=F,
                   treeheight_row=10, treeheight_col=10, annotation_legend = T, fontsize=8,
                   color=viridis(10), main='T cell spatial co-expression')

go_gs1 <- enrichGO(c(WAT_celltrek_Adipocyte_scoexp_res_cc$gs[[1]]),
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

go_gs2 <- enrichGO(c(WAT_celltrek_Adipocyte_scoexp_res_cc$gs[[2]]),
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

enrichplot::dotplot(go_gs1)
enrichplot::dotplot(go_gs2)

WAT_celltrek_Adipocyte <- AddModuleScore(WAT_celltrek_Adipocyte, features=WAT_celltrek_Adipocyte_scoexp_res_cc$gs, name='gene_set_score_', nbin=10, ctrl=50, seed=42)
## First we look into the coexpression module based on the scRNA-seq embedding
FeaturePlot(WAT_celltrek_Adipocyte, grep('gene_set_score_', colnames(WAT_celltrek_Adipocyte@meta.data), value=T))



SpatialFeaturePlot(WAT_celltrek_Adipocyte, grep('gene_set_score_', colnames(WAT_celltrek_Adipocyte@meta.data), value=T))



st_s49 <- AddModuleScore(st_s49, features=WAT_celltrek_Adipocyte_scoexp_res_cc$gs, name='gene_set_score_', nbin=10, ctrl=50, seed=42)

FeaturePlot(st_s49, grep('gene_set_score_', colnames(WAT_celltrek_Adipocyte@meta.data), value=T))

SpatialFeaturePlot(st_s49, grep('gene_set_score_', colnames(st_s49@meta.data), value=T),alpha = c(0.1,1))
