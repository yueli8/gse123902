#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("scran")

library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(clustree)
library(dplyr)
library(Nebulosa)
setwd("~/gse123902")

#primary tumor samples
#input, CreateSeuratObject, filter, save
a<-"GSM3516662_MSK_LX653_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
#构建Seurat对象，这里会有个初筛，保证所有基因在至少3个细胞中表达（0.1%细胞数），保证每个细胞至少能检测到200个基因。
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#用subset函数，质控：筛选检测到基因数目超过2500或低于200的细胞，单个细胞中线粒体基因数目占比超过>5%
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
#数据标准化：默认使用数据标准化方法是LogNormalize, 每个细胞总的表达量都标准化到10000，然后log取对数
 pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#变化基因鉴定：鉴定在细胞间表达高度变化的基因，后续研究需要集中于这部分基因，首先计算每一个基因的均值和方差，并且直接模拟其关系。默认返回2000个基因。
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#数据缩放 线性转换缩放数据，ScaleData()函数可以实现此功能。最终每个基因均值为0，方差为1。结果存放于pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
#设置参数features是因为ScaleData默认处理前面鉴定的差异基因。这一步怎么做都不会影响到后续pca和聚类，但是会影响做热图
pbmc <- ScaleData(pbmc, features = all.genes)
#移除影响方差的因素
#pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")
#对缩放后的数据进行PCA分析，默认使用前面鉴定表达变化大的基因。使用features参数可以重新定义数据集。
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#Determine the ‘dimensionality’ of the dataset
pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
#確定下面的dims
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "p62.rds")

a<-"GSM3516663_MSK_LX661_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "p63.rds")

a<-"GSM3516665_MSK_LX675_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "p65.rds")

a<-"GSM3516667_MSK_LX676_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "p67.rds")

a<-"GSM3516669_MSK_LX679_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "p69.rds")

a<-"GSM3516670_MSK_LX680_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "p70.rds")

a<-"GSM3516672_MSK_LX682_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "p72.rds")

a<-"GSM3516674_MSK_LX684_PRIMARY_TUMOUR_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "p74.rds")

p62<-readRDS(file="p62.rds")
p63<-readRDS(file="p63.rds")
p65<-readRDS(file="p65.rds")
p67<-readRDS(file="p67.rds")
p69<-readRDS(file="p69.rds")
p70<-readRDS(file="p70.rds")
p72<-readRDS(file="p72.rds")
p74<-readRDS(file="p74.rds")

#setup tech and celltype
p62<-RenameCells(p62,add.cell.id="p62",for.merge=T)
p62@meta.data$tech<-"primary"
p62@meta.data$celltype<-"primary_62"

p63<-RenameCells(p63,add.cell.id="p63",for.merge=T)
p63@meta.data$tech<-"primary"
p63@meta.data$celltype<-"primary_63"

p65<-RenameCells(p65,add.cell.id="p65",for.merge=T)
p65@meta.data$tech<-"primary"
p65@meta.data$celltype<-"primary_65"

p67<-RenameCells(p67,add.cell.id="p67",for.merge=T)
p67@meta.data$tech<-"primary"
p67@meta.data$celltype<-"primary_67"

p69<-RenameCells(p69,add.cell.id="p69",for.merge=T)
p69@meta.data$tech<-"primary"
p69@meta.data$celltype<-"primary_69"

p70<-RenameCells(p70,add.cell.id="p70",for.merge=T)
p70@meta.data$tech<-"primary"
p70@meta.data$celltype<-"primary_70"

p72<-RenameCells(p72,add.cell.id="p72",for.merge=T)
p72@meta.data$tech<-"primary"
p72@meta.data$celltype<-"primary_72"

p74<-RenameCells(p74,add.cell.id="p74",for.merge=T)
p74@meta.data$tech<-"primary"
p74@meta.data$celltype<-"primary_74"

#merge
p62_63<-merge(p62,p63)
p65_67<-merge(p65,p67)
p69_70<-merge(p69,p70)
p72_74<-merge(p72,p74)

p62_63_65_67<-merge(p62_63,p65_67)
p69_70_72_74<-merge(p69_70,p72_74)
p8<-merge(p62_63_65_67,p69_70_72_74)
saveRDS(p8,file="p6_before_integrate.rds")
mpn<-readRDS(file="p6_before_integrate.rds")

#before integrate
mpn[["percent.mt"]] <- PercentageFeatureSet(mpn, pattern = "^Mt-")
VlnPlot(mpn, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = mpn, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
plot_grid(p1,p2)

#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(
    pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
    verbose = FALSE)
}
reference.list <- pancreas.list[c("primary_62","primary_63","primary_65","primary_67",
                                  "primary_69","primary_70","primary_72","primary_74")]

#setup k.anchor and k.filter correctly
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list)
#看a<-pancreas.anchors@object.list,挑選最小的數字作爲k.weight = 79.
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, k.weight = 79)#留心k.weight,默認100.
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "p6_after_integrated.rds")

hms_individual_integrated<-readRDS(file="p6_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))

#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE)
saveRDS(hms_cluster, file = "p8_cluster_test_1.2.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.001, 
                                logfc.threshold = 0.001
)
write.table(scRNA.markers,file="p8_cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="p8_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="p8_genes_1.2.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="p8_cluster_test_1.2.rds")
DimPlot(hms_cluster, label=TRUE,reduction = "umap")

new.cluster.ids <- c("Effector_cd4","B","nkt","Effector_cd8","Effector_cd8","Effector_cd8","macrophage",
                    "Thelper","Thelper","nk","Thelper","B","progenitor","epithelial","nk","mast",
                    "macrophage","macrophage","progenitor","epithelial","endothelial","B","B","B","progenitor",
                    "B","Pre_exhasted_cd8")
                   
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "p8_cluster_id.rds")

setwd("~/gse123902/CytoTRACE")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="p8_cluster_id.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4,subsamplesize = 1000)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p

#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
"0"="Effector_cd4","1"="B","2"="nkt","3"="Effector_cd8","4"="Effector_cd8","5"="Effector_cd8","6"="macrophage",
"7"="Thelper","8"="Thelper","9"="nk","10"="Thelper","11"="B","12"="progenitor","13"="epithelial","14"="nk",
"15"="mast","16"="macrophage","17"="macrophage","18"="progenitor","19"="epithelial","20"="endothelial",
"21"="B","22"="B","23"="B","24"="progenitor","25"="B","26"="Pre_exhasted_cd8" ) 
colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)
plotCytoTRACE(result01, phenotype = phe)

#metastasis samples
a<-"GSM3516664_MSK_LX666_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "m64.rds")

a<-"GSM3516668_MSK_LX255B_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "m68.rds")

a<-"GSM3516671_MSK_LX681_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "m71.rds")

a<-"GSM3516677_MSK_LX699_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "m77.rds")

a<-"GSM3516678_MSK_LX701_METASTASIS_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#pbmc <- JackStraw(pbmc, num.replicate = 100,dims = 40)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:40)
#JackStrawPlot(pbmc, dims = 1:40)
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "m78.rds")

#input readRDS
m64<-readRDS(file="m64.rds")
m68<-readRDS(file="m68.rds")
m71<-readRDS(file="m71.rds")
m77<-readRDS(file="m77.rds")
m78<-readRDS(file="m78.rds")

#setup tech and celltype
m64<-RenameCells(m64,add.cell.id="m64",for.merge=T)
m64@meta.data$tech<-"meta"
m64@meta.data$celltype<-"meta_64"

m68<-RenameCells(m68,add.cell.id="m68",for.merge=T)
m68@meta.data$tech<-"meta"
m68@meta.data$celltype<-"meta_68"

m71<-RenameCells(m71,add.cell.id="m71",for.merge=T)
m71@meta.data$tech<-"meta"
m71@meta.data$celltype<-"meta_71"

m77<-RenameCells(m77,add.cell.id="m77",for.merge=T)
m77@meta.data$tech<-"meta"
m77@meta.data$celltype<-"meta_77"

m78<-RenameCells(m78,add.cell.id="m78",for.merge=T)
m78@meta.data$tech<-"meta"
m78@meta.data$celltype<-"meta_78"

#merge
m64_68<-merge(m64,m68)
m71_77<-merge(m71,m77)
m64_68_71_77<-merge(m64_68,m71_77)
m5<-merge(m64_68_71_77,m78)
saveRDS(m5,file="m5_before_integrate.rds")
mpn<-readRDS(file="m5_before_integrate.rds")

#before integrate
mpn[["percent.mt"]] <- PercentageFeatureSet(mpn, pattern = "^Mt-")
VlnPlot(mpn, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = mpn, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
plot_grid(p1,p2)
#before integrate
mpn[["percent.mt"]] <- PercentageFeatureSet(mpn, pattern = "^Mt-")
VlnPlot(mpn, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = mpn, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
plot_grid(p1,p2)

#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(
    pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
    verbose = FALSE)
}
reference.list <- pancreas.list[c("meta_64","meta_68","meta_71","meta_77","meta_78" )]

#setup k.anchor and k.filter correctly
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list)
#看a<-pancreas.anchors@object.list,挑選最小的數字作爲k.weight = 79.
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, k.weight = 79)#留心k.weight,默認100.
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "m5_after_integrated.rds")

hms_individual_integrated<-readRDS(file="m5_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))

#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE)
saveRDS(hms_cluster, file = "m5_cluster_test_1.2.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.001, 
                                logfc.threshold = 0.001
)
write.table(scRNA.markers,file="m5_cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="m5_cell_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="m5_genes_1.2.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="m5_cluster_test_1.2.rds")
DimPlot(hms_cluster, label=TRUE,reduction = "umap")

new.cluster.ids <- c("cytotoxic_cd8", "Recently_activated_cd4","cytotoxic_cd8",
                     "epithelial", "Effector_memory_cd4","macrophage", "Treg",
                     "epithelial", "neutrophil", "B","macrophage","B","macrophage",
                     "macrophage", "progenitor") 

names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "m5_cluster_id.rds")

setwd("~/gse123902/CytoTRACE")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="m5_cluster_id.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4,subsamplesize = 1000)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p

#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
    "0"="cytotoxic_cd8","1"="Recently_activated_cd4","2"="cytotoxic_cd8",
    "3"="epithelial", "4"="Effector_memory_cd4","5"="macrophage","6"="Treg",
    "7"="epithelial", "8"="neutrophil","9"="B","10"="macrophage","11"="B",
    "12"="macrophage","13"="macrophage", "14"="progenitor" ) 
colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)
plotCytoTRACE(result01, phenotype = phe)

#normal samples
a<-"GSM3516666_MSK_LX675_NORMAL_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "n66.rds")

a<-"GSM3516673_MSK_LX682_NORMAL_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "n73.rds")

a<-"GSM3516675_MSK_LX684_NORMAL_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "n75.rds")

a<-"GSM3516676_MSK_LX685_NORMAL_dense.csv"
a1<-data.frame(fread(a),check.names=FALSE, row.names=1)
a2<-t(a1)
pbmc<- CreateSeuratObject(counts = a2, project = "p1", min.cells = 3, min.features = 200)
pbmc
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")#线粒体基因占比计算
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
pbmc <- subset(pbmc, subset = nFeature_RNA > 200  & percent.mt < 5)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc,ndims = 40)
pbmc <- RunUMAP(pbmc, dims = 1:30)
saveRDS(pbmc, file = "n76.rds")

n66<-readRDS(file="n66.rds")
n73<-readRDS(file="n73.rds")
n75<-readRDS(file="n75.rds")
n76<-readRDS(file="n76.rds")

#setup tech and celltype
n66<-RenameCells(n66,add.cell.id="n66",for.merge=T)
n66@meta.data$tech<-"normal"
n66@meta.data$celltype<-"normal_66"

n73<-RenameCells(n73,add.cell.id="n73",for.merge=T)
n73@meta.data$tech<-"normal"
n73@meta.data$celltype<-"normal_73"

n75<-RenameCells(n75,add.cell.id="n75",for.merge=T)
n75@meta.data$tech<-"normal"
n75@meta.data$celltype<-"normal_75"

n76<-RenameCells(n76,add.cell.id="n76",for.merge=T)
n76@meta.data$tech<-"normal"
n76@meta.data$celltype<-"normal_76"

#merge
n66_73<-merge(n66,n73)
n75_76<-merge(n75,n76)
n4<-merge(n66_73,n75_76)

saveRDS(n4,file="n4_before_integrate.rds")
mpn<-readRDS(file="n4_before_integrate.rds")

#before integrate
mpn[["percent.mt"]] <- PercentageFeatureSet(mpn, pattern = "^Mt-")
VlnPlot(mpn, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
pancreas <- NormalizeData(object = mpn, normalization.method = "LogNormalize", scale.factor = 1e4)
pancreas <- FindVariableFeatures(pancreas, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
pancreas <- ScaleData(pancreas, verbose = FALSE)
pancreas <- RunPCA(pancreas, npcs = 30, verbose = FALSE)
pancreas <- RunUMAP(pancreas, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
DimPlot(pancreas, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE)
plot_grid(p1,p2)

#integrate
pancreas.list <- SplitObject(pancreas, split.by = "celltype")
for (i in 1: length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(
    pancreas.list[[i]], selection.method = "vst", nfeatures = 2000, 
    verbose = FALSE)
}
reference.list <- pancreas.list[c("normal_66","normal_73","normal_75","normal_76" )]

#setup k.anchor and k.filter correctly
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list)
#看a<-pancreas.anchors@object.list,挑選最小的數字作爲k.weight = 79.
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, k.weight = 79)#留心k.weight,默認100.
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
saveRDS(pancreas.integrated, file = "n4_after_integrated.rds")

hms_individual_integrated<-readRDS(file="n4_after_integrated.rds")
p1 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(hms_individual_integrated, reduction = "umap", group.by = "celltype")
plot_grid(p1,p2)
hms_neighbor<- FindNeighbors(hms_individual_integrated, dims = 1:30)
obj <- FindClusters(hms_neighbor, resolution = seq(0.5,1.2,by=0.1))

#resolution設置在0.4-1.2之間,越大clusters越多,查看拐點
clustree(obj)
hms_cluster <- FindClusters( hms_neighbor, resolution = 1.2)
head(Idents(hms_cluster), 5)
hms_cluster<- RunUMAP(hms_cluster, dims = 1:30)
DimPlot(hms_cluster, reduction = "umap",label = TRUE)
saveRDS(hms_cluster, file = "n4_cluster_test.rds")

scRNA.markers <- FindAllMarkers(hms_cluster, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.25, 
                                logfc.threshold = 0.25
)
write.table(scRNA.markers,file="n4_cellMarkers.txt",sep="\t",row.names=F,quote=F)

#挑选每个细胞亚群中特意高表达的20个基因
top20 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
write.csv(file="n4_markers.csv",top20)
#整理成表格，只显示基因名字
top20_table=unstack(top20, gene ~ cluster)
names(top20_table)=gsub("X","cluster",names(top20_table))
write.csv(file="n4_genes.csv",top20_table,row.names=F)

#细胞及细胞中基因与RNA数量
hms_cluster<-readRDS(file="n4_cluster_test.rds")
DimPlot(hms_cluster, reduction = "umap",label = TRUE)

slotNames(hms_cluster)
#assay
hms_cluster@assays
dim(hms_cluster@meta.data)
View(hms_cluster@meta.data)

hms_cluster<-readRDS(file="n4_cluster_test.rds")
DimPlot(hms_cluster, label=TRUE,reduction = "umap")
new.cluster.ids <- c("Thelper", "Cytotoxic_cd8","Naive_cd4","Cytotoxic_cd8","Macrophage","Macrophage",
 "Cytotoxic_cd8","Progenitor_cd8 ","Monocyte","Monocyte","Macrophage","Macrophage","Cytotoxic_cd8",
 "Cytotoxic_cd8","Macrophage","Epithelial", "Effector_memory_cd8", "Epithelial","Mast","Progenitor",
 "Epithelial","B","Chronic_activation_cd4","Macrophage","Epithelial","B","Macrophage")
names(new.cluster.ids) <- levels(hms_cluster)
hms_cluster_id<- RenameIdents(hms_cluster, new.cluster.ids)
DimPlot(hms_cluster_id, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(hms_cluster_id, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(hms_cluster_id, file = "n4_cluster_id.rds")

hms_cluster_id<-readRDS("n4_cluster_id.rds")

setwd("~/gse123902/CytoTRACE")
#devtools::install_local('CytoTRACE_0.3.3.tar.gz')
library(CytoTRACE)
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
#library(monocle)#can not use this package.error will come out.
seurat <- readRDS(file="n4_cluster_id.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
a<-as.matrix(data)
result01<-CytoTRACE(a,ncores=4,subsamplesize = 1000)#ncores only can be 4 , can not be 8.
plotCytoGenes(result01, numOfGenes = 10)#first figure
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p1
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p

#assigned_cell_type
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
"0"="Thelper", "1"="Cytotoxic_cd8","2"="Naive_cd4","3"="Cytotoxic_cd8","4"="Macrophage","5"="Macrophage",
"6"="Cytotoxic_cd8","7"="Progenitor_cd8 ","8"="Monocyte","9"="Monocyte","10"="Macrophage","11"="Macrophage",
"12"="Cytotoxic_cd8","13"="Cytotoxic_cd8","14"="Macrophage","15"="Epithelial","16"="Effector_memory_cd8",
"17"="Epithelial","18"="Mast","19"="Progenitor","20"="Epithelial","21"="B","22"="Chronic_activation_cd4",
"23"="Macrophage","24"="Epithelial","25"="B","26"="Macrophage")
colnames(colData(cds))
head(colData(cds))
a<-colData(cds)
write.csv(a,"a")
table(colData(cds)$assigned_cell_type)
phe<-colData(cds)$assigned_cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat@meta.data)
plotCytoTRACE(result01, phenotype = phe)

#Dotplot
setwd("~/gse123902/CytoTRACE")
hms_cluster<-readRDS("n4_cluster_id.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("CD79A","IL32","TRAC","NKG7","PRF1","FGFBP2","GZMH", "IL7R","TPM2",
             "KRT19", "CAV1", "FTL", "GATA2", "LST1", "LDHB", "CD8A", "RGCC","
           CD3D","CCL5","CST7", "KLDR1","GZMA","GZMB", "CD247","CTSW","HLA-DRB1"
             ,"GPX1", "VIM","CD74","S100A6", "CST3","HLA-DRA","S10A11","FTH1")
a<-DotPlot(hms_cluster,features=features)+RotatedAxis() 
b<-a$data
write.csv(b,"b")
FeaturePlot(hms_cluster,features=c("NKG7","FTH1","PRF1", "FTL","FGFBP2",
                                  "HLA-DRA"), cols=c("white","blue"),
            min.cutoff=2,ncol = 2)

setwd("~/gse123902/CytoTRACE")
hms_cluster<-readRDS("p8_cluster_id.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("CD79A","CCL5", "TPM2",  "KRT19", "CAV1", "FTL",  "CD8A",  "ANXA1", 
           "KLRD1","PGAP1", "AIF1", "SCGB1A1", "CXCL13", "GNLY", "IFI6",
            "ISG15", "CD69",  "DUSP4","TMSB10","ACTB","SERF2","S100A6","MYL6","S100A11",
           "CD63","ARPC3","GAPDH","TUBA4A","KLRB1","FYN","ZFP36L2","CNOT6L",
           "SYTL3","BTG1","IL7R","CXCR4")
a<-DotPlot(hms_cluster,features=features)+RotatedAxis() 
b<-a$data
write.csv(b,"b")
features = c("IL7R", "BTG1", "CXCR4","CCL5","FTL","TMSB10")
FeaturePlot(hms_cluster,features=features,ncol=2, cols=c("white","blue"),
            min.cutoff=1.5)

setwd("~/gse123902/CytoTRACE")
hms_cluster<-readRDS("m5_cluster_id.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("COL1A2","CD74","IL32",
             "NKG7","IL7R","KRT19",
   "GNLY","CCL5","S100A9","LYZ","IGHA1","HLA-DRA",
   "RGS1","PLCG2","CXCR4","JUNB","ZFP36L2" ,"FTH1","FTL","ATP5E","RPL8","MYL6",
   "H3F3A","RPL37A","SERF2","RPL39","PTMA","DERL3","CD79A","AL928768.3","TNFRSF17","IGKV1-16",
   "MZB1","IGHA2","IGKC","JCHAIN")
a<-DotPlot(hms_cluster,features=features)+RotatedAxis() 
b<-a$data
write.csv(b,"b")
features = c("FTH1", "FTL","CXCR4","CCL5")
FeaturePlot(hms_cluster,features=features,ncol=2, cols=c("white","blue"),
            min.cutoff=2)



















#DimPlotX
DimPlot(B_cell, reduction = "umap", split.by = "tech")
DimPlot(Natural_Killer_T, reduction = "umap", split.by = "tech")
DimPlot(Dendritic_cell, reduction = "umap", split.by = "tech")
DimPlot(CD8_T, reduction = "umap", split.by = "tech")
DimPlot(Progenitor, reduction = "umap", split.by = "tech")
DimPlot(Mast, reduction = "umap", split.by = "tech")
DimPlot(Exhausted_CD8, reduction = "umap", split.by = "tech")
DimPlot(Microglial, reduction = "umap", split.by = "tech")
DimPlot(Exhausted_CD48, reduction = "umap", split.by = "tech")
DimPlot(Naive_CD8, reduction = "umap", split.by = "tech")

RidgePlot(CD8_T, features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"),cols = c("green3","orangered","red"), group.by="tech", ncol = 4) + theme(axis.title.y = element_blank())

markers.to.plot<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2","TMSB4X")
DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
DoHeatmap(subset(Exhausted_CD8,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
DoHeatmap(subset(Exhausted_CD48,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,label = FALSE)
DoHeatmap(subset(Naive_CD8,downsample=50000),features = markers.to.plot,size=5,label = FALSE)


genes11_heatmap<-DotPlot(hms_cluster_id,features = c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2"))+RotatedAxis()
genes11_heatmap<-genes11_heatmap$data
write.csv(genes11_heatmap,"genes11_heatmap")

#expression level in each patients
setwd("~/geo/gse123902")
CD8_T<-readRDS("CD8_T.rds")
a<-DoHeatmap(subset(CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_cd8T")

hms_cluster_id<-readRDS("hms_cluster_id_test.rds")
a<-DoHeatmap(subset(hms_cluster_id,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_hms")

ex_CD8_T<-readRDS("Exhausted_CD8.rds")
a<-DoHeatmap(subset(ex_CD8_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_ex_cd8T")

ex_CD48_T<-readRDS("Exhausted_CD48.rds")
a<-DoHeatmap(subset(ex_CD48_T,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_ex_cd48T")

Naive_CD8<-readRDS("Naive_CD8.rds")
a<-DoHeatmap(subset(Naive_CD8,downsample=50000),features = markers.to.plot,size=5,group.by="celltype",label = FALSE)
a1<-a$data
write.table(a1,"a1_Naive_CD8")