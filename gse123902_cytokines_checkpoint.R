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

setwd("~/gse123902/CytoTRACE")
hms_cluster<-readRDS("n4_cluster_id.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("CXCL15","CXCL14","CXCL12","CXCL10","CXCL9","CXCL5","CXCL3","CXCL2","CXCL1",
             "CCL1","CCL5","CCL4","CCL3","CCL2","CCL7",
             "IL17A", "IL12A", "IL10", "IL6", "IL4", "IL2","IL1B",
             "TNF","IFNG","TGFB1","CSF2","CSF1","CD40","NOS2","ARG1","S100A8","S100A9")
DotPlot(hms_cluster,features=features)+RotatedAxis()  

features = c("CXCR4","CXCR3","CXCR2","CCR1", "CCR5", "CCR4", "CCR2", 
             "IL17RC", "IL17RA", "IL12RB2","IL12RB1", "IL10RB","IL10RA",
             "IL6RA","IL4RA","IL2RB","IL2RG","IL2RA",
             "TNFRSF1B", "TNFRSF1A","IFNGR2","IFNGR1",
             "TGFBR3","TGFBR2" ,"TGFBR1","CSF3R","CSF2RA","CSF1R","IL1R1")
DotPlot(hms_cluster,features=features)+RotatedAxis()  

features = c("CD27","CD28","CD40","TNFRSF4","TNFRSF9","TNFRSF18","ICOS","CD266",
             "ADORA2A","BTLA","CTLA4","LAG3","PDCD1","HAVCR2","VSIR","TIGIT","CD70",
             "CD80","CD86","CD4oLG","TNFSF9","TNFSF4","ICOSL","NECTIN2",
             "PVR","LGALS3","CD274","PDCD1IG2","CEACAM1","LGALS9")
DotPlot(hms_cluster,features=features)+RotatedAxis()  

setwd("~/gse123902/CytoTRACE")
hms_cluster<-readRDS("p8_cluster_id.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("CXCL15","CXCL14","CXCL12","CXCL10","CXCL9","CXCL5","CXCL3","CXCL2","CXCL1",
             "CCL1","CCL5","CCL4","CCL3","CCL2","CCL7",
             "IL17A", "IL12A", "IL10", "IL6", "IL4", "IL2","IL1B",
             "TNF","IFNG","TGFB1","CSF2","CSF1","CD40","NOS2","ARG1","S100A8","S100A9")
DotPlot(hms_cluster,features=features)+RotatedAxis()  

features = c("CXCR4","CXCR3","CXCR2","CCR1", "CCR5", "CCR4", "CCR2", 
             "IL17RC", "IL17RA", "IL12RB2","IL12RB1", "IL10RB","IL10RA",
             "IL6RA","IL4RA","IL2RB","IL2RG","IL2RA",
             "TNFRSF1B", "TNFRSF1A","IFNGR2","IFNGR1",
             "TGFBR3","TGFBR2" ,"TGFBR1","CSF3R","CSF2RA","CSF1R","IL1R1")
DotPlot(hms_cluster,features=features)+RotatedAxis()  

features = c("CD27","CD28","CD40","TNFRSF4","TNFRSF9","TNFRSF18","ICOS","CD266",
             "ADORA2A","BTLA","CTLA4","LAG3","PDCD1","HAVCR2","VSIR","TIGIT","CD70",
             "CD80","CD86","CD4oLG","TNFSF9","TNFSF4","ICOSL","NECTIN2",
             "PVR","LGALS3","CD274","PDCD1IG2","CEACAM1","LGALS9")
DotPlot(hms_cluster,features=features)+RotatedAxis()  

setwd("~/gse123902/CytoTRACE")
hms_cluster<-readRDS("m5_cluster_id.rds")
DimPlot(hms_cluster, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(hms_cluster, reduction = "umap", label = FALSE, pt.size = 0.5) 
features = c("CXCL15","CXCL14","CXCL12","CXCL10","CXCL9","CXCL5","CXCL3","CXCL2","CXCL1",
             "CCL1","CCL5","CCL4","CCL3","CCL2","CCL7",
             "IL17A", "IL12A", "IL10", "IL6", "IL4", "IL2","IL1B",
             "TNF","IFNG","TGFB1","CSF2","CSF1","CD40","NOS2","ARG1","S100A8","S100A9")
DotPlot(hms_cluster,features=features)+RotatedAxis()  

features = c("CXCR4","CXCR3","CXCR2","CCR1", "CCR5", "CCR4", "CCR2", 
             "IL17RC", "IL17RA", "IL12RB2","IL12RB1", "IL10RB","IL10RA",
             "IL6RA","IL4RA","IL2RB","IL2RG","IL2RA",
             "TNFRSF1B", "TNFRSF1A","IFNGR2","IFNGR1",
             "TGFBR3","TGFBR2" ,"TGFBR1","CSF3R","CSF2RA","CSF1R","IL1R1")
DotPlot(hms_cluster,features=features)+RotatedAxis()  

features = c("CD27","CD28","CD40","TNFRSF4","TNFRSF9","TNFRSF18","ICOS","CD266",
             "ADORA2A","BTLA","CTLA4","LAG3","PDCD1","HAVCR2","VSIR","TIGIT","CD70",
             "CD80","CD86","CD4oLG","TNFSF9","TNFSF4","ICOSL","NECTIN2",
             "PVR","LGALS3","CD274","PDCD1IG2","CEACAM1","LGALS9")
DotPlot(hms_cluster,features=features)+RotatedAxis()  


