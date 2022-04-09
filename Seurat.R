
library(Seurat)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
###########seurat standard pipeline
S2.counts <- Read10X(data.dir = "/Users/weihonggu/Documents/single cell sequencing/Thymus/Thymus_dim16/Thymus_filtered_feature_bc_matrix")
Thy<- CreateSeuratObject(counts = S2.counts, project = "thymus", min.cells = 3, min.features = 200)
dim(Thy_pig)

####QC
Thy[["percent.mt"]] <- PercentageFeatureSet(Thy, features=c("ND1","ND2", "ND3", "ND4",
                                                            "ND4L", "ND5","ND6","COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8"))
VlnPlot(Thy_c, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
VlnPlot(Thy_c, features = c("nFeature_RNA", "percent.mt"), ncol = 2)
plot1 <- FeatureScatter(Thy, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Thy_pig, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot1+plot2
FeatureScatter(Thy_c, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(Thy_cl, feature1 = "nFeature_RNA", feature2 = "percent.mt")
Thy <- subset(Thy, subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & percent.mt < 11)
dim(Thy)
Thy <- NormalizeData(Thy)
Thy <- FindVariableFeatures(Thy, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Thy)
Thy <- ScaleData(Thy, features = all.genes)

####regress out cellcyle 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Thy <- CellCycleScoring(Thy, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Thy<-ScaleData(Thy, vars.to.regress = c("S.Score","G2M.Score"), features = rownames(Thy))
Thy <- RunPCA(Thy, npcs = 30, verbose = FALSE)
ElbowPlot(Thy)

####dimensional reduction, clustering

Thy <- FindNeighbors(Thy, dims = 1:19)
Thy <- FindClusters(Thy, resolution = 0.5)
Thy <- RunUMAP(Thy, dims = 1:19)

DimPlot(Thy, reduction = "umap", label = T)

c_cols <- c('0'='#F8766D','1'='#E68613','2'='#CD9600','3'='#ABA300','4'='#7CAE00',
            '5'='#0CB702','6'='#D7EF16','7'='#00C19A','8'='#910C00','9'='#ED68ED',
            '10'='#00A9FF','11'='#8494FF','12'='#FFBDFF','13'='#00B8E7', '14'='#FF61CC','15'='#FF68A1')

#######gene expression
FeaturePlot(Thy, features = c("CD3D"), label =F)
genes<-c("PTPRC", "CD3D","CD3E", "CD4", "CD40LG", "CD8A", "CD8B", "RAG1","RAG2","PCNA","PTCRA", "IL7R","RORC","ENSSSCG00000029596","SOX13", "CCR9", "CCR7","CTLA4","TNFRSF9","FOXP3","ZBTB7B","RUNX3","SELL", "KLF2", "ISG15","ZBTB16","NKG7","KLRB1","KLRK1", "EOMES","CD19", "PAX5", "CD40")

DotPlot(Thy, features = genes) + RotatedAxis()+coord_flip()+theme(text = element_text(face="bold"))

markers <- FindAllMarkers(Thy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, "Thymus.markers.csv")
markers_5<-markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
write.csv(Thy.cycle.markers_20,'Table S2.csv' )
DoHeatmap(Thy, features = markers_5$gene, group.colors =  c_cols, size = 3)+ theme(text = element_text(size = 8, face="bold"))+NoLegend()

#############mature T cells
ma_cols <- c('5'='#0CB702','6'='#D7EF16','8'='#910C00','9'='#ED68ED',
             '10'='#00A9FF','11'='#8494FF','12'='#FFBDFF','13'='#00B8E7', '14'='#FF61CC')
Thy_ma<-subset(Thy_c, idents = c(5,6,8:14))
DimPlot(Thy_ma, cols =ma_cols , label.size = 6)
cluster.averages_ma <- AverageExpression(Thy_ma, return.seurat=TRUE) 

FeatureScatter(cluster.averages_ma, feature1 = "KLF2", feature2 = "S1PR1", 
               pt.size = 3, cols = ma_cols)+geom_abline()+expand_limits(x = c(0), y = c(0))+geom_text_repel(aes(label = levels(Thy_ma)),
                                                                                                            hjust = 0, nudge_x = 0, size = 6) + NoLegend() # line as y=x

#######subset

thy_unc<-subset(Thy_c, idents = c(5,6,8,9,10,11,12,13,14))
levels(thy_unc)<-c("8", "9","13", "14", "5","6")
DimPlot(thy_unc, reduction = "umap", label = TRUE, cols = c('#0CB702','#D7EF16','#910C00','#ED68ED','#FFBDFF','#00B8E7','#FF61CC'))
DotPlot(thy_unc, features = genes,dot.scale = 5) + RotatedAxis() +coord_flip()+theme(axis.text.y = element_text(size = 9,face="bold"))

#####volano plot
treg_marker<-FindMarkers(thy_unc, ident.1 = '8', ident.2 = 9, min.pct = 0.25 )
EnhancedVolcano(treg_marker,
                lab = rownames(treg_marker),
                x = 'avg_logFC',
                y = 'p_val',
                gridlines.minor=F,
                gridlines.major=F,
                title = 'Trge1 vs Treg2',
                legendPosition = 'top',
                FCcutoff =1.4,
                legendLabSize = 10,
                labSize = 3,
                axisLabSize = 9,
                borderWidth = 0.7,
                drawConnectors = T,
                labFace = 'bold',
                border='full',)

#############integration######

Thy1<-CreateSeuratObject(Thymus1counts_34negative_geneconverted, min.cells = 3, min.features = 200, project = "human")
Thy1[["percent.mt"]] <- PercentageFeatureSet(Thy1, features=c("ND1","ND2", "ND3", "ND4",
                                                              "ND4L", "ND5","ND6","COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8"))
VlnPlot(Thy1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(Thy1, feature1 = "nFeature_RNA", feature2  = "percent.mt")
Thy1<-subset(Thy1, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)
Thy1 <- NormalizeData(Thy1)
Thy1<-FindVariableFeatures(Thy1, selection.method = "vst", nfeatures = 2000)
dim(Thy1)

####pig
S2.counts <- Read10X(data.dir = "/Users/weihonggu/Documents/single cell sequencing/Thymus/Thymus_dim16/Thymus_filtered_feature_bc_matrix")
Thy_pig <- CreateSeuratObject(counts = S2.counts, project = "pig", min.cells = 3, min.features = 200)
Thy_pig[["percent.mt"]] <- PercentageFeatureSet(Thy_pig, features=c("ND1","ND2", "ND3", "ND4",
                                                                    "ND4L", "ND5","ND6","COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8"))
VlnPlot(Thy_pig, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(Thy_pig, features = c("percent.mt"))

FeatureScatter(Thyp, feature2 = "nFeature_RNA", feature1 = "percent.mt")
Thy_pig <- subset(Thy_pig, subset = nFeature_RNA > 800 & nFeature_RNA < 5000 & percent.mt < 11)
dim(Thy_pig)

Thy_pig <- NormalizeData(Thy_pig)
Thy_pig <- FindVariableFeatures(Thy_pig, selection.method = "vst", nfeatures = 2000)
Thy.list <- list(Thy1, Thy_pig)

Thy.anchors <- FindIntegrationAnchors(object.list = Thy.list, dims = 1:30)
Thy.integrated <- IntegrateData(anchorset = Thy.anchors, dims = 1:30)
DefaultAssay(Thy.integrated) <- "integrated"

Thy.integrated <- ScaleData(Thy.integrated)
thy_hp<-CellCycleScoring(Thy.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = T )
thy_hp<-ScaleData(thy_hp, vars.to.regress = c("S.Score","G2M.Score"), features = rownames(thy_hp))
thy_hp <- RunPCA(thy_hp, features = VariableFeatures(thy_hp))
ElbowPlot(thy_hp)
thy_hp <- FindNeighbors(thy_hp, dims = 1:19)
thy_hp <- FindClusters(thy_hp, resolution = 0.5)
thy_hp <- RunUMAP(thy_hp, dims = 1:19)

DimPlot(thy_hp, reduction = "umap", label = T)
DimPlot(thy_hp, reduction = "umap", group.by = "orig.ident")
DimPlot(thy_hp, reduction = "umap", split.by = "orig.ident",label = T)

DefaultAssay(thy_hp) <- "RNA"
FeaturePlot(thy_hp, features = c("CD8A"), split.by = "orig.ident", 
            cols = c("grey", "red"), label = F) #CD1A

thy_hp.markers<- FindAllMarkers(thy_hp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

CD8.markers <- FindConservedMarkers(thy_hp, ident.1 = 13, grouping.var = "orig.ident", verbose = FALSE)
NK8.markers <- FindConservedMarkers(thy_hp, ident.1 = 14, grouping.var = "orig.ident", verbose = FALSE)
write.csv(CD8.markers, "CD8aa_conserved_pig_human.csv")
write.csv(NK8.markers, "NK_like_conserved_pig_human.csv")

####Add individual clusters between species 
thy_hp$celltype<-paste(Idents(thy_hp), thy_hp$orig.ident, sep="_")
Idents(thy_hp)<-"celltype"

####DEGs
hp13<-FindMarkers(thy_hp, ident.1 = "13_human", ident.2 = "13_pig", min.pct = 0.25, logfc.threshold = 0.25)
