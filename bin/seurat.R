require(scales)
library(Seurat)
require(sctransform)
library(glmGamPoi)

library(GenomeInfoDb)

library(ggplot2)
require(reshape)
library(patchwork)
require(data.table)
require(stringr)
require(GenomicRanges)

set.seed(1234)


## For PNG file
options(bitmapType="cairo")

set.seed(1234)

## Reading input from cellranger
counts = Read10X("outs/raw_feature_bc_matrix")
heart <- CreateSeuratObject(counts = counts, project = "EB", min.cells = 10, min.features = 300)
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "MT-")
##saveRDS(heart,file="heart.preQC.RDS")

##options(bitmapType='cairo')
##png("violinplot_QC.png",width=4000,height=2000,res=300)
##VlnPlot(heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
##dev.off()


##png("featureScatter_QC.png",width=4000,height=2000,res=300)
##plot1 <- FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "percent.mt")
##plot2 <- FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
##plot1 + plot2
##dev.off()

## Filtering according uniform parameters compared across multiple runs
##heart = readRDS("heart.preQC.RDS")
heart <- subset(heart, subset = nFeature_RNA > 3000 & nFeature_RNA < 15000 & nCount_RNA>5000 & nCount_RNA<120000 & percent.mt < 25)
##saveRDS(heart, file="heart.postQC.RDS")

## Normalize data using SCTransform
heart <- SCTransform(heart, method="glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
heart <- RunPCA(heart)

## diagnostic
pc = as.data.frame(heart@reductions$pca@cell.embeddings)

meta.data = heart@meta.data
meta.data$nCount_RNA = log10(meta.data$nCount_RNA)
meta.data$nFeature_RNA = log10(meta.data$nFeature_RNA)
meta.data$nCount_SCT = log10(meta.data$nCount_SCT)
meta.data$nFeature_SCT = log10(meta.data$nFeature_SCT)
## remove orig.ident
meta.data = meta.data[,-1]


cor.pc.feature = cor(pc,meta.data)
subset.cor.pc.feature = head(cor.pc.feature,10)
cordata = melt(as.matrix(subset.cor.pc.feature))
names(cordata) = c("PC","Feature","r")

cordata = cordata[ cordata$PC!="PC_10", ]
png("PC_correlation_heatmap.png",width=2500,height=2300,res=300)
p1 <- ggplot(cordata,aes(x=PC,y=Feature)) + geom_tile(aes(fill=r)) + theme_bw() + theme(panel.grid=element_blank()) + geom_text(aes(x=PC,y=Feature,label=format(round(r,2),2))) + scale_fill_gradientn(values=rescale(c(-0.5,0.5)),colors=c("red","white","blue"),limits=c(-1,1))
print(p1)
dev.off()

png("top_PC.png",width=2000,height=2000,res=300)
VizDimLoadings(heart, dims = 1:2, reduction = "pca")
dev.off()

png("PCA.png",width=5000,height=5000,res=300)
DimPlot(heart, reduction = "pca")
dev.off()

png("topgene_pca_heatmap.png",width=5000,height=5000,res=300)
DimHeatmap(heart, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

png("elbow.png",width=2000,height=2000,res=300)
ElbowPlot(heart,ndims=50)
dev.off()

##saveRDS(heart,file="heart.rds")

## Find clusters, and UMAP of cuz
heart <- FindNeighbors(heart, dims=1:30)
heart <- FindClusters(heart, resolution = 0.5)
heart <- RunUMAP(heart,dims=1:30)

## FINAL R object for integration
saveRDS(heart, file = "heart.rds")
umap.coord = as.data.frame(heart@reductions$umap@cell.embeddings)
umap.coord$group = as.factor(Idents(heart))
png("UMAP_SCTransform.png",width=3000,height=3000,res=300)
p1 <- ggplot(umap.coord,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=group),size=1.2)+theme_bw()+theme(panel.grid=element_blank())
##p1 <- p1 + scale_color_manual(values=c("darkred","red","orange","yellow","turquoise","blue","midnightblue","black","brown","purple","pink","salmon","bisque","grey30","grey60","green","darkgreen","yellowgreen"))
print(p1)
dev.off()

