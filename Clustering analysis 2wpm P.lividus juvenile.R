library(Seurat)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(dplyr)
library(viridis)
library(ape)
library(matrixStats)
library(ggtree)
library(treeio)
library(tidytree)
library(bioDist)
library(svglite)
set.seed(255)

#Create Seurat objects

#1
pljuv_1.data <- Read10X(data.dir = "./pljuv_1")
pljuv_1_barcodes <- read.delim(file = "./pljuv_1/barcodes.tsv", stringsAsFactors = F, header = F)
pljuv_1 <- CreateSeuratObject(pljuv_1.data, project = "pljuv_1", min.cells = 3, min.features = , meta.data = rownames(pljuv_1_barcodes$v1)) #check values for cells and features

VlnPlot(pljuv_1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
pljuv_1- subset(pljuv_1, subset = nFeature_RNA > 350 & nFeature_RNA < 5000)

#2
pljuv_2.data <- Read10X(data.dir = "./pljuv_2")
pljuv_2_barcodes <- read.delim(file = "./pljuv_2/barcodes.tsv", stringsAsFactors = F, header = F)
pljuv_2 <- CreateSeuratObject(pljuv_2.data, project = "pljuv_2", min.cells = 3, min.features = , meta.data = rownames(pljuv_2_barcodes$v1)) #check values for cells and features

VlnPlot(pljuv_2, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
pljuv_2- subset(pljuv_2, subset = nFeature_RNA > 500 & nFeature_RNA < 5000)

#3
pljuv_3.data <- Read10X(data.dir = "./pljuv_3") #They call genes.tsv file barcode file for some reason, if its not present
pljuv_3_barcodes <- read.delim(file = "./pljuv_3/barcodes.tsv", stringsAsFactors = F, header = F)
pljuv_3 <- CreateSeuratObject(pljuv_3.data, project = "pljuv_3", min.cells = 3, min.features = , meta.data = rownames(pljuv_3_barcodes$v1)) #check values for cells and features

VlnPlot(pljuv_3, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
pljuv_3- subset(pljuv_3, subset = nFeature_RNA > 500 & nFeature_RNA < 5000)

pljuv_1 <- RenameCells(pljuv_1, add.cell.id = "pljuv_1")
pljuv_2 <- RenameCells(pljuv_2, add.cell.id = "pljuv_2")
pljuv_3 <- RenameCells(pljuv_3, add.cell.id = "pljuv_3")

#Normalization, anchoring and data integration
pljuv_list <- list(pljuv_1, pljuv_2, pljuv_3)

for (i in 1:length(pljuv_list)) {
  
  pljuv_list[[i]] <- NormalizeData(pljuv_list[[i]], 
                                   verbose = FALSE)
  pljuv_list[[i]] <- FindVariableFeatures(pljuv_list[[i]], 
                                          selection.method = "vst",
                                          nfeatures = 2000,
                                          verbose = FALSE)
  
}
pljuv_anchors <- FindIntegrationAnchors(object.list = pljuv_list,
                                        dims = 1:50)

pljuv <- IntegrateData(anchorset = pljuv_anchors, dims = 1:50)

#Clustering                                                                                         
pljuv_integrated <- ScaleData(object = pljuv)
pljuv_integrated <- RunPCA(object = pljuv_integrated, features = VariableFeatures(object = pljuv_integrated))
pljuv_integrated <- FindNeighbors(object = pljuv_integrated, dims = 1:50)
pljuv_integrated<- FindClusters(object = pljuv_integrated, resolution = 1)

#UMAP
pljuv_integrated_umap <- RunUMAP(object = pljuv_integrated, dims = 1:50)
DimPlot(pljuv_integrated_umap, group.by="orig.ident")
DimPlot(object = pljuv_integrated_umap, reduction = "umap", label = TRUE)
DimPlot(pljuv_integrated_umap, group.by="orig.ident", split.by = "orig.ident")

#renaming and reorganizing
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `41` = "48 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `12` = "47 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `5` = "46 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `15` = "45 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `6` = "44 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `17` = "43 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `33` = "42 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `4` = "41 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `7` = "40 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `36` = "39 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `25` = "38 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `30` = "37 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `1` = "36 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `2` = "35 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `16` = "34 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `11` = "33 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `3` = "32 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `8` = "31 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `9` = "30 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `0` = "29 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `10` = "28 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `13` = "27 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `14` = "26 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `18` = "25 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `19` = "24 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `20` = "23 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `21` = "22 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `22` = "21 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `23` = "20 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `24` = "19 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `26` = "18 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `27` = "17 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `28` = "16 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `29` = "15 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `35` = "14 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `37` = "13 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `39` = "12 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `42` = "11 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `43` = "10 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `45` = "9 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `47` = "8 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `44` = "7 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `40` = "6 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `38` = "5 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `34` = "4 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `31` = "3 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `32` = "2 ")
pljuv_integrated_umap<- RenameIdents(object = pljuv_integrated_umap, `46` = "1 ")

#marker identidication
markers <- FindAllMarkers(object = pljuv_integrated_umap, assay = "RNA", only.pos = TRUE, min.pct = 0.01 )

#labelling of cell type groups 
neurons <- WhichCells(pljuv_integrated_umap, idents = c("1 ", "2 ", "3 ", "4 ", "5 ", "6 ", "7 ", "8 ", "9 ", "10 ", "11 ", "12 ", "13 ", "14 ", "15 ", "16 ", "17 ", "18 ", "19 ", "20 ", "21 ", "22 ", "23 ", "24 ", "25 ", "26 ", "27 ", "28 ", "29 "))                    
podia_epidermins <- WhichCells(pljuv_integrated_umap, idents = c("30 ","31 ", "32 ", "33 ", "33 "))
epidermis<- WhichCells(pljuv_integrated_umap, idents = c("34 ", "35 ", "36 ", "37 ", "38"))
digestive_tract<- WhichCells(pljuv_integrated_umap, idents = c("39 ", "40 ", "41 "))
wvs<- WhichCells(pljuv_integrated_umap, idents = c("42 ", "43 ", "44 "))
skeleton<- WhichCells(pljuv_integrated_umap, idents = "45 ")
muscles<-WhichCells(pljuv_integrated_umap, idents = "46 ")
coelomocytes<-WhichCells(pljuv_integrated_umap, idents = c("47 ","48 "))
DimPlot(pljuv_integrated_umap, label=TRUE, , cells.highlight= list(neurons, podia_epidermins, epidermis, digestive_tract, wvs, skeleton, muscles, coelomocytes), cols.highlight = c("orange", "salmon","pink","#b3a0006d","blue","#3da1ff78","#00bf7661","magenta"), cols= "#be80ff6d")


#ClusterTree

pljuv_integrated_umap[["major_clusters"]] <- Idents(object = pljuv_integrated_umap)

# Neighbor-joining Tree Analysis

## prep expression matrices

pl_juv_allcells_avgmatrix <- AverageExpression(object = pljuv_integrated_umap, return.seurat = F)


Identify genes expressed in at least one cluster or subcluster above 3 transcripts per million:

expressed_genes <- rownames(pl_juv_allcells_avgmatrix$RNA)[which(rowMaxs(as.matrix(pl_juv_allcells_avgmatrix$RNA)) > 0.03)]


clusters <- levels(pljuv_integrated_umap)
clusters<-paste0("g", clusters)

#major clusters expressed genes matrix
cluster_matrix_allgenes_expressed <- as.matrix(pl_juv_allcells_avgmatrix$RNA[expressed_genes,clusters])


## NJ tree analyses of major clusters using all expressed genes

#all genes expressed
nj_majorclusters_spearman_TPM_allgenes <- nj(spearman.dist(t(cluster_matrix_allgenes_expressed)))
nj_majorclusters_spearman_TPM_allgenes_boot <- boot.phylo(phy = nj_majorclusters_spearman_TPM_allgenes, 
                                                          t(cluster_matrix_allgenes_expressed),
                                                          FUN = function(xx) nj(spearman.dist(xx)),
                                                          B = 10000)
nj_majorclusters_spearman_TPM_allgenes$node.label <- nj_majorclusters_spearman_TPM_allgenes_boot/100
nj_majorclusters_spearman_TPM_allgenes$node.label <- round(nj_majorclusters_spearman_TPM_allgenes$node.label)


# Tree Visualization

pl_juv_meta_tibble <- tibble(cellID = colnames(pljuv_integrated_umap), 
                             nCount = pljuv_integrated_umap@meta.data$nCount_RNA, 
                             nFeature = pljuv_integrated_umap@meta.data$nFeature_RNA,
                             major_cluster = pljuv_integrated_umap@meta.data$major_clusters)

pl_juv_meta_tibble_major_clusters <- pl_juv_meta_tibble %>%
  filter(major_cluster %in% clusters) %>%
  group_by(major_cluster) %>%
  summarize(cell_n = n(), nCount_avg = mean(nCount), nFeature_avg = mean(nFeature))


# generate ggtree object
gt <- ggtree(nj_majorclusters_spearman_TPM_allgenes, layout = "rectangular")



# generate tree with black tips
gt3 <- gt %<+% pl_juv_meta_tibble_major_clusters + 
  geom_text2(aes(subset = !isTip, label=label), nudge_x = -.01, nudge_y = .2, size = 2) +
  geom_tiplab(size = 4, colour = "black", linesize = 2) +
  ggplot2::xlim(0, 0.5)


# generate tree with metadata annotations
gt4 <- gt %<+% pl_juv_meta_tibble_major_clusters +
  geom_tiplab(size = 4, linesize = 2, offset = .02) +
  geom_tippoint(aes(size = nFeature_avg, colour = nCount_avg)) +
  #scale_size_continuous(range = c(1,4)) +
  scale_colour_gradientn(colours = plasma(n=7, direction = 1)) +
  ggplot2::xlim(0, 0.5)


#subset PRCs
nervous_system<-subset(pljuv_integrated_umap, idents=c("1 ", "2 ", "3 ", "4 ", "5 ", "6 ", "7 ", "8 ", "9 ", "10 ", "11 ", "12 ", "13 ", "14 ", "15 ", "16 ", "17 ", "18 ", "19 ", "20 ", "21 ", "22 ", "23 ", "24 ", "25 ", "26 ", "27 ", "28 ", "29 "))
photoreceptors<-subset(nervous_system, `Opsin3.2-3--8080`>1| `Opsin2--17910`>1 |`Opsin3.2-2--25765`>1|`Opsin3.2--25743`>1|`UnkProt-4571--18324`>1|`Opn5L--17817`>1|`Opn5L-2--17895`>1)
DefaultAssay(photoreceptors) <- "RNA"
photoreceptors <- FindVariableFeatures(photoreceptors, selection.method = "vst", nfeatures = 2000)
DefaultAssay(photoreceptors) <- "integrated"
photoreceptors <- ScaleData(photoreceptors, verbose = FALSE, vars.to.regress = c( "nCount_RNA"))
photoreceptors <- RunPCA(photoreceptors, npcs = 50, verbose = FALSE)
photoreceptors_js<-JackStraw(photoreceptors, reduction = "pca", assay = NULL, dims = 50,
                             num.replicate = 100, prop.freq = 0.01, verbose = TRUE,
                             maxit = 1000)
ScoreJackStraw(photoreceptors_js, reduction = "pca", dims = 1:50,
               score.thresh = 1e-05, do.plot = TRUE)
photoreceptors <- FindNeighbors(photoreceptors, dims = 1:19)
photoreceptors <- FindClusters(object = photoreceptors, resolution = 1)
photoreceptors <- RunUMAP(photoreceptors, dims = 1:19)
DefaultAssay(photoreceptors) <- "RNA"
DimPlot(object = photoreceptors, reduction = "umap", label = FALSE)

photoreceptors<- RenameIdents(object = photoreceptors, `14` = "PRCs (15)")
photoreceptors<- RenameIdents(object = photoreceptors, `13` = "PRCs (14)")
photoreceptors<- RenameIdents(object = photoreceptors, `12` = "PRCs (13)")
photoreceptors<- RenameIdents(object = photoreceptors, `11` = "PRCs (12)")
photoreceptors<- RenameIdents(object = photoreceptors, `10` = "PRCs (11)")
photoreceptors<- RenameIdents(object = photoreceptors, `9` = "PRCs (10)")
photoreceptors<- RenameIdents(object = photoreceptors, `8` = "PRCs (9)")
photoreceptors<- RenameIdents(object = photoreceptors, `7` = "PRCs (8)")
photoreceptors<- RenameIdents(object = photoreceptors, `6` = "PRCs (7)")
photoreceptors<- RenameIdents(object = photoreceptors, `5` = "PRCs (6)")
photoreceptors<- RenameIdents(object = photoreceptors, `4` = "PRCs (5)")
photoreceptors<- RenameIdents(object = photoreceptors, `3` = "PRCs (4)")
photoreceptors<- RenameIdents(object = photoreceptors, `2` = "PRCs (3)")
photoreceptors<- RenameIdents(object = photoreceptors, `1` = "PRCs (2)")
photoreceptors<- RenameIdents(object = photoreceptors, `0` = "PRCs (1)")
markers <- FindAllMarkers(object = photoreceptors, assay = "RNA", only.pos = TRUE, min.pct = 0.01 )
photoreceptors<-ScaleData(photoreceptors, rownames(photoreceptors))
markers %>% group_by(cluster) %>% dplyr::filter(avg_log2FC > 1) %>% slice_head(n = 10) %>% ungroup() -> top10
DoHeatmap(photoreceptors, features = top10$gene) + NoLegend()+scale_fill_viridis()

#Dotplot for genes of interest
DotPlot(object =pljuv_integrated_umap , features = c(genes of interests), col.min=0, scale.by = "size", cols = c("lightgrey", "blue")) +theme_bw() + theme(axis.text.y  = element_text(angle = 0, hjust = 1, size = 15)) +theme(axis.text.x  = element_text(angle = 0, hjust = 0.5, size = 15))+ scale_x_discrete(labels=c(gene names))+coord_flip()
