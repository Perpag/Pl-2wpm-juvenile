#For Samap analysis
#In R studio

library(reshape2)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

pl<-pljuv_integrated_umap
sp<-sp_3dpf_integrated_umap


pl[["RNA3"]] <- as(object = pl[["RNA"]], Class = "Assay")
DefaultAssay(pl) <- "RNA3"
pl[["RNA"]] <- NULL
pl <- RenameAssays(object = pl, RNA3 = 'RNA')
SaveH5Seurat(pl, filename = "pljuv.h5Seurat")
Convert("pljuv.h5Seurat", dest = "h5ad")

sp[["RNA3"]] <- as(object = sp[["RNA"]], Class = "Assay")
DefaultAssay(sp) <- "RNA3"
sp[["RNA"]] <- NULL
sp <- RenameAssays(object = sp, RNA3 = 'RNA')
SaveH5Seurat(sp, filename = "sp3dpf.h5Seurat")
Convert("sp3dpf.h5Seurat", dest = "h5ad")

#In python

cd samap_directory

python3

from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd

% bash map_genes.sh --tr1 /directory of P.lividus proteome fasta file --t1 prot --n1 pl --tr2 //directory of S. purpuratus proteome fasta file --t2 prot --n2 sp

fn1 = "directory of pljuv.h5ad"
fn2 = "directory of sp3dpf.h5ad"

filenames = {'pl':fn1,'sp':fn2} 

sm = SAMAP(
  filenames,
  f_maps = 'data/maps/',
  save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
)

sm.run(pairwise=True)
samap = sm.samap # SAM object with 2 species stitched together
keys = {'pl':'seurat_clusters','sp':'seurat_clusters'}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
D.head()
df = MappingTable

print(df)

df.to_csv('mapping table Pl 2wpm juvenile_sp_3dpf_larva.csv')



#SAMAP for Pl object with merged neuronal clusters into one

cd samap_directory
python3

from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd


fn1 = "/Users/periklespaganos/Downloads/Juvenile_paper submission/SAMap_Pl_juvenile vs Sp_larva files/h5ad_files/pljuv_merged_neurons.h5ad"
fn2 = "/Users/periklespaganos/Downloads/Juvenile_paper submission/SAMap_Pl_juvenile vs Sp_larva files/h5ad_files/sp3df.h5ad"
filenames = {'pl':fn1,'sp':fn2}

sm = SAMAP(
  filenames,
  f_maps = 'example_data/maps/',
  save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
)

sm.run(pairwise=True)
samap = sm.samap # SAM object with 2 species stitched together

keys = {'sp':'seurat_clusters','pl':'merged_clusters'}

D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
D.head()
df = MappingTable

print(df)

df.to_csv('mapping table Pl merged_neurons_sp_3dpf_larva.csv')


# Plotting of the mapping tables as heatmaps In R Studio
library(ggplot2)
library(reshape2)

#load mapping table (df)

data1 <- melt(df)
data1$X <- as.factor(data1$X)
levels(data1$X)
data1$X <- factor(data1$X, levels = levels(data1$X)[c(7,4,1, 13, 21, 14, 9, 8, 11, 5, 16, 3, 12, 15, 17, 18, 19, 2, 6, 10, 20)])

ggplot(data1, aes(x = X, y = variable, fill = value)) +geom_tile(color = "black") +scale_fill_gradient2(low = "blue", mid = "#FFFFCC",high = "#FF0000", midpoint = 0.25)+coord_fixed() +guides(fill = guide_colourbar(title = "Alignment score"))+labs(x = "Strongylocentrotus purpuratus 3 dpf larva", y = "Paracentrotus lividus 2wpm") +RotatedAxis()
ggsave("asadultntop0.tiff", units="in", width=9, height=9, dpi=300, compression = 'lzw')
