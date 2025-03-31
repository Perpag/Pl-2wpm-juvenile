library(Seurat)
library(pheatmap)
library(tidyverse)
library(cowplot)
library(dplyr)
library(viridis)
library(ggalluvial)
library(networkD3)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

pm<-pm_ovary_integrated_umap
sp<-sp_ovary
spl<-sp_3dpf_integrated_umap
dm<-dm_ovary
dr<-dr_ovary
hs<-hs_ovary
mm<-mm_ovary
mm_pg<-pg_integrated

pm[["RNA3"]] <- as(object = pm[["RNA"]], Class = "Assay")
DefaultAssay(pm) <- "RNA3"
pm[["RNA"]] <- NULL
pm <- RenameAssays(object = pm, RNA3 = 'RNA')
SaveH5Seurat(pm, filename = "pm_ovary.h5Seurat")
Convert("pm_ovary.h5Seurat", dest = "h5ad")

sp[["RNA3"]] <- as(object = sp[["RNA"]], Class = "Assay")
DefaultAssay(sp) <- "RNA3"
sp[["RNA"]] <- NULL
sp <- RenameAssays(object = sp, RNA3 = 'RNA')
SaveH5Seurat(sp, filename = "sp_adult_ovary.h5Seurat")
Convert("sp_adult_ovary.h5Seurat", dest = "h5ad")

spl[["RNA3"]] <- as(object = spl[["RNA"]], Class = "Assay")
DefaultAssay(spl) <- "RNA3"
spl[["RNA"]] <- NULL
spl <- RenameAssays(object = spl, RNA3 = 'RNA')
SaveH5Seurat(spl, filename = "sp_3df_larva.h5Seurat")
Convert("sp_3df_larva.h5Seurat", dest = "h5ad")

mm[["RNA3"]] <- as(object = mm[["RNA"]], Class = "Assay")
DefaultAssay(mm) <- "RNA3"
mm[["RNA"]] <- NULL
mm <- RenameAssays(object = mm, RNA3 = 'RNA')
SaveH5Seurat(mm, filename = "mm_adult_ovary.h5Seurat")
Convert("mm_adult_ovary.h5Seurat", dest = "h5ad")

dr[["RNA3"]] <- as(object = dr[["RNA"]], Class = "Assay")
DefaultAssay(dr) <- "RNA3"
dr[["RNA"]] <- NULL
dr <- RenameAssays(object = dr, RNA3 = 'RNA')
SaveH5Seurat(dr, filename = "dr_adult_ovary.h5Seurat")
Convert("dr_adult_ovary.h5Seurat", dest = "h5ad")

dm[["RNA3"]] <- as(object = dm[["RNA"]], Class = "Assay")
DefaultAssay(dm) <- "RNA3"
dm[["RNA"]] <- NULL
dm <- RenameAssays(object = dm, RNA3 = 'RNA')
SaveH5Seurat(dm, filename = "dm_adult_ovary.h5Seurat")
Convert("dm_adult_ovary.h5Seurat", dest = "h5ad")

hs[["RNA3"]] <- as(object = hs[["RNA"]], Class = "Assay")
DefaultAssay(hs) <- "RNA3"
hs[["RNA"]] <- NULL
hs <- RenameAssays(object = hs, RNA3 = 'RNA')
SaveH5Seurat(hs, filename = "hs_fetal_ovary.h5Seurat")
Convert("hs_fetal_ovary.h5Seurat", dest = "h5ad")

mm[["RNA3"]] <- as(object = mm[["RNA"]], Class = "Assay")
DefaultAssay(mm) <- "RNA3"
mm[["RNA"]] <- NULL
hs <- RenameAssays(object = mm, RNA3 = 'RNA')
SaveH5Seurat(mm, filename = "mm_pituitary_gland.h5Seurat")
Convert("mm_pituitary_gland.h5Seurat", dest = "h5ad")


#In python

python3

from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd

% bash map_genes.sh --tr1 /directory of P.miniata proteome fasta file --t1 prot --n1 pm --tr2 //directory of species tested proteome fasta file --t2 prot --n2 xx #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison

fn1 = "directory of pm_ovary.h5ad"
fn2 = "directory of species compared to"

filenames = {'pm':fn1,'xx':fn2} #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison

sm = SAMAP(
  filenames,
  f_maps = 'data/maps/',
  save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
)

sm.run(pairwise=True)
samap = sm.samap # SAM object with 2 species stitched together
keys = {'pm':'seurat_clusters','xx':'seurat_clusters'} #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
D.head()
df = MappingTable

print(df)

df.to_csv('mapping_table_sea_star_vs_x_species.csv')

#In R Studio

df <- load mapping scores
data1 <- melt(df)
data1$X <- as.factor(data1$X)

threshold <- 0.1  # threshold value
data1 <- filter(data1, data1$value >= threshold)

# Transform the data for ggalluvial
alluvial_data <- data1 %>%
  mutate(species1 = factor(data1$X, levels = unique(data1$X)),
         species2 = factor(data1$variable, levels = unique(data1$variable)))

ggplot(alluvial_data,
       aes(axis1 = species1, axis2 = species2, y = data1$value)) +
  geom_alluvium(aes(fill = data1$value), width = 0.25) +
  geom_stratum(width = 0.25, fill = "white", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Species 1", "Species 2"), expand = c(0.15, 0.05)) +
  scale_fill_gradient(low = "lightblue", high = "blue", name = "Alignment Score")+
  theme_minimal()+theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )+theme(
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    axis.text.y = element_blank(),  # Remove y-axis labels
    axis.ticks.y = element_blank()  # Remove y-axis ticks
  )

ggsave("Sankey_plot_pm_ovary_vs_x species.tiff", units="in", width=12, height=13
       , dpi=300, compression = 'lzw')

