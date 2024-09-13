```{r}
library(Seurat)
library(Matrix)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)
library(plyr)
library(SeuratData)
library(cowplot)
library(harmony)
library(presto)
library(dittoSeq)
#library(multtest)
library(tidyr)
library(corrplot)

setwd("~/YMJ/Cseq/228a")
```

```{r}
library(RColorBrewer)
display.brewer.all()
c1 =brewer.pal(n = 8, name = "Set1")
c2 =brewer.pal(n = 8, name = "Dark2")
c3 =brewer.pal(n = 8, name = "Paired")
c4 =brewer.pal(n = 8, name = "Accent")

Cols = c(c1,c2, c3, c4)
```

#analysis CITE-seq data

#RNA
```{r}
DefaultAssay(jata) <- 'RNA'
jata <- NormalizeData(jata, normalization.method = "LogNormalize", scale.factor = 10000)
jata <- FindVariableFeatures(jata, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(jata)
jata <- ScaleData(jata, features = all.genes)
jata <- RunPCA(jata, features = VariableFeatures(object = jata))
ElbowPlot(jata)
```
![image](https://github.com/user-attachments/assets/40614cd0-a9e4-498b-8581-d7d1b010408a)


```{r, fig.width=5, fig.height=5}
DefaultAssay(jata) <- 'RNA'
jata <- FindNeighbors(jata, reduction = "pca", dims = 1:8)
jata <- FindClusters(jata,  reduction = "pca", resolution = 0.8)
jata <- RunUMAP(jata, reduction = "pca", dims = 1:15, reduction.name = "rna.umap")
DimPlot(jata, reduction = "rna.umap", label = T,repel = T, group.by = "RNA_snn_res.0.8")
DimPlot(jata, reduction = "rna.umap", label = T, repel=T, group.by = "celltype.l2", cols=Cols) + NoLegend()
DimPlot(jata, reduction = "rna.umap", label = T, repel=T, group.by = "celltype.l1", cols=Cols) + NoLegend()
```
![image](https://github.com/user-attachments/assets/3befb319-82e7-42ce-bd79-d7a1b9377510)
Azimuth cell annotation
![image](https://github.com/user-attachments/assets/7aac1d5a-4881-4ef5-83fd-da8c72fb0bf5)

![image](https://github.com/user-attachments/assets/b55c73d0-9cb1-4db2-bde5-37cfcafd487f)


```{r , fig.width=12, fig.height=4}
DefaultAssay(jata) <- 'RNA'
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(jata), 10)
# plot variable features with and without labels
plot3 <- VariableFeaturePlot(jata)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot3 + plot4
```
![image](https://github.com/user-attachments/assets/20a7ce1a-1dc8-4ef0-8190-7c740ec6615b)


#ADT
```{r}
DefaultAssay(jata) <- 'ADT'
VariableFeatures(jata) <- rownames(jata[["ADT"]])
jata <- NormalizeData(jata, normalization.method = "LogNormalize", scale.factor = 10000)
all.ADT <- rownames(jata[["ADT"]])
jata <- ScaleData(jata, features = all.ADT)
jata <- RunPCA(jata, features = VariableFeatures(jata), reduction.name = "apca")
ElbowPlot(jata, reduction = "apca")
```
![image](https://github.com/user-attachments/assets/69803fd8-b759-4641-82ca-5dda421ab2b0)


```{r, fig.width=5, fig.height=5}
DefaultAssay(jata) <- 'ADT'
jata <- FindNeighbors(jata, reduction = "apca", dims = 1:10)
jata <- FindClusters(jata,  reduction = "apca", resolution = 1)
jata <- RunUMAP(jata, reduction = "apca", dims = 1:20, reduction.name = "adt.umap")
DimPlot(jata, reduction = "adt.umap", label = T)
DimPlot(jata, reduction = "adt.umap", label = T, repel=T, group.by = "celltype.l2", cols=Cols2) + NoLegend()
```
![image](https://github.com/user-attachments/assets/8ec8f119-ddde-4e3d-bd9f-d11c0f73293e)

![image](https://github.com/user-attachments/assets/c25dd8dc-b2c1-46be-a3c3-0248ea249a5b)


#WNN
```{r, fig.width=7, fig.height=3}
DefaultAssay(jata) <- 'RNA'
jata <- FindMultiModalNeighbors(jata, reduction.list = list("pca", "apca"), dims.list = list(1:8, 1:8), modality.weight.name = "RNA.weight")

jata <- FindClusters(jata, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)
jata <- RunUMAP(jata, nn.name = "weighted.nn", reduction.name = "wnn.umap")
```

```{r, fig.width=20, fig.height=16}

Idents(jata) <- jata$wsnn_res.0.5
plo1 <- DimPlot(jata, reduction = 'wnn.umap', label = TRUE, repel = T, group.by = "celltype.l2") + NoLegend()
plo2 <- DimPlot(jata, reduction = 'wnn.umap', label = TRUE, repel = T) + NoLegend()
plo3 <- DimPlot(jata, reduction = 'wnn.umap', label = TRUE, group.by = "Batch") + NoLegend()

Idents(jata) <- jata$RNA_snn_res.0.8
plo4 <- DimPlot(jata, reduction = 'rna.umap', label = TRUE, repel = T, group.by = "celltype.l2") + NoLegend()
plo5 <- DimPlot(jata, reduction = 'rna.umap', label = TRUE, repel = T) + NoLegend()
plo6 <- DimPlot(jata, reduction = 'rna.umap', label = TRUE, group.by = "Batch") + NoLegend()

Idents(jata) <- jata$ADT_snn_res.1
plo7 <- DimPlot(jata, reduction = 'adt.umap', label = TRUE, repel = T, group.by = "celltype.l2") + NoLegend()
plo8 <- DimPlot(jata, reduction = 'adt.umap', label = TRUE, repel = T) + NoLegend()
plo9 <- DimPlot(jata, reduction = 'adt.umap', label = TRUE, group.by = "Batch") + NoLegend()

plo1 + plo2 + plo3 + plo4 + plo5 + plo6 + plo7 + plo8 + plo9
```
![image](https://github.com/user-attachments/assets/888f19b7-a44c-4010-aa08-55f5531ca91e)
