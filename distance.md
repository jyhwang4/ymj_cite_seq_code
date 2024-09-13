
```{r}
cluster_info <- data.frame(Cluster = jata$celltype.l2,
                           UMAP1 = jata@reductions[["rna.umap"]]@cell.embeddings[,1],
                           UMAP2 = jata@reductions[["rna.umap"]]@cell.embeddings[,2])
                           
library(data.table)
cluster_data <- cluster_info[cluster_info$Cluster %in% c("CD8 Naive", "CD8 TEM"), ]
centroid_CD8_Naive <- colMeans(cluster_data[cluster_data$Cluster == "CD8 Naive", c("UMAP1", "UMAP2")])
centroid_CD8_TEM <- colMeans(cluster_data[cluster_data$Cluster == "CD8 TEM", c("UMAP1", "UMAP2")])
centroid_distance <- sqrt(sum((centroid_CD8_Naive - centroid_CD8_TEM)^2))
cat("mRNA Centroid distance between CD8 Naive and CD8 TEM:", centroid_distance, "\n")



cluster_info <- data.frame(Cluster = jata$celltype.l2,
                           UMAP1 = jata@reductions[["adt.umap"]]@cell.embeddings[,1],
                           UMAP2 = jata@reductions[["adt.umap"]]@cell.embeddings[,2])
                           
library(data.table)
cluster_data <- cluster_info[cluster_info$Cluster %in% c("CD8 Naive", "CD8 TEM"), ]
centroid_CD8_Naive <- colMeans(cluster_data[cluster_data$Cluster == "CD8 Naive", c("UMAP1", "UMAP2")])
centroid_CD8_TEM <- colMeans(cluster_data[cluster_data$Cluster == "CD8 TEM", c("UMAP1", "UMAP2")])
centroid_distance <- sqrt(sum((centroid_CD8_Naive - centroid_CD8_TEM)^2))
cat("Protein Centroid distance between CD8 Naive and CD8 TEM:", centroid_distance, "\n")



```
mRNA Centroid distance between CD8 Naive and CD8 TEM: 12.26637 
Protein Centroid distance between CD8 Naive and CD8 TEM: 3.528299 

```{r}
cluster_info <- data.frame(Cluster = jata$celltype.l2,
                           UMAP1 = jata@reductions[["rna.umap"]]@cell.embeddings[,1],
                           UMAP2 = jata@reductions[["rna.umap"]]@cell.embeddings[,2])
                           
library(data.table)
cluster_data <- cluster_info[cluster_info$Cluster %in% c("CD4 Naive", "CD4 TEM"), ]
centroid_CD8_Naive <- colMeans(cluster_data[cluster_data$Cluster == "CD4 Naive", c("UMAP1", "UMAP2")])
centroid_CD8_TEM <- colMeans(cluster_data[cluster_data$Cluster == "CD4 TEM", c("UMAP1", "UMAP2")])
centroid_distance <- sqrt(sum((centroid_CD8_Naive - centroid_CD8_TEM)^2))
cat("mRNA Centroid distance between CD4 Naive and CD4 TEM:", centroid_distance, "\n")



cluster_info <- data.frame(Cluster = jata$celltype.l2,
                           UMAP1 = jata@reductions[["adt.umap"]]@cell.embeddings[,1],
                           UMAP2 = jata@reductions[["adt.umap"]]@cell.embeddings[,2])
                           
library(data.table)
cluster_data <- cluster_info[cluster_info$Cluster %in% c("CD4 Naive", "CD4 TEM"), ]
centroid_CD8_Naive <- colMeans(cluster_data[cluster_data$Cluster == "CD4 Naive", c("UMAP1", "UMAP2")])
centroid_CD8_TEM <- colMeans(cluster_data[cluster_data$Cluster == "CD4 TEM", c("UMAP1", "UMAP2")])
centroid_distance <- sqrt(sum((centroid_CD8_Naive - centroid_CD8_TEM)^2))
cat("Protein Centroid distance between CD4 Naive and CD4 TEM:", centroid_distance, "\n")
```
mRNA Centroid distance between CD4 Naive and CD4 TEM: 6.299101 
Protein Centroid distance between CD4 Naive and CD4 TEM: 3.303466 

```{r}
#CD8
distmat <- matrix(0, nrow = 2, ncol = 3)
distmat <- as.data.frame(distmat)

distmat[1,1] <- "mRNA"
distmat[2,1] <- "Protein"
distmat[1,2] <- 12.26637
distmat[2,2] <- 3.528299
distmat[1,3] <- "CD8+ T cells"
distmat[2,3] <- "CD8+ T cells"

#CD4
distmat2 <- matrix(0, nrow = 2, ncol = 3)
distmat2 <- as.data.frame(distmat2)

distmat2[1,1] <- "mRNA"
distmat2[2,1] <- "Protein"
distmat2[1,2] <- 6.299101
distmat2[2,2] <- 3.303466
distmat2[1,3] <- "CD4+ T cells"
distmat2[2,3] <- "CD4+ T cells"
```

```{r, fig.height=3, fig.width=2}
pcol1 <- c("#FD8008")

ggplot(distmat, aes(x=V1, y=V2, group = V3)) +
 geom_line(color=pcol1, size = 1) +
 geom_point(shape = 15, color=pcol1, size = 3) +
 scale_y_continuous(limits = c(0, 13), breaks = seq(0, 13, by = 2)) +
 theme(panel.grid = element_blank()) +
 theme(axis.line = element_line(color = "black", size = 1)) +
  theme(axis.text.x = element_text(color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
 labs(x = "Assay", y = "Centroid distance between 
      CD8 Naive and CD8 TEM clusters")
```
![image](https://github.com/user-attachments/assets/72071857-fd5f-45c9-941a-e840a2360aa3)


```{r, fig.height=3, fig.width=2}
pcol2 <- c("#105B08")

ggplot(distmat2, aes(x=V1, y=V2, group = V3)) +
 geom_line(color=pcol2, size = 1) +
 geom_point(shape = 15, color=pcol2, size = 3) +
 scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2)) +
 theme(panel.grid = element_blank()) +
 theme(axis.line = element_line(color = "black", size = 1)) +
  theme(axis.text.x = element_text(color = "black")) +
  theme(axis.text.y = element_text(color = "black")) +
 labs(x = "Assay", y = "Centroid distance between 
      CD4 Naive and CD4 TEM clusters")
```

![image](https://github.com/user-attachments/assets/fc11fc82-4294-48ae-ae30-2cc9e377a1d5)
