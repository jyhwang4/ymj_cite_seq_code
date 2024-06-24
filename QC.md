# Qulity Control


```{r , fig.width=9, fig.height=6}
jata[["percent.mt"]] <- PercentageFeatureSet(jata, pattern = "^MT-")

VlnPlot(jata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "orig.ident")
```
![image](https://github.com/jyhwang4/ymj_cite_seq_code/assets/59998490/b5ecbfd1-4525-461d-aac3-5629faf20c5e)

```{r , fig.width=8, fig.height=3}
plot1 <- FeatureScatter(jata, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "orig.ident")
plot2 <- FeatureScatter(jata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")
plot1 + plot2
```
![image](https://github.com/jyhwang4/ymj_cite_seq_code/assets/59998490/df1f856a-082d-42bc-9d87-fcad72e365e3)

```{r}
jata <- subset(jata, subset = nFeature_RNA < 5000 & nCount_RNA < 30000 & percent.mt < 17)
```
