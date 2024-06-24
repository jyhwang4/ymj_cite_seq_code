# ymj_cite_seq_code

# load data
```{r}
barcode_file <- "~/YMJ/Cseq/228a/3P_RNA/barcodes.tsv.gz"
feature_file <- "~/YMJ/Cseq/228a/3P_RNA/features.tsv.gz"
matrix_file <- "~/YMJ/Cseq/228a/3P_RNA/matrix.mtx.gz"
data <- Read10X(data.dir = "~/YMJ/Cseq/228a/3P_RNA")


barcode_file <- "~/YMJ/Cseq/228a/3P_ADT/barcodes.tsv.gz"
feature_file <- "~/YMJ/Cseq/228a/3P_ADT/features.tsv.gz"
matrix_file <- "~/YMJ/Cseq/228a/3P_ADT/matrix.mtx.gz"
adata <- Read10X(data.dir = "~/YMJ/Cseq/228a/3P_ADT")

all.equal(colnames(data), colnames(adata), colnales(bdata))
```
# Create Suerat object
```{r}
jata <- CreateSeuratObject(counts = data, project = "3P-RNA", min.cells = 200, min.features = 500)

adt_assay <- CreateAssayObject(counts = adata)

jata[["ADT"]] <- adt_assay

Assays(jata)

rownames(jata[["ADT"]])

mdata <- read.csv(file = "~/YMJ/Cseq/228a/GSE164378_sc.meta.data_3P.csv.gz", row.names = 1)

jata <- AddMetaData(jata, metadata = mdata)

DefaultAssay(jata)
```

# Data assays check
```{r}
# Switch the default to ADT
DefaultAssay(jata) <- "ADT"
DefaultAssay(jata)
```

```{r}
# default assay is RNA
DefaultAssay(jata) <- "RNA"
DefaultAssay(jata)
```
# Batch correction
```{r}
jata <- subset(jata, Batch %in% c("Batch2"))

set.seed(1230)  
```
