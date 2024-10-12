if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#Paquetes
list.of.packages = c("tximeta", "tximport", "limma", "edgeR", "tidyverse", "org.Mm.eg.db", "statmod", "pheatmap", "ggplotify", "ggrepel", "SummarizedExperiment", "patchwork", "xlsx", "ragg", "OmnipathR")
#Instalación por CRAN o Bioconductor
new.packages = list.of.packages[!(list.of.packages %in% installed.packages())]
if(length(new.packages)> 0) {
  for(i in new.packages) {
    if(i %in% available.packages()[,1]){ #Chequea si el paquete está en un repositorio CRAN y lo instala
      install.packages(i,dependencies=TRUE)
    }else {BiocManager::install(i, update = TRUE, ask = FALSE, version = BiocManager::version()) #Instala por BiocManager si no está en un repositorio CRAN
    }}
}
invisible(lapply(list.of.packages, FUN=library, character.only=TRUE))
rm(list.of.packages, new.packages)
setwd("/home/guille/RNA/RNA_2024/EXPA/experimento_A")
files = data.frame(read.table(file="Metadata/metadata.txt", header = TRUE, stringsAsFactors = F))[c(18:35),]
suppressPackageStartupMessages(library(SummarizedExperiment))
se = tximeta(files)                   
se <- addExons(se)
gse <- summarizeToGene(se, assignRanges="abundant") #Asigna rangos según la isoforma más abundante de los transcritos en vez de desde la isoforma más en 5' a la isoforma más en 3'. Los creadores lo recomiendan.
gse <- addIds(gse, "GENENAME", gene=TRUE)
gse <- addIds(gse, "SYMBOL", gene=TRUE)
gse <- addIds(gse, "ENTREZID", gene = TRUE)
y <- makeDGEList(gse)
rm(se)
sampleinfo <- read.delim(file ="Metadata/sampleinfo_colon.txt", sep = " ", header = T)
group = sampleinfo[,2]
group = factor(group)
group = relevel(group, ref = "WT")
y$samples$group = group
y$samples
design = model.matrix(~0 + group)
colnames(design) = levels(group)
print(design)
keep <- filterByExpr(y, design)
print(summary(keep))
y <- y[keep, ]
points <- c(1:2) # Formas
colors <- c(1:2) #Colores
mds = plotMDS(y, col=colors[group], pch=points[group])
mds_data <- data.frame(
  Dim1 = mds$x,
  Dim2 = mds$y,
  Sample = colnames(mds$distance.matrix.squared),
  Group = as.factor(y$samples$group)
)
axis_label = round(mds$var.explained[1:2]*100)
library(ggConvexHull)
mds_data %>% ggplot(aes(x = Dim1, y = Dim2, col = Group, shape = Group)) +
  geom_point(size = 2) +
  labs(title = "MDS Analysis", x = paste0("Dimension 1: ", axis_label[1], "%"), y = paste0("Dimension 2: ", axis_label[2], "%"), ) +
  theme_minimal() +
  scale_shape_manual(values = rep(15:17, len = 7)) + 
  geom_convexhull(aes(fill=Group, 
                      colour= Group),
                  alpha = 0.2) + 
  theme(plot.title = element_text(hjust = 0.5))
y <- estimateDisp(y, design, method = "ROBUST")
plotBCV(y)
fit <- glmQLFit(y, design, method = "ROBUST")
plotQLDisp(fit)
contrast = makeContrasts(
 "KO - WT",
  levels=design)
res <- glmQLFTest(fit, contrast = contrast)
res_corrected = topTags(res, n = Inf)
head(res_corrected)
data = res_corrected$table
data$DE <- "NO"
data$DE[data$logFC > 1 & data$PValue < 0.05] <- "UP"
data$DE[data$logFC < -1 & data$PValue < 0.05] <- "DOWN"
ggplot(data = data, aes(x=logFC, y =-log10(PValue), col=DE, label = ifelse(abs(logFC) >= 2 & PValue <0.05, as.character(SYMBOL),  ''))) + geom_point() +
  geom_text_repel(hjust = 0, nudge_x = 0.1, color = "black") +
  theme_minimal() +
  labs(title = as.character(res$comparison)) +
  geom_hline(yintercept=-log10(0.05)) +
  geom_vline(xintercept=1) +
  geom_vline(xintercept=-1)
