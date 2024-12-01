if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#Paquetes
list.of.packages = c("tximeta", "tximport", "limma", "edgeR", "tidyverse", "org.Mm.eg.db", "statmod", "pheatmap", "ggplotify", "ggrepel", "SummarizedExperiment", "patchwork", "xlsx", "ragg", "OmnipathR", "clusterProfiler")
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
files = data.frame(read.table(file="Metadata/metadata.txt", header = TRUE, stringsAsFactors = F))[c(18:35),]
se = tximeta(files)                   
se <- addExons(se)
gse <- summarizeToGene(se, assignRanges="abundant") #Asigna rangos según la isoforma más abundante de los transcritos en vez de desde la isoforma más en 5' a la isoforma más en 3'. Los creadores lo recomiendan.
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
data$DE[data$logFC > 0 & data$PValue < 0.05] <- "UP"
data$DE[data$logFC < 0 & data$PValue < 0.05] <- "DOWN"
ggplot(data = data, aes(x=logFC, y =-log10(PValue), col=DE, label = ifelse(abs(logFC) >= 2 & PValue <0.05, as.character(symbol),  ''))) + geom_point() +
  geom_text_repel(hjust = 0, nudge_x = 0.1, color = "black") +
  theme_minimal() +
  geom_hline(yintercept=-log10(0.05))
setwd("/home/guille/RNA/RNA_2024/EXPA/experimento_A/Results/")
experimento = "exp_A_colon"
dir.create(experimento, showWarnings = FALSE)
if (!file.exists(paste0(experimento, "/",experimento, ".xlsx"))) {
  write.xlsx2(x = res_corrected[!is.na(res_corrected$table$gene_name) & res_corrected$table$PValue <= 0.05 ,c(1,7,13,16,17)], file = paste0(experimento, "/",experimento, ".xlsx"), col.names = T, row.names = F, append = TRUE)
} else {
  print("Ya existe")
}
res_corrected[1,c(1,7,13,16,17)]
setwd("/home/guille/RNA/RNA_2024/EXPA/experimento_A")

logcpm = cpm(y, log=TRUE)
rownames(logcpm) = y$genes$symbol #Nombre de genes
colnames(logcpm) =  paste(y$samples$group, y$samples$names, sep = "-")

DEG = res_corrected$table[res_corrected$table$PValue < 0.01 & abs(res_corrected$table$logFC) > 0 , ]
DEG_selection = logcpm[na.omit(DEG$symbol),]


q =  as.ggplot(pheatmap(DEG_selection, scale = "row", 
                        clustering_method = "complete",
                        display_numbers = F,
                        border_color = NA, cluster_cols = T, cutree_cols = 2, cutree_rows = 2, show_rownames = T,
                        #annotation_col = annotation, show_rownames = F, annotation_names_col = F,
                        #annotation_row = setNames(data.frame(Cluster = as.factor(cutree(q$tree_row, k=2))), "Cluster"), 
                        annotation_names_row = F,
                        legend_labels = F,))
ggsave(filename="Results/exp_A_colon/heatmap.png", dpi = "print", bg = "white", scale = 1.80)

options(enrichplot.colours = c("#FF3030CC", "#1E90FFCC"))

original_gene_list <- res_corrected$table$logFC
names(original_gene_list) <- res_corrected$table$symbol
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=na.omit(gene_list[!duplicated(names(gene_list))]), 
             ont ="BP", 
             keyType = "SYMBOL",
             nPermSimple = 1000,
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")
dotplot(gse,title = "Biological Process", split=".sign", ) + facet_grid(.~.sign)
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 10)
ridgeplot(gse, showCategory = 15) + labs(x = "Enrichment Distribution", title = "Biological Process") + theme(plot.title = element_text(hjust = 0.5)) 

#KEGG
kegg_gene_list <- res_corrected$table$logFC
names(kegg_gene_list) <- res_corrected$table$entrezid
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list[!duplicated(names(kegg_gene_list))],
               organism     = kegg_organism,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")
kk3 = setReadable(kk2, org.Mm.eg.db, keyType = "ENTREZID")

dotplot(kk2,  title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
cnetplot(kk3, categorySize="pvalue", foldChange=gene_list, )
ridgeplot(kk2, showCategory = 15) + labs(x = "Enrichment Distribution", title = "KEGG Pathways") + theme(plot.title = element_text(hjust = 0.5))
