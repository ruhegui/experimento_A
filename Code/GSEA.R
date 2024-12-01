library("enrichplot")
library("GOSemSim")
library("clusterProfiler")

GSEA_data <- res_corrected$table %>% 
  dplyr::filter(!is.na(FDR)) %>%
  dplyr::arrange(desc(logFC))
log2FC = GSEA_data %>% pull(logFC, name = gene_name)
log2FC = log2FC[unique(names(log2FC))]
rm(GSEA_data)
set.seed(123) #RANDOM ORDERING => USE THE SAME SEED
gsea_go <- gseGO(geneList = log2FC,
                 OrgDb = org.Mm.eg.db,
                 ont = "BP",
                 keyType = "SYMBOL",
                 seed=TRUE, 
                 eps = 0)
s_gsea = clusterProfiler::simplify(gsea_go)
s_gsea@result %>% head()
s_gsea@result %>%
  dplyr::select(ID, Description, enrichmentScore, NES, pvalue, p.adjust) %>%
  write.xlsx2(x = . , file = paste0("Results", "/GSEA_yeyuno.xlsx"), col.names = T, row.names = F)
ggplot(s_gsea, showCategory=20, aes(NES, fct_reorder(Description, NES),
                                     fill=p.adjust)) +
  geom_col() +
  geom_vline(xintercept=0, linetype="dashed", color="black",size=1)+
  scale_fill_gradientn(colours=c("#b3eebe","#46bac2", "#371ea3"),
                       guide=guide_colorbar(reverse=TRUE))+
  scale_x_continuous(expand=c(0,0))+
  scale_y_discrete(labels = scales::label_wrap(40)) +
  theme_bw() + 
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("GSEA with GO")
ggsave(paste0("Results/GSEA_colon.pdf"), bg = "white", scale = 1.3)
pairwise_termsim(s_gsea, method = "Wang", semData = godata('org.Mm.eg.db', ont="BP")) %>%
  treeplot(showCategory = 20, label_format_tiplab = 50)
ggsave(paste0("Results/GSEA_tree_colon.pdf"), bg = "white", scale = 1.5)
