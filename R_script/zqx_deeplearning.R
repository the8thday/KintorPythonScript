#! /urs/bin/env Rscript
# Title     : zqx protein data from deeplearning
# Objective : data analysis
# Created by: congliu
# Created on: 2021/9/2

library(tidyverse)
# library(GSVA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
# library(Pi)


gcnnet <- read_csv('/Users/congliu/OneDrive/kintor/qianx/deep_dock/deep_learning_to_kintor/model_GCNNet_chembl_enzyme_200k_all_protein.csv')
gcnet_immune <- read_delim('/Users/congliu/OneDrive/kintor/qianx/deep_dock/deep_learning_to_kintor/model_GCNNet_chembl_enzyme_immune.txt',
                           delim = '\t')
gcnet_immune_filter <- read_delim('/Users/congliu/OneDrive/kintor/qianx/deep_dock/deep_learning_to_kintor/model_GCNNet_chembl_enzyme_immune_cutoff_-0.97.txt',
                           delim = '\t')
docking <- read_delim('/Users/congliu/OneDrive/kintor/qianx/deep_dock/reverse_docking_to_kintor/docking_results_protein.txt',
                      delim = '\t')
gcn_54_path <- '/Users/congliu/OneDrive/kintor/qianx/deep_dock/deep_learning_to_kintor/model_GCNNet_chembl_enzyme_immune_cutoff_-0.97_train_491_protein.txt'
gcn_54 <- read_delim(gcn_54_path, delim = '\t')

# uniprt_id <- read_delim('/Users/congliu/OneDrive/kintor/qianx/deep_dock/deep_learning_to_kintor/uniprot-immune-id.txt',
#                         delim = '\t')
uniprt_id <- read_delim('/Users/congliu/OneDrive/kintor/qianx/deep_dock/deep_learning_to_kintor/all_uniprot_symbol.txt',
                        delim = '\t') %>% 
  # dplyr::mutate(Symbol = str_split(`Gene names`, ' ')[[1]][[1]])
  separate(col = `Gene names`, into = c('Symbol', 'rest'), sep = ' ', remove = F, extra = "merge")
two_19 <- read_delim('/Users/congliu/OneDrive/kintor/qianx/deep_dock/reverse_docking_to_kintor/two_method_overlap.txt',
                                  delim = '\t')

# different cohert
cohort <- 'deep2w'
outpath <- '/Users/congliu/OneDrive/kintor/qianx/deep_dock/'

get_gene <- function(dataset){
  dataset %>%
    left_join(uniprt_id, by = c('uniprot_id'='Entry')) %>% 
    separate(col = symbol, into = c('symbol', NA), sep = '_') %>%
    arrange(desc(score)) %>%
    dplyr::select(Symbol) %>% distinct() %>% pull(Symbol)
}
genelist <- get_gene(gcnnet %>% dplyr::filter(score>-0.97)) # 注意修改
# genelist <- two_19 %>% dplyr::select(Symbol) %>% distinct() %>% pull(Symbol)

gene.df <- bitr(genelist,
                fromType = 'SYMBOL',
                toType = c('ENTREZID', 'ENSEMBL'),
                OrgDb = org.Hs.eg.db,
                drop = TRUE
                )

# go 富集分析
enrich.go <- enrichGO(
  gene = genelist,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 1,
  qvalueCutoff  = 1
)
# enrich.go@pvalueCutoff <- 1
# enrich.go@qvalueCutoff <- 1

egosim <- simplify(enrich.go,
                   cutoff = 0.7,
                   by = 'p.adjust',
                   select_fun = min,
                   measure = 'Wang'
)

pp <- dotplot(enrich.go,
              font.size = 12,
              color = 'p.adjust',
              showCategory=20,
              title = 'GO(BP) Analysis'
              ) +
  # scale_color_gradient2() +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

pp2 <- barplot(enrich.go,
        font.size = 12
        ) +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

(pp3 <- cnetplot(enrich.go, 
                 node_label = 'all'))

ego2 <- enrichplot::pairwise_termsim(enrich.go)
(pp4 <- enrichplot::treeplot(ego2,
                             hclust_method = "ward.D"
                             ))

write_delim(enrich.go@result, 
            file = file.path(outpath, paste0(cohort, '_goResult.txt')), 
            delim = '\t')
write_delim(egosim@result, 
            file = file.path(outpath, paste0(cohort, '_goSim.txt')), 
            delim = '\t')
ggsave(file = file.path(outpath, paste0(cohort, '_goResult.pdf')), 
       plot = pp, width = 10, height = 10)
ggsave(file = file.path(outpath, paste0(cohort, '_goCnet.pdf')), 
       plot = pp3, width = 10, height = 10)
ggsave(file = file.path(outpath, paste0(cohort, '_goTree.pdf')), 
       plot = pp4, width = 10, height = 10)

# kegg
kegg <- enrichKEGG(
  gene = gene.df$ENTREZID,
  keyType = 'ncbi-geneid',
  organism = 'hsa',
  pAdjustMethod = 'fdr',
  pvalueCutoff = 0.5,
  qvalueCutoff = 0.5,
  use_internal_data = T
)
kegg@pvalueCutoff <- 1
kegg@qvalueCutoff <- 1

keggx <- setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')

kegg_p <- barplot(kegg,
        font.size = 15,
        title = 'KEGG Analysis'
        )
write_delim(kegg@result, 
            file = file.path(outpath, paste0(cohort, '_KEGG.txt')), delim = '\t')
write_delim(keggx@result, 
            file = file.path(outpath, paste0(cohort, '_KEGGX.txt')), delim = '\t')
ggsave(file.path(outpath, paste0(cohort, '_KEGG.pdf')), 
       plot = kegg_p, width = 10, height = 10)

# browseKEGG(kegg, pathID = 'hsa04080')

# WikiPathways, Reactome
ewiki <- enrichWP(gene.df$ENTREZID, 
         organism = "Homo sapiens",
         pvalueCutoff = 0.5,
         qvalueCutoff = 0.5
         )
ewiki_p <- dotplot(ewiki,
                   font.size = 12,
                   color = 'p.adjust',
                   showCategory=20,
                   title = 'WikiPathway ORA Analysis'
) +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})
ggsave(file.path(outpath, paste0(cohort, '_WikiPathway.pdf')), 
       plot = ewiki_p, width = 10, height = 10)
ewiki <- setReadable(ewiki, OrgDb = org.Hs.eg.db)
write_delim(ewiki@result, 
            file = file.path(outpath, paste0(cohort, '_WikiPathway.txt')), delim = '\t')

eReactome <- ReactomePA::enrichPathway(gene=gene.df$ENTREZID, 
                                       pvalueCutoff = 0.05, 
                                       readable=TRUE,
                                       organism = "human",
                                       pAdjustMethod = "BH",
                                       qvalueCutoff = 0.2
                                       )
react_p <- dotplot(eReactome,
              font.size = 12,
              color = 'p.adjust',
              showCategory=20,
              title = 'Reactome Analysis'
) +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})

# viewPathway("E2F mediated regulation of DNA replication", 
#             readable = TRUE, 
#             foldChange = geneList)
ggsave(file.path(outpath, paste0(cohort, '_Reactome.pdf')), 
       plot = react_p, width = 10, height = 10)
write_delim(eReactome@result, 
            file = file.path(outpath, paste0(cohort, '_Reactome.txt')), delim = '\t')


y <- gsePathway(gene.df$ENTREZID, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)

# DOSE
enrich.do <- DOSE::enrichDO(
  gene = gene.df$ENTREZID,
  ont = 'DO',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

do_p <- dotplot(enrich.do,
                   font.size = 12,
                   color = 'p.adjust',
                   showCategory=20,
                   title = 'DO Analysis'
) +
  scale_y_discrete(labels = function(x) {str_wrap(x, width = 40)})
ggsave(file.path(outpath, paste0(cohort, '_DO.pdf')), 
       plot = do_p, width = 10, height = 10)
write_delim(enrich.do@result, 
            file = file.path(outpath, paste0(cohort, '_DO.txt')), delim = '\t')


# gsea 富集 免疫通路
get_id <- function(dataset){
  dataset %>%
    left_join(uniprt_id, by = c('uniprot_id'='Entry')) %>% 
    dplyr::select(score, Symbol) %>% 
    left_join(gene.df, by = c('Symbol'='SYMBOL')) %>% 
    dplyr::select(ENTREZID, score) %>% 
    arrange(desc(score)) %>% distinct()
}
ss <- get_id(gcnnet %>% dplyr::filter(score>-0.97)) # 注意修改数据集
gene_sort <- ss %>% pull(score)
names(gene_sort) <- ss$ENTREZID

c7 <- msigdbr::msigdbr(species = "Homo sapiens", category = 'C7') %>% 
  dplyr::select(gs_name, entrez_gene) %>% as.data.frame()
colnames(c7)<-c("ont","gene")

gsea <- GSEA(gene_sort, 
             TERM2GENE = c7, 
             pvalueCutoff = 0.5)
em2 <- setReadable(gsea, 
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID")
em2.df <- as.data.frame(em2)
em2.df <- em2.df %>%
  dplyr::select(ID,NES,setSize,p.adjust)
GSEA_result2 <- merge(c7,em2.df, by.x='ont', by.y = 'ID')
id <- 1
anno <- em2[id, c("NES", "pvalue", "p.adjust")]
lab <- paste0(names(anno), "=",  round(anno, 3), collapse="\n")

gsea_p1 <- enrichplot::gseaplot(gsea, 1,
                                title = gsea$Description[1]
                                ) +
  annotate("text", 0.8, 0.85,
           label = lab, hjust=0, vjust=0)

gsea_p <- enrichplot::gseaplot2(gsea, 1:2,
                                title = 'GSEA Analysis'
                                )
ridgeplot(gsea)
ggsave(file = file.path(outpath, paste0(cohort, '_GSEA.pdf')), 
       plot = gsea_p, width = 10, height = 10)

# gsea kegg
gseKEGG.res <- gseKEGG(gene_sort, 
                       organism = "hsa",
                       keyType = "kegg",
                       minGSSize = 5, 
                       maxGSSize = 500,
                       pvalueCutoff=1
                       )
gseKEGG_p <- enrichplot::gseaplot2(gseKEGG.res ,1:5, 
                      pvalue_table = FALSE
                      )

ggsave(file = file.path(outpath, paste0(cohort, '_gseKEGG.pdf')), 
       plot = gseKEGG_p, width = 10, height = 10)
write_delim(gseKEGG.res@result, 
            file = file.path(outpath, paste0(cohort, '_gseKEGG.txt')), delim = '\t')

# gsea go
gsea.go <- gseGO(gene_sort,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 pvalueCutoff = 1,
                 ont='BP',
                 pAdjustMethod = 'BH'
)
(gseGo_p <- enrichplot::gseaplot2(gsea.go, 1:5, 
                                  pvalue_table = T))
# gseaplot(gsea.go, geneSetID="GO:0006629")
ggsave(file = file.path(outpath, paste0(cohort, '_gseGO.pdf')), 
       plot = gseGo_p, width = 10, height = 10)
write_delim(gsea.go@result, 
            file = file.path(outpath, paste0(cohort, '_gseGO.txt')), delim = '\t')



