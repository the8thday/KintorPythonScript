#! /urs/bin/env Rscript
# Title     :
# Objective : cluster analysis
# Created by: congliu
# Created on: 2021/8/17

library(tidyverse)

tfre <- readxl::read_excel('/Users/congliu/OneDrive/kintor/qianx/Kintor-myc_Protac-质谱检测结果-20210812.xlsx',
                           sheet='TFRE')

profiling <- readxl::read_excel('/Users/congliu/OneDrive/kintor/qianx/Kintor-myc_Protac-质谱检测结果-20210812.xlsx',
                                sheet='profiling')

group_info <- readxl::read_excel('/Users/congliu/OneDrive/kintor/qianx/Kintor-myc_Protac-质谱检测结果-20210812.xlsx',
                                sheet='sample_info')

# profiling <- readxl::read_excel('/Users/congliu/OneDrive/kintor/qianx/MYC_profiling.xlsx')
tfre <- profiling
tfre_df <- tfre %>%
  select(process_code, gene_id, gene_symbol, ifot) %>%
  pivot_wider(names_from = process_code, values_from = ifot)

# get id, 注意修改不同的sheet
cohert <- 'TFRE_3nm'
path <- '/Users/congliu/OneDrive/kintor/qianx/drop_na'
group_id <- group_info %>% filter(methods=='profiling') %>%
  filter(cellProcess %in% c('ctrl','3nM')) %>% pull(processcode)

tfre_1 <- tfre_df %>%
  select(gene_id,gene_symbol,{{ group_id }}) %>%
  drop_na()
# tfre_1 <- tfre_1 %>% replace(is.na(.),0)
# tfre_1 <- tfre_1[!rowSums(tfre_1==0)>3,]
# tfre_1 <- tfre_1 %>% select(where(~sum(is.na(.x))>3))

Pvalue <- rep(0, nrow(tfre_1))
log2_FC <- rep(0, nrow(tfre_1))
mean_test <- rep(0, nrow(tfre_1))
mean_control <- rep(0, nrow(tfre_1))

for(i in seq_len(nrow(tfre_1))){
        if(sd(tfre_1[i,3:5])==0 && sd(tfre_1[i,6:8])==0){
                Pvalue[i] <- "NA"
                log2_FC[i]<- "NA"
        }else{
                y <- t.test(as.numeric(tfre_1[i, 3:5]), as.numeric(tfre_1[i, 6:8]))
                Pvalue[i] <- y$p.value
                log2_FC[i] <- log2((mean(as.numeric(tfre_1[i,6:8])))/(mean(as.numeric(tfre_1[i,3:5]))))
                mean_control[i] <- y$estimate[[1]]
                mean_test[i] <- y$estimate[[2]]
        }
}

# 对p value进行FDR校正
fdr <- p.adjust(Pvalue, "BH")
# 在原文件后面加入log2FC，p value和FDR,共3列；
out <- cbind(tfre_1,mean_control, mean_test, log2_FC,Pvalue,fdr) %>%
  arrange(fdr, log2_FC)

write_delim(out, file.path(path, glue::glue('{cohert}_control.txt')),
            delim = '\t')
out <- as_tibble(out) %>%
  mutate(log2_FC = as.numeric(log2_FC))
# plot: volcano+MAplot
degs <- out %>% filter(fdr <= 0.05) %>%
  filter(abs(log2_FC)>=1) %>%
  select(gene_symbol, mean_control, mean_test,log2_FC,Pvalue,fdr)

write_delim(degs, file.path(path, glue::glue('{cohert}_degs.txt')),
            delim = '\t')

df_p <- out %>% mutate(baseMean = log2(sqrt(mean_control * mean_test))) %>%
  mutate(label = case_when(
  log2_FC >= 1 & fdr <= 0.05 ~ 'up',
  log2_FC < -1 & fdr <= 0.05 ~ 'down',
  TRUE ~ 'noSig'
)) %>% arrange(fdr)

top10_up <- df_p %>% filter(label == 'up') %>%
  slice_head(n=10)
top10_down <- df_p %>% filter(label == 'down') %>%
  slice_head(n=10)

p <- ggplot(df_p, aes(x = log2_FC, y = -log10(fdr), color = label)) +
  geom_point(size = 2) +
  scale_colour_manual(values  = c("#B31B21", "#1465AC", "darkgray"), limits = c('up', 'down', 'noSig')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'),
        legend.position = c(0.9, 0.93)) +
  geom_vline(xintercept = c(-1, 1), color = 'gray', size = 0.3) +
  geom_hline(yintercept = -log(0.01, 10), color = 'gray', size = 0.3) +
  # xlim(-5, 5) +
  # ylim(0, 6) +
  labs(x = '\nLog2 Fold Change', y = '-log10(FDR)', color = '', title = 'Treat vs Control\n') +
  ggrepel::geom_text_repel(data = bind_rows(top10_up, top10_down), aes(x = log2_FC, y = -log10(fdr),
                                                                       label = gene_symbol),
        size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)
ggsave(plot = p, filename = file.path(path, glue::glue('{cohert}_volcano.pdf')),
       width = 10, height = 10, dpi = 500)

# MAplot
p2 <- ggplot(df_p, aes(x = baseMean, y = log2_FC, color = label)) +
  geom_point(size = 2) +
  scale_colour_manual(values  = c("#B31B21", "#1465AC", "darkgray"), limits = c('up', 'down', 'noSig')) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2),
      panel.background = element_rect(color = 'black', fill = 'transparent'),
      legend.key = element_rect(fill = 'transparent'), legend.position = 'top') +
  labs(x = '\nLog2 Mean Expression', y = 'Log2 Fold Change\n', color = '') +
  geom_hline(yintercept = c(0, -log2(2), log2(2)), linetype = c(1, 2, 2),
               color = c("black", "black", "black")) +
  ggrepel::geom_text_repel(data = bind_rows(top10_up, top10_down), aes(x = baseMean, y = log2_FC,
                                                                       label = gene_symbol),
        size = 3,box.padding = unit(0.5, 'lines'), segment.color = 'black', show.legend = FALSE)

ggsave(plot = p2, filename = file.path(path, glue::glue('{cohert}_MAplot.pdf')),
       width = 10, height = 10, dpi = 500)


# intresting gene
i_gene <- Hmisc::Cs(CRL4,IKZF1,IKZF3,GSPT1,ZFP91,ZNF98,FBW7)

tfre_1 %>% filter(gene_symbol %in% i_gene) %>%
  write_delim(file.path(path, glue::glue('{cohert}_Target.txt')),
            delim = '\t')
