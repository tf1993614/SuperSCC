library(tidyverse)
library(ggrepel)

# read file
files = list.files(path = "/mnt/disk5/zhongmin/superscc/结果位置/免疫内皮细胞",
                   pattern = ".+csv",
                   full.names = T,
                   recursive = T)

dge_df = map(files, ~ read.csv(.x))[1:4]

# get the DE genes with customed cutoff
dge_df2 = map(
  dge_df,
  function(x) {
    x %>%
      subset(group == "target") %>%
      subset(pvals_adj < 0.05 & abs(logfoldchanges) > 1)
  }
)

# convert gene ID to symbol
id2symbol = AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = c("ALIAS", "ENSEMBL"))
first = function(x){x[[1]]}
dge_df2[[3]]$names = AnnotationDbi::mapIds(org.Hs.eg.db, keys = dge_df2[[3]]$names, column = c("ALIAS"), keytype = "ENSEMBL", multiVals = first)

# organize the dataframe for visuliaztion
final_df = dge_df[[2]] %>% subset(group == "target")  %>% subset(! is.na(logfoldchanges))
final_df$DE = final_df$pvals_adj < 0.05 & abs(final_df$logfoldchanges) > 1

# get the highlight text table
label_df = final_df %>%
  subset(str_detect(names, "^HLA.+") | str_detect(names, "CD74|IL3RA")) %>%
  subset(DE == T)

# set colors
colors = c("#E64B35FF", "#3C5488FF")
names(colors) = c(TRUE, FALSE)

# volcano plot
p1 = final_df  %>% 
  ggplot(aes(logfoldchanges, -log10(pvals_adj), colour = DE)) +
  geom_point() +
  coord_cartesian(xlim = c(-10, 10)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_label_repel(data = label_df, aes(logfoldchanges, -log10(pvals_adj), label = names)) +
  theme_bw() +
  scale_color_manual(values = colors) +
  theme(text = element_text(colour = "black"),
        axis.title.x.bottom = element_text(size =13, face = "bold"),
        axis.text.x.bottom = element_text(face = "bold",colour = "black"),
        axis.title.y.left = element_text(size =13, face = "bold"),
        axis.text.y.left = element_text(face = "bold", colour = "black"),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, vjust = 3,
                                  face = "bold", size = 15))

ggsave(plot = p1,
       filename = "volcano_plot_endothelial_Dong_2020_1.pdf",
       units = "px",
       width = 2500,
       height = 2000)
