library(tidyverse)
library(ggrepel)

setwd("~/jupyter_notebooks/working_script/gene_module/F_cluster")

# draw stack plot to show the frequency of stress B cell
# for clinical meta feature
meta = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/应激_免疫内皮metadata/Laughney_2020_原始metadata.csv")
meta = meta %>% 
  dplyr::select(treated_naive, 
                AJCC_stage, 
                sample_primary_met,
                patient,
                Fcluster) %>%
  subset(Fcluster %in% c(6, 7, 10)) %>%
  mutate(stress = if_else(Fcluster == 6, "yes", "no")) 

sub_meta1 = meta %>% group_by(stress) %>% count(treated_naive) %>% mutate(type = "Treatment") 
sub_meta2 = meta %>% group_by(stress) %>% count(AJCC_stage) %>% mutate(type = "AJCC_stage") %>% dplyr::rename(meta = "AJCC_stage")
sub_meta3 = meta %>% group_by(stress) %>% count(sample_primary_met) %>% mutate(type = "Sample_primary_meta") %>% dplyr::rename(meta = "sample_primary_met")
sub_meta4 = meta %>% group_by(stress) %>% count(patient) %>% mutate(type = "Patient")%>% dplyr::rename(meta = "patient")
sub_others = list(sub_meta2, sub_meta3, sub_meta4)

colors = c("#E64B35FF", "#3C5488FF")
names(colors) = c("yes", "no")

p2 = sub_meta1 %>%
  ggplot(aes(treated_naive, n, fill = stress)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = colors, labels = ~ if_else(.x == "yes", "Stress B", "Normal B")) +
  theme_bw() +
  ylab("Frequency") +
  facet_wrap(~type) +
  theme(text = element_text(colour = "black"),
        axis.title.x.bottom = element_blank(),
        axis.text.x.bottom = element_text(colour = "black"),
        axis.title.y.left = element_text(size =13),
        axis.text.y.left = element_text(colour = "black"),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.title = element_blank(),
        strip.background = element_rect(fill = "black", color = "black"),
        strip.text = element_text(color = "white")
        )

ggsave(p2, filename = "Treatment_stress_B.pdf", width = 2000, height = 1500, units = "px")

p3 = imap(
  sub_others,
   function(x, y) {
     fig = x %>%
    ggplot(aes(meta, n, fill = stress)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = colors, labels = ~ if_else(.x == "yes", "Stress B", "Normal B")) +
    theme_bw() +
    ylab("Frequency") +
    facet_wrap(~type) +
    theme(text = element_text(colour = "black"),
          axis.title.x.bottom = element_blank(),
          axis.text.x.bottom = element_text(colour = "black", angle = 45, hjust = 1, vjust = 1),
          axis.title.y.left = element_text(size =13),
          axis.text.y.left = element_text(colour = "black"),
          panel.border = element_rect(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          legend.title = element_blank(),
          strip.background = element_rect(fill = "black", color = "black"),
          strip.text = element_text(color = "white"),
    )
   }
) 

p3[[1]] = p3[[1]] + theme(legend.position = "none")
p3[[2]] = p3[[2]] + theme(legend.position = "none", axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank())
p3[[3]] = p3[[3]] + theme(axis.title.y.left = element_blank(), axis.text.y.left = element_blank(), axis.ticks.y.left = element_blank())

p3 = p3[[1]] + p3[[2]] + p3[[3]]

ggsave(p3, filename = "Other_clinical_features_stress_B.pdf", width = 3000, height = 1500, units = "px")


# draw volcano plot to show the DE genes 
# between "Stressed B-cell-associated T cells" and 
# "other T cells"
volcano_df = read.csv("/mnt/disk5/zhongmin/superscc/结果位置/Data_Laughney2020_Lung/3个亚群的marker.csv")

volcano_df1  = volcano_df %>% subset(cluster == "Stressed B-cell-associated T cells") %>%
  subset(avg_log2FC > 0)

volcano_df2  = volcano_df %>% subset(cluster == "Stressed B-cell-associated T cells") %>%
  subset(avg_log2FC < 0)

volcano_df3 = bind_rows(volcano_df1, volcano_df2)
  
volcano_df3 = volcano_df3 %>% 
  mutate(
    DE = case_when(
      avg_log2FC > 0.5 & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -0.5 & p_val_adj < 0.05 ~ "Down",
      .default = "No difference"
    )
  )

label_df = volcano_df3 %>% 
  subset(gene %in% c(
    "HSPA1A", "HSP90AA1", "HSP90AB2",
    "HSPD1", "HSPA1B", "HSPE1", "HSPH1",
    "HSPA8", "HSPB1", "HSPA6", "HSPA4",
    "HSPA2", "HSPA1L"
  ))

colors = c("#E64B35FF", "#3C5488FF", "lightgrey")
names(colors) = c("Up", "Down", "No difference")

p4 = volcano_df3  %>% 
  ggplot(aes(avg_log2FC, -log10(p_val_adj), colour = DE)) +
  geom_point() +
  coord_cartesian(xlim = c(-2.5, 5)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_label_repel(data = label_df, aes(avg_log2FC, -log10(p_val_adj), label = gene)) + 
  theme_bw() +
  scale_color_manual(values = colors) +
  theme(text = element_text(colour = "black"),
        axis.title.x.bottom = element_text(size =13),
        axis.text.x.bottom = element_text(,colour = "black"),
        axis.title.y.left = element_text(size =13),
        axis.text.y.left = element_text(colour = "black"),
        panel.border = element_rect(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 3,
                                  size = 15)) 

ggsave(p4, filename = "DE_Stressed B-cell-associated T cells.pdf", width = 2000, height = 1500, units = "px")

