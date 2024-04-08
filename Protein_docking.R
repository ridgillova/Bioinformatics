library(ggplot2)
library(tidyverse)
library(ggrepel)

docking = read.table("protein_docking.csv", header = TRUE, sep = ",")
docking = docking[-118,1:3]
docking$Score = as.numeric(docking$Score)
docking$Label = case_when((docking$Ligand == "ZINC00646292" & docking$Pose == 3) ~ "ZINC00646292",
                         (docking$Ligand == "ZINC00833965" & docking$Pose == 3) ~ "ZINC00833965",
                         (docking$Ligand == "ZINC00547875" & docking$Pose == 1) ~ "ZINC00547875",
                         (docking$Ligand == "Zanamivir" & docking$Pose == 2) ~ "Zanamivir")
ggplot(docking, aes(x = Score)) + 
  geom_density(fill = "#E41A1C", alpha = 0.5) +
  labs(title = "Density Distribution of Energy Scores", 
       y = "Density") +
  geom_vline(xintercept = -7, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = -7.4, colour = "#4DAF4A") +
  geom_vline(xintercept = -8.1, colour = "#FF7F00") +
  geom_vline(xintercept = -6.1, colour = "#377EB8") +
  geom_vline(xintercept = -7.1, colour = "#984EA3") +
  geom_text_repel(aes(label = Label, colour = Label), y = 0.78, nudge_x = 0.2, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("#4DAF4A","#FF7F00","#377EB8","#984EA3")) +
  xlim(-8.1,-5.5) +
  theme_bw() + 
  theme(plot.title = element_text(size = 15, face = "bold"), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10))
ggsave("density_plot3.svg", width = 6.5, height = 4)
