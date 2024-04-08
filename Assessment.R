
#### LOADING LIBRARIES ####
library(dplyr)
library(ggplot2)
library(amap)
library(reshape2)
library(ggrepel)
library(cowplot)
library(clusterProfiler)
library(org.Hs.eg.db)


#### PARAMETERS ####

organism = org.Hs.eg.db

# significance
sig_cutoff = 0.001
fold_cutoff = 2
qvalue_cutoff = 0.01

# colours - significance
colour_sig = "#D6604D"
colour_sig_down = "#4393C3"
colour_nonsig = "grey18"
# colours - groups
RColorBrewer::brewer.pal(3,"Set1") # to get the code for colours from different palettes
colour_groups =  c("#E41A1C", "#377EB8", "#4DAF4A") # change?
colour_SvP = c("#E41A1C", "#377EB8")
colour_SMvS = c("#377EB8", "#4DAF4A")
colour_SMvP = c("#E41A1C", "#4DAF4A")

# colours heatmap
# colours
colours_hm = colorRampPalette(c("#2166AC", "mistyrose", "#B2182B"))(100)
colours_hm2 = colorRampPalette(c("#053061", "mistyrose", "#67001F"))(100)

# colours for pathway analysis barplots
colours_p = c("#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D")
colours_p = colorRampPalette(colours_p)(100)

# colours for fold vs fold scatterplot
colours_sig_groups = c("#E41A1C", "#377EB8", "#4DAF4A", "grey50")

# labels/levels
labels_groups = c("Proliferating", "Senescent", "Senescent_MtD")
labels_SvP = c("Proliferating", "Senescent")
labels_SMvS = c("Senescent", "Senescent_MtD")
labels_SMvP = c("Proliferating", "Senescent_MtD")
levels_groups = c("Prolif", "Senes", "Senes_MtD")
levels_SvP = c("Prolif", "Senes")
levels_SMvS = c("Senes", "Senes_MtD")
levels_SMvP = c("Prolif", "Senes_MtD")

#### INITIAL LOADING AND PARSING OF DATA ####

# Load data
annotations = read.table("Human_Background_GRCh38.p13.csv", header = TRUE, row.names = 1, sep = "\t")
em = em = read.table("EM.csv", header = TRUE, row.names = 1, sep = "\t") # 24858
ss = read.table("sample_sheet.csv", header = TRUE, sep = "\t")
de_SMvP = parse_de("DE_Senes_MtD_vs_Prolif.csv") # NA omitted: 24375 (from 24858)
de_SMvS = parse_de("DE_Senes_MtD_vs_Senes.csv") # NA omitted: 24857 (from 24858)
de_SvP = parse_de("DE_Senes_vs_Prolif.csv") # NA omitted: 24857 (from 24858)

# check whether observations that have no p.adj could be significant
nrow(subset(de_SMvP, is.na(p.adj) & p < 0.001 & abs(log2fold) > 2)) # 0
nrow(subset(de_SMvS, is.na(p.adj) & p < 0.001 & abs(log2fold) > 2)) # 0
nrow(subset(de_SvP, is.na(p.adj) & p < 0.001 & abs(log2fold) > 2)) # 0
# can omit NAs
de_SMvP = na.omit(de_SMvP)
de_SMvS = na.omit(de_SMvS)
de_SvP = na.omit(de_SvP)

# Annotate em 
em_annotated = annotate_em(em, annotations) 
nrow(em) # 24858
nrow(em_annotated) # 24600
# 258 observed genes have no annotations - probably not really important genes (pseudogenes/no known function etc.)

#### DELETE?? ####
# genes missing in annotations
em_annotated_all = merge(em, annotations, by.x = 0, by.y = 0, all.x = TRUE)
genes_missing = subset(em_annotated_all, is.na(SYMBOL))[,1] # 258
nrow(subset(de_SMvP[genes_missing,], p.adj < 0.001 & abs(log2fold) > 2)) # 16
nrow(subset(de_SMvS[genes_missing,], p.adj < 0.001 & abs(log2fold) > 2)) # 13
nrow(subset(de_SvP[genes_missing,], p.adj < 0.001 & abs(log2fold) > 2)) # 14
# relatively small part of that would be significant and probably not important genes anyway..

# Merge and parse de files
master_SMvP = merge_de(de_SMvP, em_annotated) # 24124 
master_SMvS = merge_de(de_SMvS, em_annotated) # 24599
master_SvP = merge_de(de_SvP, em_annotated) # 24599

# set levels for groups
ss$SAMPLE_GROUP = factor(ss$SAMPLE_GROUP, levels = levels_groups)

# create em with symbols for possible future use
em_symbols = em_annotated[,ss$SAMPLE]
row.names(em_symbols) = em_annotated[,"gene_name"]

# scale em
em.s = scale_em(em_annotated[,ss$SAMPLE])

# get significant genes (table and vector with names of significant genes)
master_SMvP_sig = filter_sig(master_SMvP)
sig_genes_SMvP = master_SMvP_sig$ENSEMBL_ID # 1798
master_SMvS_sig = filter_sig(master_SMvS) 
sig_genes_SMvS = master_SMvS_sig$ENSEMBL_ID # 1571
master_SvP_sig = filter_sig(master_SvP)
sig_genes_SvP = master_SvP_sig$ENSEMBL_ID # 1291

# get up- and down-regulated genes from master tables
sig_genes_SMvP_up = subset(master_SMvP, change == "Up")$ENSEMBL_ID # 907
sig_genes_SMvS_up = subset(master_SMvS, change == "Up")$ENSEMBL_ID # 891
sig_genes_SvP_up = subset(master_SvP, change == "Up")$ENSEMBL_ID # 466
sig_genes_SvP_down = subset(master_SvP, change == "Down")$ENSEMBL_ID # 429
sig_genes_SMvP_down = subset(master_SMvP, change == "Down")$ENSEMBL_ID # 1106
sig_genes_SMvS_down = subset(master_SMvS, change == "Down")$ENSEMBL_ID # 862


#### PLOTS ####

# Summary of significant genes in each category ####
summary_sig = data.frame(group = c("Senes vs Prolif", "Senes vs Prolif", "Senes_MtD vs Senes", "Senes_MtD vs Senes", "Senes_MtD vs Prolif", "Senes_MtD vs Prolif"),
                         change = c("up", "down", "up", "down", "up", "down"),
                         count = c(length(sig_genes_SvP_up), length(sig_genes_SvP_down), length(sig_genes_SMvS_up), length(sig_genes_SMvS_down), length(sig_genes_SMvP_up), length(sig_genes_SMvP_down)))
summary_sig$group = factor(summary_sig$group, levels = c("Senes vs Prolif", "Senes_MtD vs Senes", "Senes_MtD vs Prolif"))
ggp_summary = ggplot(summary_sig, aes(x = group, 
                       y = count, 
                       fill = change)) + 
  geom_col(position = position_dodge()) + 
  geom_label(aes(label = count, group = change), position = position_dodge(width = 0.9), fill = "white", size = 6) +
  labs(title = "Summary of significant Genes",
       y = "Number of Genes") +
  scale_fill_manual(name = "",
                    values = c(colour_sig, colour_sig_down),
                    labels = c("Up-regulated", "Down-regualted")) +
  theme_bw() +
  my_theme + theme(axis.title.x = element_blank())
ggsave("summary.svg", height = 5, width = 6)


# Principal Component Analysis ####
ggp_pca = pca_plot(em.s, scale = FALSE, samples_sheet = ss, colours = colour_groups, labels = labels_groups)
ggsave("pca.svg", height = 7, width = 8)

# Volcano plots ####
# replacing p.adjust = 0  (to avoid infinite values for mlog10p)
master_SvP.a = replace_inf(master_SvP, delete = TRUE)
master_SMvP.a = replace_inf(master_SMvP, delete = TRUE)
master_SMvS.a = replace_inf(master_SMvS, delete = TRUE)

# senes vs prolif 
ggp_volcano1 = volcano_plot(master_SvP.a, n_labeled = 10, title = "Senescent vs Proliferating Cells")
ggsave("volcano1.svg", height = 8, width = 10)

# senes_mtd vs prolif
ggp_volcano2 = volcano_plot(master_SMvP.a, n_labeled = 10, title = "Mitochondria Depleted Senescent vs Proliferating Cells")
ggsave("volcano2.svg", height = 8, width = 10)

# senes_mtd vs senes
ggp_volcano3 = volcano_plot(master_SMvS.a, n_labeled = 10, title = "Mitochondria Depleted Senescent vs Senescent Cells")
ggsave("volcano3.svg", height = 8, width = 10)

# MA plots ####
## senes vs prolif
ggp_ma1 = ma_plot(master_SvP, title = "Senescent vs Proliferating Cells", add_lines = TRUE)
ggsave("ma1.svg", height = 6, width = 7)
## senes_mtd vs prolif
ggp_ma2 = ma_plot(master_SMvP, title = "Mitochondria Depleted Senescent\nvs Proliferating Cells", add_lines = TRUE)
ggsave("ma2.svg", height = 6, width = 7)
## senes_mtd vs senes
ggp_ma3 = ma_plot(master_SMvS, title = "Mitochondria Depleted Senescent\nvs Senescent Cells",add_lines = TRUE)
ggsave("ma3.svg", height = 6, width = 7)

# Expression density per sample ####
ggp_exp = expression_density_plot(data = em, sample_sheet = ss, scale = FALSE, levels = levels_groups,
                                  colours = colour_groups, labels = labels_groups)
ggsave("exp_density.svg", height = 6, width = 7)
# scaled expression
ggp_exp.s = expression_density_plot(data = em, sample_sheet = ss, scale = TRUE, levels = levels_groups,
                                  colours = colour_groups, labels = labels_groups)
ggsave("exp_density_scaled.svg", height = 6, width = 7)


# Multiple gene violin box plots ####
ggp_violin_boxplot1 = multigene_boxplot_facet(get_em_subset(master_SvP_sig)[,-(7:9)], title = "Senescent vs Proliferating Cells",
                                              levels = levels_SvP, labels = labels_SvP, colours = colour_SvP)
ggsave("violin_boxplot1.svg", height = 7, width = 8.5)

ggp_violin_boxplot2 = multigene_boxplot_facet(get_em_subset(master_SMvP_sig)[,-(4:6)], title = "Mitochondria Depleted Senescent\nvs Proliferating Cells",
                                              levels = levels_SMvP, labels = labels_SMvP, colours = colour_SMvP)
ggsave("violin_boxplot2.svg", height = 7, width = 8.5)

ggp_violin_boxplot3 = multigene_boxplot_facet(get_em_subset(master_SMvS_sig)[,-(1:3)], title = "Mitochondria Depleted Senescent\nvs Senescent Cells",
                                              levels = levels_SMvS, labels = labels_SMvS, colours = colour_SMvS)
ggsave("violin_boxplot3.svg", height = 7, width = 8.5)


#### PATHWAY ANALYSIS ####
# SvP
ora_results_SvP = ora(sig_genes_SvP) # 13.79%
ora_results_SvP_up = ora(sig_genes_SvP_up) # 16.71%
ora_results_SvP_down = ora(sig_genes_SvP_down) # 7.93%
# SMvS
ora_results_SMvS = ora(sig_genes_SMvS) # 15.27%
ora_results_SMvS_up = ora(sig_genes_SMvS_up) # 15.02%
ora_results_SMvS_down = ora(sig_genes_SMvS_down) # 15.37%
# SMvP
ora_results_SMvP = ora(sig_genes_SMvP) # 14.4%
ora_results_SMvP_up = ora(sig_genes_SMvP_up) # 18.96%
ora_results_SMvP_down = ora(sig_genes_SMvP_down) # 9.76%

ora_table_SvP = get_ora_results(ora_results_SvP)
ora_table_SvP_up = get_ora_results(ora_results_SvP_up) 
ora_table_SvP_down = get_ora_results(ora_results_SvP_down) 
# SMvS
ora_table_SMvS = get_ora_results(ora_results_SMvS) 
ora_table_SMvS_up = get_ora_results(ora_results_SMvS_up) 
ora_table_SMvS_down = get_ora_results(ora_results_SMvS_down) 
# SMvP
ora_table_SMvP = get_ora_results(ora_results_SMvP)
ora_table_SMvP_up = get_ora_results(ora_results_SMvP_up)
ora_table_SMvP_down = get_ora_results(ora_results_SMvP_down)

ora_summary_all = rbind(get_summary_ora(ora_table_SvP, group = "SvP"),
                      get_summary_ora(ora_table_SvP_up, group = "SvP_up"),
                      get_summary_ora(ora_table_SvP_down, group = "SvP_down"),
                      get_summary_ora(ora_table_SMvS, group = "SMvS"),
                      get_summary_ora(ora_table_SMvS_up, group = "SMvS_up"),
                      get_summary_ora(ora_table_SMvS_down, group = "SMvS_down"),
                      get_summary_ora(ora_table_SMvP, group = "SMvP"),
                      get_summary_ora(ora_table_SMvP_up, group = "SMvP_up"),
                      get_summary_ora(ora_table_SMvP_down, group = "SMvP_down"))
ora_summary_all$group = factor(ora_summary_all$group, levels = unique(ora_summary_all$group))

# Bar plot to compare most significant ontologies between groups ####
barplot_multi = plot_multi_ora_barplots(rbind(ora_summary_all[grep("SvP", ora_summary_all$group),], ora_summary_all[grep("SMvS", ora_summary_all$group),]))
ggsave("barplot_multi.svg", height = 12, width = 20)

# get candidate genes from top ontology ####
candidate_genes_SvP = extract_enriched_genes(ora_table_SvP)
em_candidate_SvP = scale_em(na.omit(em_symbols[candidate_genes_SvP,]))
em_candidate_SvP$gene = row.names(em_candidate_SvP)
# SMvS
candidate_genes_SMvS = extract_enriched_genes(ora_table_SMvS)
em_candidate_SMvS = scale_em(na.omit(em_symbols[candidate_genes_SvP,]))
em_candidate_SMvS$gene = row.names(em_candidate_SvP)

# plot gene expression for candidate genes ####
ggp_multigene_go_SvP = multigene_boxplot(em_candidate_SvP[,-(7:9)], labels = labels_SvP, colours = colour_SvP, title = "Nuclear Division")
ggsave("multigene_SvP.svg", height = 4, width = 15)
ggp_multigene_go_SMvS = multigene_boxplot(em_candidate_SMvS[,-(1:3)], labels = labels_SMvS, colours = colour_SMvS, title = "Extracellular Matrix Organization")
ggsave("multigene_SMvS.svg", height = 4, width = 15)

#### COMPARING ALL THREE GROUPS ####

# Combining de tables and em_annotated
# only focusing on de_SvP and de_SMvS
de_all = combine_de(de1 = de_SvP, suffix_de1 = ".SvP",
                    de2 = de_SMvS, suffix_de2 = ".SMvS")
master = merge(em_annotated, de_all, by.x = 0, by.y = 0)
row.names(master) = master[,1]
master = master[,-1]

# getting all significant genes to use in heatmap
master_sig = subset(master,( sig.SMvS == TRUE | sig.SvP == TRUE))
sig_genes = row.names(master_sig)

# Calculating number of genes in each group - Venn diagram ####
n_sig_SvP = nrow(subset(master, sig.SvP == TRUE)) # 1291
n_sig_SvP_up = nrow(subset(master, (sig.SvP == TRUE & log2fold.SvP > 0))) # 862
n_sig_SvP_down = nrow(subset(master, (sig.SvP == TRUE & log2fold.SvP < 0))) # 429
n_sig_SMvS = nrow(subset(master, sig.SMvS == TRUE)) # 1572
n_sig_SMvS_up = nrow(subset(master, (sig.SMvS == TRUE & log2fold.SMvS > 0))) # 466
n_sig_SMvS_down = nrow(subset(master, (sig.SMvS == TRUE & log2fold.SMvS < 0))) # 1106
n_sig_all = nrow(master_sig) # 2520
n_sig_both_up = nrow(subset(master,(sig.SvP == TRUE & log2fold.SvP > 0 & sig.SMvS == TRUE & log2fold.SMvS > 0)))
n_sig_both_down = nrow(subset(master,(sig.SvP == TRUE & log2fold.SvP < 0 & sig.SMvS == TRUE & log2fold.SMvS < 0)))
n_sig_SvP_up_SMvS_down = nrow(subset(master,(sig.SvP == TRUE & log2fold.SvP > 0 & sig.SMvS == TRUE & log2fold.SMvS < 0)))
n_sig_SvP_down_SMvS_up = nrow(subset(master,(sig.SvP == TRUE & log2fold.SvP < 0 & sig.SMvS == TRUE & log2fold.SMvS > 0)))


# Fold vs Fold scatterplot ####
ggp_foldVfold = plot_foldVfold(master,
               x = log2fold.SvP, suffix_de1 = ".SvP", label_de1 = "Pro vs Sen",
               y = log2fold.SMvS, suffix_de2 = ".SMvS", label_de2 = "Sen vs Sen_MtD")
ggsave("FoldvFold.svg", height = 7, width = 7)

# Heat map ####
ggp_heatmap_all = heat_map(em[sig_genes,], title = "All Significant Genes")
ggsave("heatmap_all.svg", height = 10, width = 7)

# Signatures ####
# up/up - 349
signature_1 = row.names(subset(master_sig, (master_sig$log2fold.SvP > 0 & master_sig$log2fold.SMvS > 0 )))
# down/down - 296
signature_2 = row.names(subset(master_sig, (master_sig$log2fold.SvP < 0 & master_sig$log2fold.SMvS < 0 )))
# up/down - 1346
signature_3 = row.names(subset(master_sig, (master_sig$log2fold.SvP > 0 & master_sig$log2fold.SMvS < 0 )))
# down/up - 524
signature_4 = row.names(subset(master_sig, (master_sig$log2fold.SvP < 0 & master_sig$log2fold.SMvS > 0 )))

ggp_heatmap_sign1 = heat_map(em[signature_1,])
ggp_heatmap_sign2 = heat_map(em[signature_2,])
ggp_heatmap_sign3 = heat_map(em[signature_3,])
ggp_heatmap_sign4 = heat_map(em[signature_4,])

# Metagene box plot ####
ggp_metagene1 = plot_metagene(signature_1, title = "Signature 1")
ggsave("metagene1.svg", height = 6, width = 4)
ggp_metagene2 = plot_metagene(signature_2, title = "Signature 2")
ggsave("metagene2.svg", height = 6, width = 4)
ggp_metagene3 = plot_metagene(signature_3, title = "Signature 3")
ggsave("metagene3.svg", height = 6, width = 4)
ggp_metagene4 = plot_metagene(signature_4, title = "Signature 4")
ggsave("metagene4.svg", height = 6, width = 4)

# Pathway analysis for signatures ####
ora_results_sign1 = ora(signature_1, "ENSEMBL", organism_db = org.Hs.eg.db) # 15.47% gene IDs failed to map
ora_results_sign2 = ora(signature_2, "ENSEMBL", organism_db = org.Hs.eg.db) # 9.8%
ora_results_sign3 = ora(signature_3, "ENSEMBL", organism_db = org.Hs.eg.db) # 17.38%
ora_results_sign4 = ora(signature_4, "ENSEMBL", organism_db = org.Hs.eg.db) # 11.45%

# get results from the pathway analysis
ora_table_sign1 = get_ora_results(ora_results_sign1)
ora_table_sign2 = get_ora_results(ora_results_sign2)
ora_table_sign3 = get_ora_results(ora_results_sign3)
ora_table_sign4 = get_ora_results(ora_results_sign4)

# bind all summaries from ora into one table to be able to plot together
ora_summary_sign_all = rbind(get_summary_ora(ora_table_sign1, "Signature 1"),
                             get_summary_ora(ora_table_sign2, "Signature 2"),
                             get_summary_ora(ora_table_sign3, "Signature 3"),
                             get_summary_ora(ora_table_sign4, "Signature 4"))

# plot top 10 ontologies from signatures ####
barplot_multi_sign = plot_multi_ora_barplots(ora_summary_sign_all)
ggsave("barplot_multi_sign.svg", heigh = 14, width = 20)

# extract candidate genes - as symbols (from ORA) ####
sign1_candidate_genes = extract_enriched_genes(ora_table_sign1, gene_set = "response to hypoxia")
sign2_candidate_genes = extract_enriched_genes(ora_table_sign2, gene_set = "nuclear division")
sign3_candidate_genes = extract_enriched_genes(ora_table_sign3, gene_set = "positive regulation of cytokine production")
sign4_candidate_genes = extract_enriched_genes(ora_table_sign4, gene_set = "organelle fission")

em_sign1 = scale_em(na.omit(em_symbols[sign1_candidate_genes,]))
em_sign1$gene = row.names(em_sign1)

em_sign2 = scale_em(na.omit(em_symbols[sign2_candidate_genes,]))
em_sign2$gene = row.names(em_sign2)

em_sign3 = scale_em(na.omit(em_symbols[sign3_candidate_genes,]))
em_sign3$gene = row.names(em_sign3)

em_sign4 = scale_em(na.omit(em_symbols[sign4_candidate_genes,]))
em_sign4$gene = row.names(em_sign4)

# plot multi gene barplot from top ontology ####
ggp_multigenesign1 = multigene_boxplot(em_sign1, title = "Response to Hypoxia")
ggsave("multigene1.svg", height = 4, width = 16)
ggp_multigenesign2 = multigene_boxplot(em_sign2, title = "Nuclear Division")
ggsave("multigene2.svg", height = 4, width = 16)
ggp_multigenesign3 = multigene_boxplot(em_sign3, title = "Positive Regulation of Cytokine Production")
ggsave("multigene3.svg", height = 4, width = 16)
ggp_multigenesign4 = multigene_boxplot(em_sign4, title = "Organelle Fission")
ggsave("multigene4.svg", height = 4, width = 16)


# create multigene box plots from top 10 candidate genes ####
ggp_sign1_boxplot = violin_boxplot_genes(em_sign1[1:10,])
ggsave("sign1_violinbox.svg", height = 6, width = 8)

ggp_sign2_boxplot = violin_boxplot_genes(em_sign2[1:10,])
ggsave("sign2_violinbox.svg", height = 6, width = 8)

ggp_sign3_boxplot = violin_boxplot_genes(em_sign3[1:10,])
ggsave("sign3_violinbox.svg", height = 6, width = 8)

ggp_sign4_boxplot = violin_boxplot_genes(em_sign4[1:10,])
ggsave("sign4_violinbox.svg", height = 6, width = 8)

