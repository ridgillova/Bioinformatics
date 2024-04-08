
#### LOAD AND PARSE DATA ####

# Annotate EM ####
annotate_em = function(em, annotation)
{
  names(annotations) = c("gene_name", "chromosome", "start", "stop", "gene_type")
  em_annotated = merge(em,annotations, by.x = 0, by.y = 0)
#  row.names(em_annotated) = em_annotated[,"gene_name"]
#  names(em_annotated)[1] = "ENSEMBL_ID"
  row.names(em_annotated) = em_annotated[,1]
  em_annotated = em_annotated[,-1]
  return(em_annotated)
}

# Parsing DE ####

# categorizing as up/down/non-significant ####
categorise_change = function(data, sig_threshold = sig_cutoff, fold_threshold = fold_cutoff)
{
  data$change = case_when((data$p.adj < sig_threshold & data$log2fold < (-fold_threshold)) ~ "Down",
                          (data$p.adj < sig_threshold & data$log2fold > fold_threshold) ~ "Up",
                          TRUE ~ "No change")
  data$change = factor(data$change, levels = c("Up", "Down", "No change"))
  return(data)
}

# categorizing as sig/non-significant ####
categorise_sig = function(data, sig_threshold = sig_cutoff, fold_threshold = fold_cutoff)
{
  data$sig = ifelse(data$p.adj < sig_threshold & abs(data$log2fold) > fold_threshold, TRUE, FALSE)
  data$sig = factor(data$sig, levels = c("TRUE", "FALSE"))
  return(data)
}


parse_de = function(de_file, samples_sheet = ss, sig_threshold = sig_cutoff, fold_threshold = fold_cutoff)
{
  de = read.table(de_file, header = TRUE, row.names = 1, sep = "\t")
#  de = na.omit(de)
  de = categorise_sig(de, sig_threshold, fold_threshold)
  de = categorise_change(de, sig_threshold, fold_threshold)
  sorted_order = order(de[,"p.adj"], decreasing = FALSE)
  de = de[sorted_order,]
  return(de)
}


merge_de = function(de, em_annotated, samples_sheet = ss)
{
  master = merge(em_annotated, de, by.x = 0, by.y = 0)
  row.names(master) = master[,"gene_name"]
  names(master)[1] = "ENSEMBL_ID"
  master$mean_expression = rowMeans(master[,as.vector(samples_sheet$SAMPLE)])
  master$mlog10p = -log((master$p.adj),10)
  return(master)
}

replace_inf = function(master, delete = FALSE)
{
  if(delete == FALSE){
  min_p.adj = min(master$p.adj[master$p.adj>0]) # get the min p.adj that is greater than 0
  master$p.adj[master$p.adj==0] = min_p.adj # replace all p.adj = 0 by the min p.adj
  master$mlog10p = -log((master$p.adj),10)
  return(master)
  } else {
    master = subset(master, mlog10p != "Inf")
  }
  return(master)
}

# combine DE tables into one ####
combine_de = function(de1, de2, suffix_de1, suffix_de2)
{
  de_all = merge(de1,de2, by.x = 0, by.y = 0, suffixes = c(suffix_de1, suffix_de2))
  row.names(de_all) = de_all[,1]
  de_all = de_all[,-1]
  return(de_all)
}


# Filter only significant genes function ####
filter_sig = function(data, sig_threshold = sig_cutoff, fold_threshold = fold_cutoff)
{
  data_sig = subset(data, data$p.adj < sig_threshold & abs(data$log2fold) > fold_threshold)
  return(data_sig)
}

# Scale expression data ####
scale_em = function(em, samples_sheet = ss)
{
  em = em[,samples_sheet$SAMPLE]
  em.s = scale(data.frame(t(em)))
  em.s = data.frame(t(em.s))
  em.s = na.omit(em.s)
  return(em.s)
}


#### CUSTOM THEME ####

my_theme = theme(text = element_text(family = "sans"),
                 panel.grid.minor = element_blank(),
                 plot.title = element_text(face = "bold", size = 24),
                 axis.text.x = element_text(size = 14),
                 axis.text.y = element_text(size = 14),
                 axis.title.x = element_text(size = 20),
                 axis.title.y = element_text(size = 20),
                 legend.text = element_text(size = 16), 
                 legend.title = element_text(size = 20),
                 strip.text.x = element_text(size = 12, vjust=0.5), 
                 legend.position = "bottom")

#### PLOTS AND ANALYSIS ####

# Volcano plot ####
volcano_plot = function(data, sig_threshold = sig_cutoff, fold_threshold = fold_cutoff, colours = c("#E41A1C", "#377EB8", "black"), n_labeled = 5, title = "Volcano Plot", add_lines = TRUE)
{
  #data$change = case_when((data$p.adj < sig_threshold & data$log2fold < (-fold_threshold)) ~ "b",
  #                                    (data$p.adj < sig_threshold & data$log2fold > fold_threshold) ~ "a",
  #                                    TRUE ~ "c")
  # update master_sig to include change column
  master_sig = subset(data, sig == TRUE)
  # make sure master_sig is order by p.adj - probably not necessary (but can include)
  sorted_order = order(master_sig$p.adj, decreasing = FALSE)
  master_sig = master_sig[sorted_order,]
  # create tables to be used for labels
  master_sig_up = subset(master_sig, master_sig$log2fold > 0)
  master_sig_down = subset(master_sig, master_sig$log2fold < 0)
  master_up_top5 = subset(master_sig_up[1:n_labeled,])
  master_down_top5 = subset(master_sig_down[1:n_labeled,])
  master_up_down_top5 = rbind(master_down_top5, master_up_top5)
  # plot
  ggp = ggplot(data = data,
         aes(x = log2fold,
             y = mlog10p,
             colour = change)) +
    geom_point(alpha = 0.5, size = 2) + 
    geom_label_repel(data = master_up_down_top5, 
                     aes(label = gene_name,
                         colour = change),
                     size = 5,
                     show.legend = FALSE) +
    scale_colour_manual(values = colours,
                        labels = c("Up-regulated","Down-regulated", "Non-significant"),
                        name = "") +
    labs(title = title, 
         x = "Log2 Fold Change", 
         y = "-Log(adjusted p-value)") +
    theme_bw() +
    my_theme
  # add horizontal and vertical lines to indicate cutoff values if add_lines set to TRUE
  if(add_lines == TRUE)
  {
    ggp = ggp  + 
    geom_vline(xintercept = -fold_threshold,
               linetype = "dashed",
               colour = "darkgrey",
               linewidth = 0.5) +
    geom_vline(xintercept = fold_threshold,
               linetype = "dashed", 
               colour = "darkgrey", 
               linewidth = 0.5) +
    geom_hline(yintercept = -log(sig_threshold),
               linetype = "dashed", 
               colour = "darkgrey", 
               linewidth = 0.5)
  }
  return(ggp)
}

# MA plot ####
ma_plot = function(data, sig_threshold = sig_cutoff, fold_threshold = fold_cutoff, colours = c(colour_sig, colour_nonsig), title = "MA Plot", add_lines = FALSE)
{
  ggp = ggplot(data = data,
         aes(x = log10(mean_expression),
             y = log2fold)) +
  geom_point(aes(colour = sig)) +
  labs(title = title,
         x = "Log(Mean Expression)",
         y = "Log2 Fold Change") +
  scale_colour_manual(name = "",
                        values = colours,
                        labels = c("Significant", "Non-significant")) +
  theme_bw() +
  my_theme
  if (add_lines == TRUE)
  {
    ggp = ggp + geom_hline(yintercept = fold_threshold,
                           linetype = "dashed",
                           colour = "darkgrey",
                           linewidth = 0.5) +
      geom_hline(yintercept = -fold_threshold,
                 linetype = "dashed",
                 colour = "darkgrey",
                 linewidth = 0.5)
  }
  return(ggp)
}

# PCA plot ####
pca_plot = function(data, scale = TRUE, samples_sheet = ss, colours = colour_groups, labels = labels_groups, title = "PCA")
{
  if (scale == TRUE)
  {
    data = scale_em(data)
  }
  # pca
  pca = prcomp(t(as.matrix(sapply(data, as.numeric))))
  pca_coordinates = data.frame(pca$x)
  # get variance
  vars = apply(pca$x, 2, var)
  prop_x = round(vars["PC1"] / sum(vars),4) * 100
  prop_y = round(vars["PC2"] / sum(vars),4) * 100
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="")
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="")
  # plot
  ggp = ggplot(pca_coordinates, aes(x = PC1,
                                    y = PC2,
                                    colour = samples_sheet$SAMPLE_GROUP)) + 
    geom_point(size = 2) +
    geom_text_repel(aes(label = samples_sheet$SAMPLE),
                    size = 6,
                    force_pull = 2,
                    force = 2,
                    max.overlaps = 20,
#                    min.segment.length = 0,
#                    segment.color = NA,
                    show.legend = FALSE) + 
    scale_colour_manual(name = "",
                        values = colours,
                        labels = labels) +
    labs(title = title,
         x = x_axis_label,
         y = y_axis_label) + 
    theme_bw() +
    my_theme + 
    theme(legend.text = element_text(size = 14))
  return(ggp)
}

# Expression density plot ####
expression_density_plot = function(data, sample_sheet = ss, scale = FALSE, levels = levels_groups, colours = colour_groups, labels = labels_groups, title = "Distribution of Expression")
{
  if(scale == TRUE)
  {
    em.s = scale_em(data)
    em.m = reshape2::melt(em.s)
    
    em.m = merge(em.m, ss, by.x = 1, by.y = 1)
    names(em.m) = c("sample", "expression", "group")
    # change the order of levels
    em.m$Group = factor(em.m$group, levels = levels)
    # making the density plot
    ggp = ggplot(em.m, 
                 aes(x = expression, fill = group))
    x_lab = "Expression (z-scores)"
  }
  if(scale == FALSE)
  {
    em.m = reshape2::melt(data)
    
    em.m = merge(em.m, ss, by.x = 1, by.y = 1)
    names(em.m) = c("sample", "expression", "group")
    # change the order of levels
    em.m$Group = factor(em.m$group, levels = levels)
    # making the density plot
    ggp = ggplot(em.m, 
                 aes(x = log(expression+0.01), fill = group))
    x_lab = "Log(Expression)"
  }
  ggp = ggp + geom_density(alpha = 0.7) + 
    facet_wrap(~sample, ncol = 3) +
    scale_fill_manual(name = "",
                      values = colours,
                      labels = labels) +
    labs(title = title,
         x = x_lab,
         y = "Density") + 
    theme_bw() +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
    my_theme
  return(ggp)
}


# Multiple box plots ####

# Get subset of expression values to visualise ####
get_em_subset = function(data, samples_sheet = ss, n = 5) 
{
  data_sig = subset(data, sig == TRUE)
  # pick n most significant genes for both up-regulated and down-regulated
  data_up = subset(data_sig, (data_sig$log2fold) > 0)
  sorted_order_up = order(data_up$log2fold, decreasing = TRUE)
  data_up = data_up[sorted_order_up,]
  data_down = subset(data_sig, (data_sig$log2fold) < 0)
  sorted_order_down = order(data_down$log2fold, decreasing = FALSE)
  data_down = data_down[sorted_order_down,]
  data_subset = rbind(data_up[1:n,], data_down[1:n,])
  genes = row.names(data_subset)
  data_subset = data_subset[,samples_sheet$SAMPLE]
  data_subset = scale_em(data_subset)
  data_subset$gene = factor(genes, levels = genes) # to keep the most up-regulated and most down-regulated together and ordered according to log2fold change
  return(data_subset)
}

# One plot with multiple genes on x-axis ####
multigene_boxplot = function(em, scaled = TRUE, samples_sheet = ss, n = NULL, colours = colour_groups, labels = labels_groups, title = "Gene Expression")
{
  em.m = melt(em)
  em.m = merge(em.m, samples_sheet, by.x = 2, by.y = 1)
  names(em.m) = c("sample", "gene", "expression", "group")
  if (scaled == TRUE){
    ggp = ggplot(em.m, aes(x = gene,
                         y = expression, 
                         fill = group, 
                         colour = group)) + 
      geom_boxplot(alpha = 0.7)
    y_lab = "Expression\n(z-scores)"
  }
  if (scaled == FALSE){
    ggp = ggplot(em.m, aes(x = gene,
                         y = log(expression+0.01), 
                         fill = group, 
                         colour = group)) + 
      geom_boxplot(alpha = 0.7)
    y_lab = "Log(Expression)"
  }
  ggp = ggp + 
    scale_fill_manual(aesthetics = c("fill", "colour"),
                      name = "", 
                      values = colours,
                      labels = labels) +
    labs(title = title,
         x = "",
         y = y_lab) +
    theme_bw() + 
    my_theme + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(size = 13, angle = 90))
  return(ggp)
}

# Multiple genes in facets ####
multigene_boxplot_facet = function(em, scaled = TRUE, samples_sheet = ss, n = 5, points = FALSE, levels = levels_groups, labels = labels_groups, colours = colour_groups, title = "Gene Expression")
{
  # melt data
  gene_data = melt(em)
  # merge with samples_sheet to get group information
  gene_data = merge(gene_data, samples_sheet, by.x = 2, by.y = 1)
  names(gene_data) = c("sample", "gene", "expression", "group")
  gene_data$group = factor(gene_data$group,
                           levels = levels_groups)
  # plot
  if (scaled == TRUE)
  {
    ggp = ggplot(gene_data, 
                 aes(x = group, 
                     y = expression,
                     fill = group))
    y_lab = "Expression (z-scores)"
  }
  if (scaled == FALSE)
  {
    ggp = ggplot(gene_data, 
                 aes(x = group, 
                     y = log(expression+0.01),
                     fill = group))
    y_lab = "Log(Expression+0.01)"
  }
  ggp = ggp +
    geom_violin(aes(colour = group),alpha = 0.7, scale = "width", trim = FALSE) +
    geom_boxplot(width = 0.2)
  if (points == TRUE)
  {
    ggp = ggp + geom_jitter(aes(fill = group),colour = "black",pch = 21, width = 0.2)
      
  }
  ggp = ggp +
    facet_wrap(~gene, ncol = n) +
    scale_fill_manual(aesthetics = c("fill", "colour"),
                      name = "", 
                      values = colours,
                      labels = labels) +
    labs(title = title,
         x = "Sample Group",
         y = y_lab) +
    theme_bw() + 
    my_theme +
    theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  return(ggp)
}

# Heat map ####
heat_map = function(data, samples_sheet = ss, cluster_axis = "y", dist_method = "spearman", cluster_method = "average", reorder_method = "average", levels_rug = levels_groups, colours_heatmap = colours_hm, colours_rug = colour_groups, labels_rug = labels_groups, title = "Heat Map")
{
  data = data[,samples_sheet$SAMPLE] # ensure only expression values present
  data = na.omit(data) # ensure no NAs in data
  data = scale_em(data)
  heat_map.mx = as.matrix(data)
  # use grepl method to check whether x or y are present in cluster_axis
  # clustering - y
  if (any(grepl("y", cluster_axis)))
  {
    y_dist = amap::Dist(heat_map.mx, method = dist_method)
    y_cluster = hclust(y_dist, method = cluster_method)
    y_dendrogram = as.dendrogram(y_cluster)
    y_dendrogram_reorder = reorder(y_dendrogram, 0, FUN = reorder_method)
    y_order = order.dendrogram(y_dendrogram_reorder)
    heat_map.c = heat_map.mx[y_order,]
  }
  # clustering - x
  if (any(grepl("x", cluster_axis)))
  {
    x_dist = amap::Dist(t(heat_map.mx), method = dist_method)
    x_cluster = hclust(x_dist, method = cluster_method)
    x_dendrogram = as.dendrogram(x_cluster)
    x_dendrogram_reorder = reorder(x_dendrogram, 0, FUN = reorder_method)
    x_order = order.dendrogram(x_dendrogram_reorder)
    heat_map.c = heat_map.mx[,x_order]
  }
  
  heat_map.m = reshape2::melt(heat_map.c) # melted
  heat_map.m = merge(heat_map.m, samples_sheet, by.x = 2, by.y = 1)
  names(heat_map.m) = c("sample", "gene", "expression", "group")
  heat_map.m$group = factor(heat_map.m$group, levels = levels_rug)
  
  # plot
  heat_map = ggplot(data = heat_map.m) + 
    geom_tile(aes(x = sample,
                  y = gene,
                  fill = expression)) + 
    scale_fill_gradientn(colours = colours_heatmap,
                         name = "Expression\n(z-scores)") + 
    labs(title = title, 
         x = "Sample", 
         y = "Gene") + 
    theme(text = element_text(family = "sans"),
          plot.title = element_text(face = "bold", size = 24),
          axis.line = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size = 14, angle = 20, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20),
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 16))
  # make rug
  rug = ggplot(heat_map.m) + 
    geom_bar(aes(x = sample,
                 y = 1,
                 fill = group),
             stat = "identity",
             width = 1) + 
    scale_fill_manual(name = "Group",
                      values = colours_rug,
                      labels = labels_rug) +
    theme(legend.text = element_text(size = 12), 
          legend.title = element_text(size = 16))
  # save legends from heat map and rug
  legend = cowplot::plot_grid(NULL,cowplot::get_legend(heat_map), cowplot::get_legend(rug), NULL, ncol = 1, nrow = 4, axis = "tblr") # add NULL plots to get the legends closer together
  # remove legends from figures to be able to combine plots
    heat_map = heat_map + theme(panel.background = element_rect(fill = "white"),
                              plot.background = element_rect(fill = "white", colour = NA),
                              legend.position = "none")
  rug = rug + 
    theme_void() + 
    theme(panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.spacing.x = unit(1, "mm"),
          legend.position = "none")
  # combine heat map with rug and then add legend to compound plot
  plot = cowplot::plot_grid(heat_map, rug, align = "v", ncol = 1, axis = "tb", rel_heights = c(15, 0.5))
  plot_legend = cowplot::plot_grid(plot, legend, nrow = 1, align = "h", axis = "btlr", rel_widths = c(10,3), rel_heights = c(15, 0.5))
  return(plot_legend)
}

#### PATHWAY ANALYSIS ####

# ORA ####
ora = function(genes_vector, gene_ID_type = "ENSEMBL", organism_db = organism, ontology = "BP", significance = sig_cutoff, qvalue = qvalue_cutoff)
{
  # create sig genes vector
  # genes = subset(master, sig == TRUE)
  # genes_vector = row.names(sig_genes)
  sig_genes_entrez = bitr(genes_vector, 
                          fromType = gene_ID_type, toType = "ENTREZID", 
                          OrgDb = organism_db)
  
  ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID,
                         OrgDb = organism_db,
                         readable = TRUE,
                         ont = ontology,
                         pvalueCutoff = significance,
                         qvalueCutoff = qvalue)
  return(ora_results)
}

get_ora_results = function(ora_results)
{
  gene_set = ora_results$geneID
  description = ora_results$Description
  p.adj = as.numeric(ora_results$p.adjust)
  count = as.numeric(ora_results$Count)
  ora_results_table = data.frame(cbind(gene_set, p.adj, count))
  row.names(ora_results_table) = description
  ora_results_table = ora_results_table[order(ora_results_table$count, decreasing = TRUE),]
  return(ora_results_table)
}

get_summary_ora = function(ora_table, group)
{
  summary_ora = data.frame(ontology = c(row.names(ora_table)),
                           p.adj = as.numeric(ora_table$p.adj),
                           count = as.numeric(ora_table$count),
                           group = group)
  summary_ora = summary_ora[order(summary_ora$count, decreasing = TRUE),]
  summary_ora$ontology = factor(summary_ora$ontology, levels = rev(summary_ora$ontology))
  return(summary_ora)
}

extract_enriched_genes = function(ora_table, gene_set = 1)
{
  ora_table = ora_table[order(ora_table$count, decreasing = TRUE),]
  enriched_gene_set = as.character(ora_table[gene_set,"gene_set"])
  candidate_genes = unlist(strsplit(enriched_gene_set,"/"))
  return(candidate_genes)
}

# Custom bar plot to compare ontology counts ####
plot_multi_ora_barplots = function(ora_summary_all, n_ontologies = 10)
{
  groups = unique(ora_summary_all$group)
  ora_summary_subset = data.frame()
  for(i in 1:length(groups)){
    subset = subset(ora_summary_all, group == paste(groups[i]))
    if(n_ontologies > nrow(subset)) {
      subset = subset[1:nrow(subset),]
    } else{
      subset = subset[1:n_ontologies,]
    }
    subset = subset[order(subset$count, decreasing = TRUE),]
    subset$ontology = factor(subset$ontology, levels = rev(subset$ontology))
    ora_summary_subset = rbind(ora_summary_subset, subset)
  }
  ora_summary_subset = na.omit(ora_summary_subset)
  ggp_barplot = ggplot(ora_summary_subset, aes(x = count,
                                              y = ontology,
                                              fill = p.adj)) + 
    geom_col() + 
    scale_fill_gradientn(name = "Adjusted p-value",
                         colours = colours_p) + 
    geom_text(aes(label = count), nudge_x = 5, size = 7) +
    labs(title = paste("Top",n_ontologies,"Ontologies"),
         x = "Count", 
         y = "Ontology") + 
    facet_wrap(~group, scales = "free_y") +
    scale_y_discrete(labels = function(x) lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse="\n")) +
    theme_bw() +
    my_theme + 
    theme(plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
          strip.text.x = element_text(size = 18), 
          axis.text.y = element_text(size = 19),
          axis.text.x = element_text(size = 18), 
          legend.position = "right")
  return(ggp_barplot)
}

## GSEA ####

gsea = function(data, ontology = "BP", gene_ID_type = "SYMBOL", significance = sig_cutoff, p_adjust_method = "BH", min_gs_size = 3, max_gs_size = 800, organism_db = org_db)
{
  # creating gsea input
  gsea_input = data$log2fold
  names(gsea_input) = row.names(data) # vector will have the gene names associated with the values
  gsea_input = na.omit(gsea_input)
  gsea_input = sort(gsea_input, decreasing = TRUE)
  # analysis
  gsea_results = gseGO(geneList = gsea_input,
                     ont = paste(ontology),
                     keyType = paste(gene_ID_type),
                     nPermSimple = 10000,
                     minGSSize = min_gs_size,
                     maxGSSize = max_gs_size,
                     pvalueCutoff = significance,
                     verbose = TRUE,
                     OrgDb = org.Mm.eg.db,
                     pAdjustMethod = paste(p_adjust_method))
  return(gsea_results)
}


# Plot metagene ####
plot_metagene = function(genes_vector, em_scaled = em.s, samples_sheet = ss, boxplot = TRUE, mean_point = FALSE, jitter = TRUE, labels = labels_groups, colours = colour_groups, title = "Mean Expression")
{
  metagene = data.frame(colMeans(em_scaled[genes_vector,]))
  metagene = merge(metagene, ss, by.x = 0, by.y = 1)
  names(metagene) = c("sample","mean_expression","sample_group")
  
  ggp = ggplot(metagene, aes(x = sample_group,
                               y = mean_expression, 
                               fill = sample_group)) + 
    geom_violin(aes(colour = sample_group), alpha = 0.7, scale = "width", trim = FALSE)
  if (boxplot == TRUE)
  {
    ggp = ggp + geom_boxplot(width = 0.2)
  }
  if (jitter == TRUE)
  {
    ggp = ggp + geom_jitter(aes(fill = sample_group), colour = "black", pch = 21) 
  }
  if (mean_point == TRUE)
  {
    ggp = ggp + stat_summary(fun.data = "mean_cl_boot", geom = "point", colour = "goldenrod1", size = 2)
  }
  ggp = ggp + 
    scale_fill_manual(name = "", 
                      aesthetics = c("colour", "fill"), 
                      values = colours, 
                      labels = labels) + 
    labs(title = title,
         x = "Sample Group", 
         y = "Mean Expression (z-scores)") + 
    theme_bw() + 
    my_theme + 
    theme(legend.position = "none")
  return(ggp)
}

# Fold vs fold scatterplot ####
plot_foldVfold = function(data, x, y, suffix_de1, suffix_de2, colours = colours_sig_groups, label_de1, label_de2, title = "Fold vs Fold Scatterplot")
{
  data$significant = case_when(data[,paste0("sig", suffix_de1)] == TRUE & data[,paste0("sig", suffix_de2)] == TRUE ~ "both",
                               data[,paste0("sig", suffix_de1)] == TRUE & data[,paste0("sig", suffix_de2)] == FALSE ~ suffix_de1,
                               data[,paste0("sig", suffix_de1)] == FALSE & data[,paste0("sig", suffix_de2)] == TRUE ~ suffix_de2,
                               TRUE ~ "none")
  data$significant = factor(data$significant, levels = c(suffix_de1, suffix_de2, "both", "none"))
  
  #  data = subset(data, data$significant != "none")
  cor_res = cor.test(data[,paste0("log2fold", suffix_de1)], data[,paste0("log2fold", suffix_de2)])
  if(cor_res$p.value < 0.001)
  {
    p_value = paste("p < 0.001")
  } else
  {
    p_value = paste("p = ", round(cor_res$p.value,3))
  }
  
  ggp = ggplot(data, aes(x = {{x}},
                         y = {{y}})) + 
    geom_point(aes(colour = significant), alpha = 0.5)  +
    scale_colour_manual(name = "Significant",
                        values = colours, 
                        labels = c(paste(label_de1, "only"), paste(label_de2, "only"), "Both", "None")) + 
    geom_smooth(method = "lm", colour = "goldenrod") +
    annotate("text", x = -5, y = 10,
             label = paste("r = ", round(cor_res$estimate,2),",", p_value), 
             size = 6) +
    labs(title = title,
         x = paste("Log2Fold",label_de1),
         y = paste("Log2Fold",label_de2)) +
    theme_bw() + 
    my_theme +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))
  return(ggp)
}
