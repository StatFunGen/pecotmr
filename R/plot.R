#' Venn Diagram
#' @param data a list with the twas siginificant gene_id results of four method "SuSiE","Lasso","Enet" and "MR.ASH"
#' @return plot object
# @importFrom ggvenn ggvenn
# @export
venn <- function(data) {
  venn_plot <- ggvenn(data, c("SuSiE", "Lasso", "Enet", "MR.ASH"), show_percentage = TRUE, fill_color = c("red", "orange", "blue", "green"))
  return(venn_plot)
}


#' Manhattan plot
#' @param twas_results a data frame of twas results with columns "gene_name", "gene_id","chr","susie_pval","lasso_pval","enet_pval" and "mrash_pval", where twas results are the output of the twas_scan function. "gene_name" is the ensemble ID and "gene_id" is the corresponding gene name,
#'        "susie_pval", "lasso_pval","enet_pval" and "mrash_pval" are the pvalues of susie and other three competing twas method.
#' @param gene_data a data frame with columns "chr", "start", "end", and "ID", "chr" is the chromosome of gene, "start" and "end" are the position, "ID" is the gene name.
#' @return plot object
# @import ggplot2
# @importFrom ggrepel geom_label_repel
# @importFrom stringr str_sub
# @importFrom ggnewscale new_scale_color
# @export
manhattan_plot <- function(twas_results, gene_data) {
  min_pval <- apply(twas_results[, c("susie_pval", "lasso_pval", "enet_pval", "mrash_pval")], 1, function(x) min(x, na.rm = TRUE))
  data_all_gene <- twas_results %>%
    select(gene_name, chr, gene_id, susie_pval, lasso_pval, enet_pval, mrash_pval) %>%
    mutate(min_pval = min_pval) %>%
    mutate(chr = as.numeric(chr))
  gene_pos <- gene_data %>%
    mutate(chr = as.numeric(str_sub(`#chr`, 4))) %>%
    select(-`#chr`) %>%
    setNames(c("start_bp", "end_bp", "gene_name", "chr"))
  gene_pos_pval <- merge(data_all_gene, gene_pos, by = c("gene_name", "chr"))

  susie_select <- Mic_genes$gene_pq_adj %>%
    filter(susie_pval < (2.5 * 10^(-6) / 4)) %>%
    select(gene_id) %>%
    t() %>%
    as.vector()
  lasso_select <- Mic_genes$gene_pq_adj %>%
    filter(lasso_pval < (2.5 * 10^(-6) / 4)) %>%
    select(gene_id) %>%
    t() %>%
    as.vector()
  enet_select <- Mic_genes$gene_pq_adj %>%
    filter(enet_pval < (2.5 * 10^(-6) / 4)) %>%
    select(gene_id) %>%
    t() %>%
    as.vector()
  ash_select <- Mic_genes$gene_pq_adj %>%
    filter(mrash_pval < (2.5 * 10^(-6) / 4)) %>%
    select(gene_id) %>%
    t() %>%
    as.vector()
  gene_annotate <- unique(c(susie_select, lasso_select, enet_select, ash_select))

  results <- NULL
  # Define the methods to compare
  methods <- list(
    susie_select, lasso_select, enet_select, ash_select
  )

  # Generate all possible combinations of two methods
  combinations <- combn(methods, 2, simplify = FALSE)

  # Loop through combinations and store the intersections
  for (i in seq_along(combinations)) {
    results[[i]] <- intersect(combinations[[i]][[1]], combinations[[i]][[2]])
  }

  # Combine all unique intersections
  two_more <- unique(unlist(results))

  don <- gene_pos_pval %>%
    # Compute chromosome size
    group_by(chr) %>%
    summarise(chr_len = max((start_bp + end_bp) / 2)) %>%
    # Calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(gene_pos_pval, ., by = c("chr" = "chr")) %>%
    # Add a cumulative position of each SNP
    arrange(chr, (start_bp + end_bp) / 2) %>%
    mutate(BPcum = (start_bp + end_bp) / 2 + tot) %>%
    mutate(susie_highlight = ifelse(gene_id %in% susie_select, "yes", "no")) %>%
    mutate(lasso_highlight = ifelse(gene_id %in% lasso_select, "yes", "no")) %>%
    mutate(enet_highlight = ifelse(gene_id %in% enet_select, "yes", "no")) %>%
    mutate(ash_highlight = ifelse(gene_id %in% ash_select, "yes", "no")) %>%
    mutate(two_highlight = ifelse(gene_id %in% two_more, "yes", "no")) %>%
    mutate(gene_annotate = ifelse(gene_id %in% gene_annotate, "yes", "no"))

  axisdf <- don %>%
    group_by(chr) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)

  manhattan <- ggplot(don, aes(x = BPcum, y = -log10(min_pval))) +

    # Show all points
    geom_point(aes(color = as.factor(chr), show.legend = FALSE), alpha = 0.8, size = 5) +
    scale_color_manual(values = rep(c("#b6b6b6", "#9f9f9f"), 22), guide = FALSE) +
    # custom X axis:
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0)) + # remove space between plot area and x axis
    ylim(0, 90) +
    new_scale_color() +
    geom_point(data = subset(don, susie_highlight == "yes"), aes(color = "red"), size = 5) +
    geom_point(data = subset(don, lasso_highlight == "yes"), aes(color = "orange"), size = 5) +
    geom_point(data = subset(don, enet_highlight == "yes"), aes(color = "green"), size = 5) +
    geom_point(data = subset(don, ash_highlight == "yes"), aes(color = "blue"), size = 5) +
    geom_point(data = subset(don, two_highlight == "yes"), aes(color = "black"), size = 5) +
    geom_label_repel(data = subset(don, gene_annotate == "yes"), aes(label = gene_id), size = 5, max.overlaps = Inf, color = "black", alpha = 0.75) +
    scale_color_manual(values = c("red", "orange", "green", "blue", "black"), breaks = c("red", "orange", "green", "blue", "black"), name = "Methods", labels = c("SuSiE", "Lasso", "Enet", "MR.ASH", ">=2 methods")) +
    xlab("Chromosome") +
    ylab("-log10(pval)") +
    geom_hline(yintercept = -log10(2.5 * 10^(-6) / 4), linetype = "dashed", color = "red") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 45),
      legend.title = element_text(face = "bold", size = 45),
      legend.text = element_text(size = 30),
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 45),
      legend.position = "top",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  return(manhattan)
}
