correlation.dotplot <- function(otu.wilcox.biomarker,
                                metabolites.wilcox.biomarker,
                                otu,
                                metabolites,
                                title = "Correlation Dotplot with Grid") {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  
otu <- otu[, 1:(ncol(otu) - 1)]
  otu <- otu / colSums(otu) * 100
  colnames(otu.wilcox.biomarker)[1] <- "V1"
  colnames(metabolites.wilcox.biomarker)[1] <- "V1"
  otu1 <- otu[rownames(otu) %in% otu.wilcox.biomarker$V1, ]
  metabolites1 <- metabolites[rownames(metabolites) %in% metabolites.wilcox.biomarker$V1,
                              1:(ncol(metabolites) - 1)]
  common_samples <- intersect(colnames(otu1), colnames(metabolites1))
  otu1 <- otu1[, common_samples]
  metabolites1 <- metabolites1[, common_samples]
  cor_result <- cor(t(metabolites1), t(otu1), method = "spearman")
  p_result <- matrix(NA, nrow = nrow(cor_result), ncol = ncol(cor_result))
  for (i in 1:nrow(cor_result)) {
    for (j in 1:ncol(cor_result)) {
      test <- cor.test(
        x = as.numeric(metabolites1[i, ]),
        y = as.numeric(otu1[j, ]),
        method = "spearman",
        exact = FALSE  
      )
      p_result[i, j] <- test$p.value
    }
  }
  rownames(p_result) <- rownames(cor_result)
  colnames(p_result) <- colnames(cor_result)
  df <- cor_result %>%
    as.data.frame() %>%
    rownames_to_column("Metabolite") %>%
    pivot_longer(-Metabolite, names_to = "OTU", values_to = "Correlation")
  df$pvalue <- as.vector(p_result)
  sig_otus <- df %>%
    filter(pvalue < 0.05) %>%
    pull(OTU) %>%
    unique()
  df <- df %>% filter(OTU %in% sig_otus)  
  df$p_label <- sprintf("%.2f", df$pvalue)  
  df$OTU <- factor(df$OTU, levels = unique(df$OTU))
  df$Metabolite <- factor(df$Metabolite, levels = rev(unique(df$Metabolite)))
  p <- ggplot(df, aes(x = OTU, y = Metabolite)) +
    geom_tile(fill = "white", color = "grey90", linewidth = 0.2) +
    geom_point(aes(size = abs(Correlation), color = Correlation), alpha = 0.8) +
    geom_text(aes(label = p_label), size = 2.2, color = "black", hjust = 0.5, vjust = 0.5) +
    scale_color_gradient2(
      low = "#5861AC", mid = "grey", high = "#F28080",
      limits = c(-1, 1), name = "Cor"
    ) +
    scale_size_continuous(
      range = c(3, 10), name = "|Cor|"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10)
    ) +
    labs(x = NULL, y = NULL, title = title)
  return(list(cor = cor_result, p = p_result, plot = p, dataf = df))
}
