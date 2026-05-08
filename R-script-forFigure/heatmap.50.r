library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(multcompView)
library(multcomp)
library(tibble)

abundance_file <- "rpkm.ampID.xls"
group_file <- "group.xls"
topN <- 50  
abund <- read.delim(abundance_file, check.names = FALSE)
group_info <- read.delim(group_file)
rownames(abund) <- abund$ID
abund <- abund[,-1]
common_samples <- intersect(colnames(abund), group_info$SampleID)
abund <- abund[, common_samples]
group_info <- group_info %>% filter(SampleID %in% common_samples)
abund_z <- t(scale(t(as.matrix(abund))))
abund_z <- as.data.frame(abund_z)
abund_long <- abund_z %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "SampleID", values_to = "Zscore")
df <- left_join(abund_long, group_info, by = "SampleID")
top_genes <- df %>%
  group_by(Gene) %>%
  summarise(mean_z = mean(Zscore, na.rm = TRUE)) %>%
  arrange(desc(mean_z)) %>%
  slice(1:topN) %>%
  pull(Gene)
df <- df %>% filter(Gene %in% top_genes)
stat_letters <- data.frame()
for (tis in unique(df$tissu)) {
  df_sub <- df %>% filter(tissu == tis)
  for (gene in unique(df_sub$Gene)) {
    dat <- df_sub %>% filter(Gene == gene)
    dat$group <- factor(dat$group)
    if (n_distinct(dat$group) < 2 || sd(dat$Zscore) == 0) {
      levels_list <- levels(dat$group)
      temp <- data.frame(tissu = tis, Gene = gene, group = levels_list, label = rep("a", length(levels_list)),pval=NaN)
      stat_letters <- rbind(stat_letters, temp)
      next
    }
    if (nlevels(dat$group) == 2) {
      g1 <- dat %>% filter(group == levels(dat$group)[1]) %>% pull(Zscore)
      g2 <- dat %>% filter(group == levels(dat$group)[2]) %>% pull(Zscore)
      if (sd(g1) == 0 && sd(g2) == 0) {
        lbls <- rep("a", 2)
      } else {
        pval <- tryCatch(t.test(g1, g2)$p.value, error = function(e) NA)
        lbls <- if (!is.na(pval) && pval < 0.05) c("a", "b") else c("a", "a")
      }
      temp <- data.frame(tissu = tis, Gene = gene, group = levels(dat$group), label = lbls,pval=pval)
      stat_letters <- rbind(stat_letters, temp)
    } else {
      fit <- aov(Zscore ~ group, data = dat)
      tuk <- tryCatch(glht(fit, linfct = mcp(group = "Tukey")), error = function(e) NULL)
      if (is.null(tuk)) {
        letters <- rep("a", nlevels(dat$group))
        temp <- data.frame(tissu = tis, Gene = gene, group = levels(dat$group), label = letters,pval=NaN)
      } else {
        cld_res <- cld(tuk,decreasing=T)
		tuk_summary <- summary(tuk)
        pvals <- tuk_summary$test$pvalues
        temp <- data.frame(tissu = tis,
                           Gene = gene,
                           group = names(cld_res$mcletters$Letters),
                           label = cld_res$mcletters$Letters,pval=pvals)
      }
      stat_letters <- rbind(stat_letters, temp)
    }
  }
}
write.table(stat_letters,"heatmap.xls",row.names=F,col.names=T,quote=F,sep="\t")
stat_letters <- stat_letters %>%
  group_by(Gene) %>%
  filter(!(all(label == "a"))) %>%
  ungroup() %>%
  group_by(Gene, tissu) %>%
  mutate(
    label = if (all(label == "a")) "" else label
  ) %>%
  ungroup() %>% arrange(Gene,tissu)
df_plot <- df %>%
  group_by(Gene, group, tissu) %>%
  summarise(Zscore = mean(Zscore, na.rm = TRUE), .groups = "drop")
df_plot <- left_join(df_plot, stat_letters, by = c("Gene", "group", "tissu"))
df_plot$group <- factor(df_plot$group, levels = c("14days", "42days","105days"))
df_plot$tissu <-factor(df_plot$tissu,levels= c("Rumen", "Abomasum","Caecum"))
my_palette <- rev(c("#e04c35", "#f68052", "#f3e7a6", "#b9d6e1", "#98bace", "#466c9a"))
z_limit <- max(abs(df_plot$Zscore), na.rm = TRUE)
p <- ggplot(df_plot, aes(x = group, y = Gene)) +
  geom_tile(aes(fill = Zscore), colour = "white") +
  geom_text(aes(label = label), color = "black", fontface = "bold", size = 3) +
  facet_wrap(~tissu, nrow = 1) +
  scale_fill_gradientn(colours = my_palette,
                       limits = c(-z_limit, z_limit),
                       name = "Z-score") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 11),
        axis.text.y = element_text(size = 8, color = "black"),
        strip.text = element_text(face = "bold", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(size = 10)) +
  xlab("") + ylab("")
ggsave("heatmap.pdf", p, width = 10, height = 8)


