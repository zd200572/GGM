library(ggplot2)
library(ggrepel)
library(vegan)

treat_colors <- c(
  Abomasum = "#56B4E9",
  Ceacum   = "#009E73",
  Colon    = "#E69F00",
  Rumen    = "#B2182B"
)
mypch <- c(21:25, 0:14)

da <- read.table("rpkm.xls", sep = "\t", header = TRUE, check.names = FALSE, quote = "")
rownames(da) <- da[,1]
da <- da[,-1]
da <- t(da)
map <- read.table("pca.group.xls", header = TRUE, sep = "\t", check.names = FALSE, quote = "")
rownames(map) <- map[,1]
map$tissu <- factor(map$tissu, levels = c("Abomasum", "Ceacum", "Colon", "Rumen"))
map$group <- factor(map$group)
da <- da[rownames(map), ]
da <- da[, apply(da, 2, function(x) any(x > 0))]
pca <- prcomp(da, scale. = TRUE)
pc_var <- summary(pca)$importance[2, ] * 100
pc_df <- as.data.frame(pca$x[, 1:2])
pc_df$Sample <- rownames(pc_df)
pc_df$Treat  <- map[rownames(pc_df), "tissu"]  
pc_df$Group  <- map[rownames(pc_df), "group"]  
group_levels <- unique(pc_df$Group)
group_shapes <- setNames(mypch[1:length(group_levels)], group_levels)
p1=ggplot(pc_df, aes(x = PC1, y = PC2, fill = Treat, shape = Group)) +
  stat_ellipse(
    geom = "polygon",
    aes(group = Treat, fill = Treat, color = Treat),
    alpha = 0.2,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  geom_point(size = 4, stroke = 0.8, color = "black") +
  scale_fill_manual(values = treat_colors) +
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = group_shapes) +
  theme_bw() +
  labs(
    title = "PCA",
    x = paste0("PC1 (", round(pc_var[1], 2), "%)"),
    y = paste0("PC2 (", round(pc_var[2], 2), "%)")
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.box = "vertical",
    legend.position = "right"
  ) +
  guides(
    fill = guide_legend(title = "Tissue", override.aes = list(shape = 21, color = "black")),
    shape = guide_legend(title = "Group")
  )
  ggsave(p1,filename="pca.pdf")
