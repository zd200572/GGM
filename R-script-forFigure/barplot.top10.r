library(ggplot2)
library(reshape2)
library(ggalluvial)
library(stringr)
library(dplyr)

cbbPalette <- c(
  '#B2182B', '#56B4E9', '#E69F00', '#009E73', '#F0E442', '#0072B2',
  '#D55E00', '#CC79A7', '#CC6666', '#9999CC', '#66CC99', '#999999',
  '#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462',
  '#B3DE69', '#FCCDE5', '#D9D9D9', '#BC80BD', '#CCEBC5', '#FFED6F',
  '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C',
  '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928',
  '#FBB4AE', '#B3CDE3', '#CCEBC5', '#DECBE4', '#FED9A6', '#FFFFCC',
  '#E5D8BD', '#FDDAEC', '#F2F2F2', '#B3E2CD', '#FDCDAC', '#CBD5E8',
  '#F4CAE4', '#E6F5C9'
)

data = read.delim('phylum.abundance.xls', header=TRUE, check.names=FALSE, sep='\t', row.names = 1)
group_info <- read.delim("group.xls", header = TRUE)
df <- as.data.frame(t(data))    
df$SampleID <- rownames(df)
df <- merge(df, group_info, by = "SampleID")
group_avg <- df %>%
  group_by(group, tissu) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")
group_avg$GroupID <- paste0(group_avg$tissu, "_", group_avg$group)
rownames(group_avg) <- group_avg$GroupID
data_avg <- t(group_avg[, !(colnames(group_avg) %in% c("group", "tissu", "GroupID"))])
colnames(data_avg) <- rownames(group_avg)
data_avg <- t(t(data_avg) / colSums(data_avg) * 100)
f.abundance <- as.data.frame(data_avg)
sum <- apply(f.abundance, 1, sum)
f.abundance <- cbind(f.abundance, sum)
f.abundance <- f.abundance[order(f.abundance[, 'sum'], decreasing = TRUE), ]
f.abundance <- subset(f.abundance, select = -sum)
f.abundance <- f.abundance[rownames(f.abundance) != 'Unclassified', ]
f.abundance.1 <- f.abundance
if (nrow(f.abundance.1) > 10) {
  f.abundance.1 <- f.abundance.1[1:10, ]
  f.abundance.1 <- t(f.abundance.1)
  sum2 <- apply(f.abundance.1, 1, sum)
  Others <- 100.00001 - sum2
  f.abundance.1 <- as.data.frame(cbind(f.abundance.1, Others))
}
f.abundance.1$GroupID <- rownames(f.abundance.1)
taxon <- melt(f.abundance.1)
colnames(taxon) <- c("GroupID", "Taxon", "value")
taxon$tissu <- sapply(strsplit(as.character(taxon$GroupID), "_"), `[`, 1)
taxon$group <- sapply(strsplit(as.character(taxon$GroupID), "_"), `[`, 2)
taxon$group <- factor(taxon$group, levels = c("14days", "42days","105days"))
taxon$tissu <- factor(taxon$tissu, levels = c("Rumen", "Abomasum","Caecum"))
p <- ggplot(data = taxon, aes(x = group, y = value, alluvium = Taxon, stratum = Taxon)) +
  geom_alluvium(aes(fill = Taxon), alpha = .5, width = 0.7) +
  geom_stratum(aes(fill = Taxon), width = 0.7) +
  ylab("Relative abundance (%)") + xlab("") +
  scale_fill_manual(values = cbbPalette, name = 'Taxon') +
  facet_wrap(~tissu, scales = "free_x",nrow=1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    panel.background = element_rect(fill = 'transparent', color = 'black'),
    axis.text.x = element_text(colour = 'black', size = 12, face = 'bold'),
    axis.text.y = element_text(colour = 'black', size = 10),
    axis.ticks.length = unit(0.4, 'lines'),
    axis.ticks = element_line(color = 'black'),
    axis.line = element_line(colour = 'black'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    legend.text = element_text(colour = 'black', size = 12),
    legend.title = element_text(size = 14, colour = 'black', face = 'bold')
  ) +
  scale_y_continuous(limits = c(0, 100.001), expand = c(0, 0))
ggsave(plot = p, filename = "phylum_top10.pdf", width = 12, height = 6, useDingbats = FALSE)

