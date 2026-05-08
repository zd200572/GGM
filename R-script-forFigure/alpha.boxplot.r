library(vegan)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(ape)
library(picante)
library(dplyr)
library(multcomp)

group.name <- c("Group","Gut")

Loading_Color <- function() {
    cbbPalette1 <- c("#B2182B","#56B4E9","#E69F00","#009E73","#F0E442","#0072B2",
                     "#D55E00","#CC79A7","#CC6666","#9999CC","#66CC99","#999999",
                     "#ADD1E5")
    cbbPalette <- c("#54C0CC","#7EA00E","#1F4F59","#213502")
    result <- list(cbbPalette,cbbPalette1)
    return(result)
}
Loading_Group <- function(file_dir) {
    group <- read.table(file = file_dir,header = TRUE,sep = "\t")
    colnames(group) <- c("variable","Group1","Group2")
    group$Group1 <- factor(group$Group1,levels = unique(group$Group1))
    group$Group2 <- factor(group$Group2,levels = c("Rumen", "Abomasum","Caecum","Colon") )
    return(group)
}
Loading_Table <- function(table_dir) {
    bac <- read.table(file = table_dir,header = FALSE,row.names = 1,sep = "\t",quote = "")
    colnames(bac) <- bac[1,]
    bac <- bac[-1,]
    for (i in 1:(ncol(bac)-1)) {
        bac[,i] <- as.numeric(bac[,i])
    }
    return(bac)
}
cbbPalette <- Loading_Color()
opt <- list(input = "diversity.xls",group = "group.xls")
group <- Loading_Group(opt$group)
result <- Loading_Table(opt$input)
Group_numb1 <- length(unique(group[,2]))
Group_numb2 <- length(unique(group[,3]))
Sample_numb <- length(unique(group[,1]))
l.group1 <- max(str_length(levels(group$Group1)))
l.group2 <- max(str_length(levels(group$Group2)))
result[] <- lapply(result, as.numeric)
alpha1 <- result
alpha1 <- alpha1 %>% rownames_to_column(var = "Sample")

source("diff_alpha.r")
result0 <- diff.alpha(alpha1,group)
pdf("Diff_alpha.pdf",
    width = 2.5 + Group_numb1*Group_numb2*0.35, height = l.group1*0.12 + 4.5)
result0
dev.off()

