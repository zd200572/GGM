diff.alpha=function(alpha1,group){
    colnames(alpha1)[1] <- "variable"
    alpha1 <- alpha1[,c(1,2,3,4)]
    data <- merge(alpha1,group)
    aa <- c()
    for (j in 1:length(levels(group$Group2))) {
        alpha <- data[data$Group2 == levels(data$Group2)[j],]
        for (i in 2:4) {
            alpha.test <- alpha[,c(i,5)]
            colnames(alpha.test) <- c("Num","Group")
            alpha.test$Group <- factor(alpha.test$Group)
            x <- c("a","b")
            y <- c("a","a")
            if (length(levels(alpha.test$Group)) == 2) {
                fit1 <- t.test(Num~Group,data = alpha.test)
                dd <- alpha.test %>%
                    group_by(Group) %>%
                    summarise(Max = max(Num))
                test <- data.frame(Group1 = levels(alpha.test$Group),
                                   value.x = if(fit1$p.value < 0.05){
                                       x
                                   }else{
                                       y
                                   },
                                   value.y = dd$Max + max(dd$Max)*0.05,
                                   variable = rep(colnames(alpha)[i],2),
                                   Group2 = levels(data$Group2)[j])
            }else{
                fit1 <- aov(Num~Group,data = alpha.test)
                tuk1<-glht(fit1,linfct=mcp(Group="Tukey"))
                res1 <- cld(tuk1)
                dd <- alpha.test %>%
                    group_by(Group) %>%
                    summarise(Max = max(Num))
                test <- data.frame(Group1 = levels(alpha.test$Group),
                                   value.x = res1$mcletters$Letters,
                                   value.y = dd$Max,
                                   variable = rep(colnames(alpha)[i],length(levels(alpha.test$Group))),
                                   Group2 = rep(levels(data$Group2)[j],length(levels(alpha.test$Group))))
            }
            aa <- as.data.frame(rbind(aa,test))
        }
    }

    aa$value.y[aa$variable == "Richness"] <- aa$value.y[aa$variable == "Richness"] + max(aa$value.y[aa$variable == "Richness"])*0.09
    aa$value.y[aa$variable == "Shannon"] <- aa$value.y[aa$variable == "Shannon"] + max(aa$value.y[aa$variable == "Shannon"])*0.035
    aa$value.y[aa$variable == "Simpson"] <- aa$value.y[aa$variable == "Simpson"] + max(aa$value.y[aa$variable == "Simpson"])*0.035

    test <- aa
    alpha.test <- data[,-1]
    alpha.test <- melt(alpha.test)

    alpha.test$variable <- factor(alpha.test$variable,levels = c("Richness","Shannon","Simpson"))
    test$Group2 <- factor(test$Group2,levels = levels(group$Group2))
    test$Group1 <- factor(test$Group1,levels = levels(group$Group1))
    test$variable <- factor(test$variable,levels = c("Richness","Shannon","Simpson"))

    tt <- alpha.test %>%
        group_by(Group1,Group2,variable) %>%
        summarise(Mean = mean(value))

    tt$variable <- factor(tt$variable,levels = c("Richness","Shannon","Simpson"))
    p <- ggplot(alpha.test,aes(Group1,value)) +
        geom_boxplot(aes(color = Group1),width = 0.8,outlier.color = "transparent") +
        geom_jitter(aes(color = Group1),width = 0.4,size = 1.5,alpha = 0.5) +
        geom_line(data = tt,aes(Group1,Mean,group = 1),linewidth = 1) +
        geom_point(data = tt,aes(Group1,Mean),pch = 15,size = 2) +
        geom_text(data = test,aes(x = Group1,y = value.y,label = value.x),
                  size = 4.5,color = "black",fontface = "bold") +
        labs(y = "Variations in alpha diversity",
             x = "") +
        scale_color_manual(values = cbbPalette) +
        facet_grid(variable~Group2,scales = "free_y") +
        theme_bw()+
        theme(panel.grid=element_blank(),
              axis.ticks.length = unit(0.4,"lines"),
              axis.ticks = element_line(color='black'),
              axis.line = element_line(colour = "black"),
              axis.title.x=element_blank(),
              axis.title.y = element_text(face = "bold",color = "black",size = 22),
              axis.text.y=element_text(colour='black',size=10),
              axis.text.x=element_text(colour = "black",size = 12,face = "bold",
                                       angle = 45,vjust = 1,hjust = 1),
              legend.position = "none",
              strip.text = element_text(colour = "black",face = "bold",size = 14),
              strip.text.y = element_text(colour = "black",face = "bold",size = 14))
    return(p)
}

