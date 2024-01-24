library(vegan)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(hillR)
library(ggpubr)
library(ggpmisc)
library(stringr)
library(cowplot)
library(ggtext)
#references
#https://rpubs.com/an-bui/vegan-cheat-sheet
#https://sites.google.com/site/mb3gustame/home
#https://www.flutterbys.com.au/stats/tut/tut13.2.html

metadata <- read.csv("disturbance_total_metadata_no2CDay5.csv", header=TRUE)
otu <- read.csv("disturbance_total_otu_genus_no2CDay5.csv", header=TRUE, check.names = FALSE)

#metadata$Frequency <- as.character(metadata$Frequency)
#metadata$Frequency <- factor(metadata$Frequency, levels=c("7","5","3","2","1"))
metadata$Substrate <- factor(metadata$Substrate, levels=c("Cellulose","Glucose"))

#attach otu to metadata
data <- cbind(metadata,otu)
#drop columns
data <-data[-c(8)]
#data$Frequency <- as.character(data$Frequency)




######## Calculate Hill numbers as suggested by Reviewer #2 #####################

comm <- data[,c(1,8:173)]
row.names(comm) <- comm[,1]
comm <- comm[,-c(1)]

#h0 is richness
#h1 is exponent of shannon entropy -> "number of common ASVs"
#h2 is inverse simpson

h0 <- as.data.frame(hill_taxa(comm, q = 0, MARGIN = 1, base = exp(1)))
h0$Sample_Name <- row.names(h0)
h1 <- as.data.frame(hill_taxa(comm, q = 1, MARGIN = 1, base = exp(1)))
h1$Sample_Name <- row.names(h1)
h2 <- as.data.frame(hill_taxa(comm, q = 2, MARGIN = 1, base = exp(1)))
h2$Sample_Name <- row.names(h2)


hill_total <- merge(h0, h1, by = "Sample_Name")
hill_total <- merge(hill_total, h2, by = "Sample_Name")
colnames(hill_total) <- c("Sample_Name", "Hill0","Hill1", "Hill2")

hill_total <- merge(hill_total, metadata, by = "Sample_Name")

#write.csv(hill_total, "diversity_metadata.csv", row.names = FALSE)

hill_total$Frequency2 <- 1/hill_total$Frequency

hill_summary <- hill_total %>%
  group_by(Substrate, Frequency) %>%
  summarise(meanHill1 = mean(Hill1),
            sdHill1 = sd(Hill1),
            meanHill2 = mean(Hill2),
            sdHill2 = sd(Hill2),
            meanHill0 = mean(Hill0),
            sdHill0 = sd(Hill0))



# do some t.tests to claim 1/3 is highest for cellulose and 1/1 is highest for glucose

datacel <- subset(hill_total, Substrate == "Cellulose")
dataglu <- subset(hill_total, Substrate == "Glucose")

tab1<-compare_means(Hill1 ~ Frequency, data=datacel, ref.group="3", method="t.test")
tab1$Substrate <- "Cellulose"


tab2<-compare_means(Hill1 ~ Frequency, data=dataglu, ref.group="1", method="t.test")
tab2$Substrate <- "Glucose"

tab_total <- rbind(tab1,tab2)
colnames(tab_total)[1] = "DiversityIndex"
colnames(tab_total)[2] = "Reference"
colnames(tab_total)[3] = "Comparison"
colnames(tab_total)[4] = "p-value"

#write.csv(tab_total, "disturbance_diversity_stats.csv", row.names = FALSE)

tab_total2 <- tab_total %>%
  mutate(across("Reference", str_replace, "1", "1/1"),
         across("Reference", str_replace, "3", "1/3"),
         across("Comparison", str_replace, "1", "1/1"),
         across("Comparison", str_replace, "2", "1/2"),
         across("Comparison", str_replace, "3", "1/3"),
         across("Comparison", str_replace, "5", "1/5"),
         across("Comparison", str_replace, "7", "1/7"))

tab_total2 <- tab_total2[,c(9,8,2,3,6)]
colnames(tab_total2) <- c("Substrate", "Method","Reference \n Frequency","Comparison \n Frequency","p-value")

tab3 <- compare_means(Hill1 ~ Substrate, data=hill_total, group.by = "Frequency", method = "t.test")
tab3$Frequency <- c("1/1", "1/2", "1/3", "1/5", "1/7")
tab3 <- tab3[,c(1,9,3,4,7)]
tab3 <- tab3[c(5,4,3,2,1),]
colnames(tab3) <- c("Frequency", "Method","Reference","Comparison","p-value")

###### Graph Hill numbers #######

## Start with H1 


hill_total$Frequency2 <- as.factor(hill_total$Frequency2)

stats_data <- function(y) {
  return(data.frame(
    y=7.8,
    label = paste('n=', length(y), '\n')
  ))
}

pHill1 <-ggplot(hill_total, aes(x = Frequency2, y = Hill1, color = Substrate, fill = Substrate)) +
  geom_boxplot() +
  stat_summary(fun.data = stats_data,
               geom = "text",
               position = position_dodge(width = 0.9),
               size = 2.5) +
  theme_classic(base_size=10)+
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  scale_x_discrete(#breaks = c(0.14285714285,0.2,0.33333333333,0.5,1),
                     labels = c("1/7","1/5","1/3","1/2","1/1")) +
  labs(y="Number of Common ASVs", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Hill 1 Diversity \n Exponential of Shannon's Entropy")
#pHill1


p_stats <- ggtexttable(tab_total2, rows = NULL,
                            theme = ttheme("blank",  base_size = 7.5,
                                           padding = unit(c(1.5,3), "mm"))) 

p_stats <- tab_add_title(p_stats, "Comparison of Hill 1 means of \n selected frequency treatments", size = 12)

p_stats2 <- ggtexttable(tab3, rows = NULL,
                       theme = ttheme("blank",  base_size = 7.5,
                                      padding = unit(c(1.5,3), "mm"))) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
p_stats2 <- tab_add_title(p_stats2, "Comparison of Hill 1 means \n between substrate treatments", size = 12)


p1 <- pHill1 
p2 <- p_stats
p3 <- p_stats2

fig4 <- ggdraw () +
  draw_plot(p1, x = 0, y = 0, width = .55, height = 1) +
  draw_plot(p2, x = 0.65, y = 0.245, width = .3, height = 1,) +
  draw_plot(p3, x = 0.65, y = -.245, width = .3, height = 1) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.556, 0.56), y = c(1, 1, 0.44))
#fig4

ggsave("Figure_4_hill1.tiff", width = 6.875, height = 4.8, units = "in", device = "tiff", dpi = 300)
#ggsave("Figure_hill1.png", device = "png", dpi = 300)
#ggsave("Figure_hill1.pdf", device = "pdf", dpi = 300)




# graph H0 and H2 for supplemental

hill_total$Frequency2 <- as.factor(hill_total$Frequency2)

stats_data1 <- function(y) {
  return(data.frame(
    y=55,
    label = paste('n=', length(y), '\n')
  ))
}

pHill0 <-ggplot(hill_total, aes(x = Frequency2, y = Hill0, color = Substrate, fill = Substrate)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test", label = "p.signif", label.x = c(1.5))+
  stat_summary(fun.data = stats_data1,
               geom = "text",
               position = position_dodge(width = 0.9)) +
  theme_classic(base_size=15)+
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  scale_x_discrete(#breaks = c(0.14285714285,0.2,0.33333333333,0.5,1),
    labels = c("1/7","1/5","1/3","1/2","1/1")) +
  labs(y="Number of ASVs", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Hill 0 Diversity \n Richness")
pHill0

stats_data2 <- function(y) {
  return(data.frame(
    y=5.5,
    label = paste('n=', length(y), '\n')
  ))
}

pHill2 <-ggplot(hill_total, aes(x = Frequency2, y = Hill2, color = Substrate, fill = Substrate)) +
  geom_boxplot() +
  stat_compare_means(method = "t.test", label = "p.signif", label.x = c(1.5))+
  stat_summary(fun.data = stats_data2,
               geom = "text",
               position = position_dodge(width = 0.9)) +
  theme_classic(base_size=15)+
  theme(legend.position = "bottom") +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  scale_x_discrete(#breaks = c(0.14285714285,0.2,0.33333333333,0.5,1),
    labels = c("1/7","1/5","1/3","1/2","1/1")) +
  labs(y="Diversity (Hill 2, Inverse-Simpson)", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Hill 2 Diversity \n Inverse-Simpson Index")
pHill2

ggarrange(pHill0, pHill2, common.legend = TRUE, legend = "bottom", labels = "AUTO", font.label=list(size = 20))

ggsave("Supplemental_Figure_2_Hill0_Hill2.tiff", device = "tiff", width = 10.8, height = 6.5, units = "in", dpi = 300)
#ggsave("Hill0_Hill2.png", device = "png", dpi = 300)
#ggsave("Hill0_Hill2.pdf", device = "pdf", dpi = 300)


#lets do all three for fun

hill_table <- hill_total[,c(1:6,8:10)]
hill_table <- melt(hill_table, id =c("Sample_Name","Frequency","Substrate","Replicate","Time","Run"))


phill <- ggplot(hill_table, aes(x = variable, y = value, color = Substrate)) +
  geom_point() +
  facet_grid(~Frequency)

phill




## regression
## Star with H1 cause that's probably the best to compare

hill_total$Frequency2 <- as.numeric(hill_total$Frequency2)
pHill1_lin <-ggplot(hill_total, aes(x = Frequency2, y = Hill1, color = Substrate)) +
  geom_point() +
  ylim(1,9.5) +
  stat_poly_line(formula = y ~ x, linetype = "dotted", se = T) +
  stat_poly_eq(use_label(c("eq", "R2")),
               formula = y ~ x,
               label.x = "left",
               label.y = "top") +
  theme_classic(base_size=13)+
  #theme(legend.position = "none") +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  scale_x_continuous(breaks = c(0.14285714285,0.2,0.33333333333,0.5,1),
                     labels = c("1/7","1/5","1/3","1/2","1/1")) +
  labs(y="Number of common ASVs", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Hill 1, Linear Regression")
pHill1_lin


pHill1_poly <- ggplot(hill_total, aes(x = Frequency2, y = Hill1, color = Substrate)) +
  geom_point() +
  ylim(1,9.5) +
  stat_poly_line(formula = y ~ poly(x,2), linetype = "dotted", se = T) +
  stat_poly_eq(use_label(c("eq", "R2")),
               formula = y ~ poly(x, 2),
               label.x = "left",
               label.y = "top") +
  theme_classic(base_size=13)+
  #theme(legend.position = "none") +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  scale_x_continuous(breaks = c(0.14285714285,0.2,0.33333333333,0.5,1),
                     labels = c("1/7","1/5","1/3","1/2","1/1")) +
  labs(y="Number of common ASVs", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Hill 1, Quadratic Regression")
#pHill1_poly

preg <- ggarrange(pHill1_lin, pHill1_poly, common.legend = TRUE, legend = "right", labels = "AUTO", font.label=list(size = 20) )
preg


x <- c(1/1, 1/2, 1/3, 1/5, 1/7)
c_lin <- 4.31 -(1.26*x)
g_lin <- 1.51 + (1.89*x)
c_quad <- 3.8 - (4.68*x) - (5.57 * (x^2))
g_quad <- 2.25 + (6.54*x) + (1.26 * (x^2))

predicted <- data.frame(x, c_lin, g_lin, c_quad, g_quad)
predicted <- melt(predicted, id = c("x"))
Substrate_Frequency <- c("1C","2C","3C","5C","7C","1G","2G","3G","5G","7G","1C","2C","3C","5C","7C","1G","2G","3G","5G","7G")
predicted <-cbind(predicted, Substrate_Frequency)

hill_total <- merge(hill_total, predicted, by = "Substrate_Frequency")

hill_total <- hill_total %>%
  mutate(variable = str_replace(variable, "g_|c_", ""))

hill_total <- hill_total %>%
  mutate(variable = str_replace(variable, "lin", "Linear"))

hill_total <- hill_total %>%
  mutate(variable = str_replace(variable, "quad", "Quadratic"))

hill_total$residual <- (hill_total$Hill1) - (hill_total$value)

hill_total$Frequency <- factor(hill_total$Frequency, levels = c("7", "5", "3", "2", "1"))


p_resid <- ggplot(hill_total, aes(x = value, y = residual, shape = variable, color = Frequency)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_classic(base_size = 13) + 
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854")) +
  labs(y="Residuals", x="Predicted Hill 1") + 
  guides(shape=guide_legend("Regression")) +
  facet_grid(~Substrate)
  
p_resid

ggarrange(preg, p_resid, ncol = 1, labels = c("", "C"), font.label = list(size = 20))

ggsave("Supplemental_Figure_4_regression.tiff", width = 8, height = 8, units = "in", device = "tiff", dpi = 300)
#ggsave("regression.png", device = "png", dpi = 300)
#ggsave("regression.pdf", device = "pdf", dpi = 300)



###### Abundance Bar graph ######

metadata <- read.csv("disturbance_total_metadata_top15_no2CDay5.csv", header=TRUE)
otu <- read.csv("disturbance_total_otu_genus_top15_no2CDay5.csv", header=TRUE, check.names = FALSE)


data_labels <- c("1C_A_1" = "1", "1C_A_2" = "2", "1C_A_5" = "5", "1C_A_6" = "6", "1C_A_7" = "7", 
                 "1C_B_1" = "1", "1C_B_2" = "2", "1C_B_4" = "4", "1C_B_5" = "5", "1C_B_7" = "7", 
                 "1C_C_1" = "1", "1C_C_2" = "2", "1C_C_4" = "4", "1C_C_6" = "6", "1C_C_7" = "7", 
                 "1C_D_1" = "1", "1C_D_2" = "2", "1C_D_3" = "3", "1C_D_4" = "4", "1C_D_5" = "5", "1C_D_6" = "6", 
                 "1C_E_1" = "1", "1C_E_2" = "2", "1C_E_3" = "3", "1C_E_4" = "4","1C_E_5" = "5", "1C_E_6" = "6", "1G_A_2" = "2", "1G_A_4" = "4", "1G_A_5" = "5", "1G_A_6" = "6", 
                 "1G_B_2" = "2", "1G_B_4" = "4", "1G_B_5" = "5", "1G_B_6" = "6", "1G_B_7" = "7", "1G_C_4" = "4", "1G_C_5" = "5", "1G_C_7" = "7", 
                 "1G_D_1" = "1", "1G_D_2" = "2", "1G_D_4" = "4", "1G_D_6" = "6", "1G_D_7" = "7", "1G_E_1" = "1", "1G_E_2" = "2", "1G_E_4" = "4", "1G_E_5" = "6", "1G_E_6" = "6", "2C_A_2" = "2", 
                 "2C_B_2" = "2", "2C_B_3" = "3", "2C_C_1" = "1", "2C_C_3" = "3", "2C_D_1" = "1", "2C_D_2" = "2", "2C_D_3" = "3", "2C_E_1" = "1", "2C_E_2" = "2", "2G_A_1" = "1", "2G_A_2" = "2", 
                 "2G_B_1" = "1", "2G_E_3" = "3", "3C_A_7" = "7", "3C_B_7" = "7", "3C_C_7" = "7", "3C_D_7" = "7", "3C_E_6" = "6", "3C_E_7" = "7", "3G_B_7" = "7", "3G_C_7" = "7", "3G_D_7" = "7", 
                 "5C_A_1" = "1", "5C_A_2" = "2", "5C_A_3" = "3", "5C_A_4" = "4", "5C_A_5" = "5", "5C_A_6" = "6", "5C_A_7" = "7", "5C_B_2" = "2", "5C_B_3" = "3", "5C_B_4" = "4", "5C_B_5" = "5", 
                 "5C_B_6" = "6", "5C_B_7" = "7", "5C_C_1" = "1", "5C_C_2" = "2", "5C_C_3" = "3", "5C_C_4" = "4", "5C_C_5" = "5", "5C_C_6" = "6", "5C_C_7" = "7", "5C_D_1" = "1", "5C_D_2" = "2", 
                 "5C_D_3" = "3", "5C_D_4" = "4", "5C_D_5" = "5", "5C_D_6" = "6", "5C_D_7" = "7", "5C_E_1" = "1", 
                 "5C_E_2" = "2", "5C_E_3" = "3", "5C_E_4" = "4", "5C_E_5" = "5", "5C_E_6" = "6", "5C_E_7" = "7", "5G_A_1" = "1", "5G_A_2" = "2", "5G_A_3" = "3", "5G_A_4" = "4", "5G_A_5" = "5", 
                 "5G_A_6" = "6", "5G_A_7" = "7", "5G_B_1" = "1", "5G_B_2" = "2", "5G_B_3" = "3", "5G_B_4" = "4", "5G_B_5" = "5", "5G_B_6" = "6", "5G_B_7" = "7", "5G_C_1" = "1", "5G_C_2" = "2", 
                 "5G_C_3" = "3", "5G_C_4" = "4", "5G_C_5" = "5", "5G_C_6" = "6", "5G_C_7" = "7", "5G_D_1" = "1", "5G_D_2" = "2", "5G_D_3" = "3", "5G_D_4" = "4", "5G_D_5" = "5", "5G_D_6" = "6", "5G_D_7" = "7", "5G_E_1" = "1", "5G_E_2" = "2", 
                 "5G_E_3" = "3", "5G_E_4" = "4", "5G_E_5" = "5", "5G_E_6" = "6", "5G_E_7" = "7", "7C_A_1" = "1", "7C_A_2" = "2", "7C_B_1" = "1", "7C_C_1" = "1", "7C_D_1" = "1", "7C_E_1" = "1", "7G_A_1" = "1", "7G_A_2" = "2", "7G_B_1" = "1", 
                 "7G_B_2" = "2", "7G_C_2" = "2", "7G_D_2" = "2", "7G_E_1" = "1", "7G_E_2" = "2", "2C_A_4" = "4", "2C_A_6" = "6", "2C_A_7" = "7", "2C_B_4" = "4", "2C_B_6" = "6", "2C_C_4" = "4", 
                 "2C_C_6" = "6", "2C_D_4" = "4", "2C_D_6" = "6", "2C_D_7" = "7", "2C_E_3" = "3", "2C_E_6" = "6", "2C_E_7" = "7", "2G_A_3" = "3", "2G_A_4" = "4", "2G_A_6" = "6", "2G_A_7" = "7", "2G_B_3" = "3", "2G_B_4" = "4", "2G_B_6" = "6", 
                 "2G_B_7" = "7", "2G_C_2" = "2", "2G_C_3" = "2", "2G_C_4" = "4", "2G_C_6" = "6", "2G_C_7" = "7", "2G_D_2" = "2", "2G_D_3" = "3", "2G_D_4" = "4", "2G_D_5" = "5", "2G_D_6" = "6", 
                 "2G_D_7" = "7", "2G_E_2" = "2", "2G_E_4" = "4", "2G_E_5" = "5", "2G_E_6" = "6", "2G_E_7" = "7", "3C_A_1" = "1", "3C_A_2" = "2", "3C_A_3" = "3", "3C_A_4" = "4", "3C_A_5" = "5", 
                 "3C_A_6" = "6", "3C_B_1" = "1", "3C_B_2" = "2","3C_B_3" = "3", "3C_B_4" = "4", "3C_B_5" = "5", "3C_B_6" = "6", "3C_C_1" = "1", "3C_C_2" = "2", "3C_C_3" = "3", "3C_C_4" = "4", "3C_C_5" = "5", "3C_C_6" = "6", "3C_D_1" = "1", 
                 "3C_D_2" = "2", "3C_D_3" = "3", "3C_D_4" = "4", "3C_D_5" = "5", "3C_D_6" = "6", "3C_E_2" = "2", "3C_E_3" = "3", "3C_E_4" = "4", "3C_E_5" = "5", "3G_A_1" = "1", "3G_A_2" = "2", "3G_A_3" = "3", "3G_A_4" = "4", "3G_A_5" = "5", 
                 "3G_A_6" = "6", "3G_A_7" = "7", "3G_B_1" = "1", "3G_B_3" = "3", "3G_B_4" = "4", "3G_B_5" = "5", "3G_B_6" = "6", "3G_C_2" = "2", "3G_C_3" = "3", "3G_C_4" = "4", "3G_C_5" = "5", 
                 "3G_C_6" = "6", "3G_D_1" = "1", "3G_D_2" = "2", "3G_D_3" = "3", "3G_D_4" = "4", "3G_D_5" = "5", "3G_D_6" = "6", "3G_E_1" = "1", "3G_E_2" = "2", "3G_E_3" = "3", "3G_E_4" = "4", "3G_E_5" = "5", "3G_E_6" = "6", "7C_A_3" = "3",
                 "7C_A_4" = "4", "7C_A_5" = "5", "7C_A_6" = "6", "7C_A_7" = "7", "7C_B_2" = "2", "7C_B_3" = "3", "7C_B_4" = "4", "7C_B_5" = "5", "7C_B_6" = "6", "7C_B_7" = "7", "7C_C_2" = "2", 
                 "7C_C_3" = "3", "7C_C_4" = "4", "7C_C_5" = "5", "7C_C_6" = "6", "7C_C_7" = "7", "7C_D_2" = "2", "7C_D_3" = "3", "7C_D_4" = "4", "7C_D_5" = "5", "7C_D_6" = "6", "7C_D_7" = "7", "7C_E_2" = "2", "7C_E_3" = "3", "7C_E_5" = "5",
                 "7C_E_6" = "6", "7C_E_7" = "7", "7G_A_3" = "3", "7G_A_4" = "4", "7G_A_5" = "5", "7G_A_6" = "6", "7G_A_7" = "7", "7G_B_3" = "3", "7G_B_4" = "4", "7G_B_5" = "5", "7G_B_6" = "6", 
                 "7G_B_7" = "7", "7G_C_3" = "3", "7G_C_4" = "4", "7G_C_5" = "5", "7G_C_6" = "6", "7G_C_7" = "7", "7G_D_3" = "3", "7G_D_4" = "4", "7G_D_5" = "5", "7G_D_6" = "6", "7G_D_7" = "7", "7G_E_3" = "3", "7G_E_4" = "4", "7G_E_5" = "5",
                 "7G_E_6" = "6", "7G_E_7" = "7",
                 "2C_A_5" = "5","2C_B_5" = "5","2C_C_5" = "5","2C_D_5" = "5","2C_E_5" = "5")


#order samples with separate file that I worked with in excel
order <-read.csv("disturbance_total_order_top15.csv", header=TRUE)
metadata$Sample_Name <- factor(metadata$Sample_Name, levels=unique(order$Sample_Name))

# strucutre factors for easier graphing
metadata$Frequency <- as.character(metadata$Frequency)
metadata$Frequency <- factor(metadata$Frequency, levels=c("7","5","3","2","1"))
metadata$Substrate <- factor(metadata$Substrate, levels=c("inoculum","Cellulose","Glucose"))

#attach otu to metadata
data <- cbind(metadata,otu)
#drop redundant column
data <-data[-c(8)]


data <- melt(data,id=c("Sample_Name","Substrate","Frequency","Substrate_Frequency","Replicate","Time","Run"))
#data<- within(data, Time_Replicate <-paste(Time,Replicate, sep=""))
#make it a nice looking % number
data$value <- data$value*100
head(data)

colorCount <- 15
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

data_summary <- data %>%
  group_by(Substrate, Frequency, variable, Time) %>%
  summarise(count=n(),
            mean=mean(value),
            median=medium(value),
            min=min(value),
            max=max(value))
data_sum_g <- subset(data_summary, Substrate == "Glucose")
data_sum_g5 <- subset(data_sum_g, Frequency == "5")
data_sum_c <- subset(data_summary, Substrate == "Cellulose")
data_sum_c2 <- subset(data_sum_c, Frequency == "2")

datac <- subset(data, data$Substrate=="Cellulose")
datag <- subset(data, data$Substrate=="Glucose")



## Bar graphs
pg <-ggplot(datag) + 
  geom_bar(aes(y=value,x=Sample_Name,fill=variable, color = variable),stat="identity") + 
  theme_classic(base_size=13) + #make it pretty
  theme(axis.text.x=element_text(size = 8.5),
    legend.position="bottom",
    legend.text = element_markdown()) + 
  ylab("Relative Abundance (%)") + 
  xlab("Samples arranged by Time (Day Sampled) in replicates") +
  labs(fill='Genus')+ #change the legend title
  scale_fill_manual(values=getPalette(colorCount),
                    labels = c("*Pseudomonas*.1","*Asticcacaulis*","*Cellvibrio*", 
                               "*Pseudomonas*.2","*Cohnella*","*Lacunisphaera*",
                               "*Taibaiella*","*Achromobacter*","*Leadbetterella*",
                               "*Bordetella*", "*Shinella*", "*Xylophilus*", 
                               "*Methylophilus*","*Stenotrophomonas*","*Pleomorphomonas*")) +
  scale_color_manual(values=getPalette(colorCount)) +
  scale_x_discrete(labels=data_labels) +
  facet_grid( . ~ Frequency, scales = "free_x") + 
  ggtitle("Glucose Samples faceted by Disturbance Frequency (1/n Days)") +
  guides(color="none", fill=guide_legend("Genus"))

pc <-ggplot(datac) + 
  geom_bar(aes(y=value,x=Sample_Name,fill=variable, color = variable),stat="identity") + 
  theme_classic(base_size=13) + #make it pretty
  theme(axis.text.x=element_text(size = 8.5),
    legend.position="bottom",
    legend.text = element_markdown()) + 
  ylab("Relative Abundance (%)") + 
  xlab("Samples arranged by Time (Day Sampled)in replicates") +
  labs(fill='Genus')+ #change the legend title
  scale_fill_manual(values=getPalette(colorCount),
                    labels = c("*Pseudomonas*.1","*Asticcacaulis*","*Cellvibrio*", 
                               "*Pseudomonas*.2","*Cohnella*","*Lacunisphaera*",
                               "*Taibaiella*","*Achromobacter*","*Leadbetterella*",
                               "*Bordetella*", "*Shinella*", "*Xylophilus*", 
                               "*Methylophilus*","*Stenotrophomonas*","*Pleomorphomonas*")) +
  scale_color_manual(values=getPalette(colorCount)) +
  scale_x_discrete(labels=data_labels) +
  facet_grid( . ~ Frequency , scales = "free_x")  + 
  ggtitle("Cellulose Samples faceted by Disturbance Frequency (1/n Days)") +
  guides(color="none", fill=guide_legend("Genus"))


ggarrange(pc + rremove("xlab"), pg, common.legend = TRUE, ncol=1, nrow=2,legend = "bottom")

ggsave("Figure_2_abundance.tiff", width = 10.3, height = 6.8, units = "in", device = "tiff", dpi = 300)
#ggsave("abundance.png", device = "png", dpi = 300)
#ggsave("abundance.pdf", device = "pdf", dpi = 300)


######### BAR GRAPH of just DISTURBANCE FREQUENCY 2 FOR SUPPLEMENTAL ############

metadata2 <- read.csv("disturbance_total_metadata.csv", header=TRUE)
otu2 <- read.csv("disturbance_total_otu_genus_top15.csv", header=TRUE, check.names = FALSE)

metadata2$Sample_Name <- factor(metadata2$Sample_Name, levels=unique(order$Sample_Name))



metadata2$Substrate <- factor(metadata2$Substrate, levels=c("Cellulose","Glucose"))

#attach otu to metadata
data2 <- cbind(metadata2,otu2)
#drop columns
data2 <-data2[-c(8)]
data2 <- subset(data2, Frequency == "2")
data2 <- melt(data2,id=c("Sample_Name","Substrate","Frequency","Substrate_Frequency","Replicate","Time","Run"))
#data2 <- subset(data2, variable %in% tax_list)
data2$value <- data2$value * 100


p_f2 <-ggplot(data2) + 
  geom_bar(aes(y=value,x=Sample_Name,fill=variable),stat="identity") + 
  theme_classic(base_size=15) + #make it pretty
  theme(axis.text.x=element_text(size = 10),
        legend.position="bottom",
        legend.text = element_markdown()) + 
  ylab("Relative Abundance (%)") + 
  xlab("Samples arranged by Time (Day Sampled)in replicates") +
  labs(fill='Genus')+ #change the legend title
  scale_fill_manual(values=getPalette(colorCount),
                    labels = c("*Pseudomonas*.1","*Asticcacaulis*","*Cellvibrio*", 
                               "*Pseudomonas*.2","*Cohnella*","*Lacunisphaera*",
                               "*Taibaiella*","*Achromobacter*","*Leadbetterella*",
                               "*Bordetella*", "*Shinella*", "*Xylophilus*", 
                               "*Methylophilus*","*Stenotrophomonas*","*Pleomorphomonas*")) +
  scale_color_manual(values=getPalette(colorCount)) +
  scale_x_discrete(labels=data_labels) +
  facet_grid( . ~ Frequency , scales = "free_x")  + 
  ggtitle("Disturbance Frequency 2 Samples") +
  guides(color="none", fill=guide_legend("Genus")) +
  facet_grid(~Substrate, scales = "free")

p_f2


ggsave("Supplemental_Figure_1_freq2_abundance.tiff", height = 7.5, width = 10, units = "in", device = "tiff", dpi = 300)
#ggsave("freq2_abundance.png", device = "png", dpi = 300)
#ggsave("freq2_abundance.pdf", device = "pdf", dpi = 300)

## Get summary table of ASVs for easy reference
data2 <- data

#by treatment+taxa
summary1 <- data2 %>%
  group_by(Substrate, variable) %>%
  summarise(mean=mean(value)*100,
            sd=sd(value)*100,
            count=n())

#write.csv(summary1, "taxa_summary1.csv", row.names=FALSE)


#by disturbance+taxa
summary2 <- data2 %>%
  group_by(Substrate, Frequency, variable) %>%
  summarise(mean=mean(value*100),
            sd=sd(value*100),
            count=n())

#by disturbance+day+taxa
summary3 <- data2 %>%
  group_by(Substrate, Time, Frequency, variable) %>%
  summarise(mean=mean(value*100),
            sd=sd(value*100),
            count=n())

summary <- metadata %>%
  group_by(Substrate,Frequency, Time) %>%
  summarise(count=n())


###### LINE GRAPH FROM FIRST SUBMISSION, IGNORE IN FAVOR OF BAR GRAPH ###############

datag2 <- datag %>%
  group_by(variable, Frequency, Time) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            count=n())

datac2 <- datac %>%
  group_by(variable, Frequency, Time) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            count=n())

#make list of abundant+relevant taxa for each substrate treatment
tax_list <- c("Pseudomonas.1","Asticcacaulis","Cellvibrio","Pseudomonas.2","Cohnella","Lacunisphaera")



#geom_text(data=sums, aes(x=Substrate, y= percent, label=sum)) +


datag2 <- subset(datag2, variable %in% tax_list)
datac2 <- subset(datac2, variable %in% tax_list)

sumsg <- datag2[,c(3,6)]
sumsg$y <- 90

pgline <- ggplot(datag2,aes(y=mean, x= Time, color=variable)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  geom_text(data=datag2, aes(x= Time, y = 90, label = count)) +
  theme_classic(base_size=15) + #make it pretty
  theme(#axis.text.x=element_blank(),
    legend.position="bottom") + 
  #theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
  ylab("Relative Abundance (%)") + 
  xlab("Time (Day)") +
  labs(fill='Genus')+ #change the legend title
  scale_color_manual(values=getPalette(colorCount)) +
  facet_grid( . ~ Frequency, scales = "free_x") + 
  ggtitle("Glucose Samples faceted by Disturbance Frequency (1/n Days)") +
  guides(color=guide_legend("Genus"))


pgline

pcline <- ggplot(datac2,aes(y=mean, x= Time, color=variable)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
  geom_text(data=datac2, aes(x= Time, y = 90, label = count)) +
  theme_classic(base_size=15) + #make it pretty
  theme(#axis.text.x=element_blank(),
    legend.position="bottom") + 
  #theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
  ylab("Relative Abundance (%)") + 
  xlab("Time (Day)") +
  labs(fill='Genus')+ #change the legend title
  scale_color_manual(values=getPalette(colorCount)) +
  facet_grid( . ~ Frequency, scales = "free_x") + 
  ggtitle("Cellulose Samples faceted by Disturbance Frequency (1/n Days)") +
  guides(color=guide_legend("Genus"))

pcline


ggarrange(pcline + rremove("xlab"), pgline, common.legend = TRUE, ncol=1, nrow=2,legend = "bottom")



### OLD ALPHA DIVERSITY ANALYSIS, Can be ignored ###



#simple counts. Ignore this first one, not as useful as Menhinick's
#count <- ddply(data,~Sample_Name,function(x){
#  data.frame(Count=sum(x[,8:173]>0))
#})
#count_list<-list(metadata,richness)


#Menhinick's index - 'species' richness
# number of species divided by square root of number of individuals in sample
menhinick <- function(x) {
  sum(x>0)/sqrt(sum(x))
}


richness <- ddply(data,~Sample_Name,function(x){
  data.frame(Richness=menhinick(x[,8:173]>0))
})


# Shannon
#Assumes all species are represented and randomly sampled, balance with second
shan <- ddply(data,~Sample_Name,function(x){
  data.frame(Shannon=diversity(x[,8:173],index="shannon"))
})
shan
#shan_list<-list(metadata,shan)


# Inverse Simpon
# Provides more weight to common or dominant species, a few rare species will not affect. balance with second index
invsimpson <- ddply(data,~Sample_Name,function(x){
  data.frame(InverseSimpson=diversity(x[,8:173],index="invsimpson"))
})
invsimpson
#in_list<-list(metadata,invsimpson)

#Pieloue's Index- evennnes
# Represents ratio between overserved Shannon, and greateast possible shannon for that sample
pielou <-ddply(data,~Sample_Name,function(x){
  data.frame(Pielou=exp(diversity(x[,8:173], index="shannon"))/sum(x[,8:173]>0))
})
pielou
#pie_list<-list(metadata,pielou)

metadata <- merge(metadata, richness, by="Sample_Name")
metadata <- merge(metadata, shan, by="Sample_Name")
metadata <- merge(metadata, invsimpson, by="Sample_Name")
metadata <- merge(metadata, pielou, by="Sample_Name")
metadata <- merge(metadata, count, by="Sample_Name")

#write.csv(metadata, "disturbance_metadata_diversity_indices_no2CDay5.csv", row.names = FALSE)


###### Graphs######

#metadata$Frequency <- as.factor(metadata$Frequency)

#my_comparisons <- list(c("3","7"),c("3","5"),c("3","2"),c("3","1"), c("1","2"),c("1","3"),c("1","5"),c("1","7"))


#I should just write a function lol oh well

pshan <-ggboxplot(metadata, x="Frequency",y="Shannon", color="Substrate", fill="Substrate", lwd=1) +
  theme_classic(base_size=22)+
  theme(aspect.ratio=1) +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  labs(y="Shannon Index", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Shannon Diversity") + 
  stat_compare_means(method="t.test", aes(group=Substrate), label="p.signif")
pshan

pinv<-ggboxplot(metadata, x="Frequency",y="InverseSimpson", color="Substrate", fill="Substrate", lwd=1) +
  theme_classic(base_size=22)+
  theme(aspect.ratio=1) +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  labs(y="Inverse Simpson Index", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Inverse Simpson Diversity") + 
  stat_compare_means(method="t.test", aes(group=Substrate), label="p.signif")

ppie<-ggboxplot(metadata, x="Frequency",y="Pielou", color="Substrate", fill="Substrate", lwd=1) +
  theme_classic(base_size=22)+
  theme(aspect.ratio=1) +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  labs(y="Pielou's Index", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Pielou's Evenness") + 
  stat_compare_means(method="t.test", aes(group=Substrate), label="p.signif")

prich <-ggboxplot(metadata, x="Frequency",y="Richness", color="Substrate", fill="Substrate", lwd=1) +
  theme_classic(base_size=22)+
  theme(aspect.ratio=1) +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  labs(y="Menhinick's Diversity Index", x="Disturbance Frequency (1/n days)") + 
  ggtitle("Richness") + 
  stat_compare_means(method="t.test", aes(group=Substrate), label="p.signif")

pshan
pinv
ppie
prich


#ggarrange(pshan + rremove("xlab"), pinv + rremove("xlab"), ppie,prich, common.legend = TRUE, ncol=2, nrow=2,legend = "right")

ggarrange(pshan, pinv, common.legend = TRUE, ncol=2,legend = "right")

ggarrange(ppie, prich, common.legend = TRUE, ncol=2,legend = "right")




#stats 
datacel <- subset(metadata, Substrate == "Cellulose")
dataglu <- subset(metadata, Substrate == "Glucose")

#run a bunch of t tests
tab1<-compare_means(Shannon ~ Frequency, data=datacel, ref.group="3", method="t.test")
tab2<-compare_means(InverseSimpson ~ Frequency, data=datacel, ref.group="3", method="t.test")
tab3<-compare_means(Shannon ~ Frequency, data=dataglu, ref.group="1", method="t.test")
tab4<-compare_means(InverseSimpson ~ Frequency, data=dataglu, ref.group="1", method="t.test")

tab_total <- rbind(tab1,tab2,tab3,tab4)
colnames(tab_total)[1] = "DiversityIndex"
colnames(tab_total)[2] = "Reference"
colnames(tab_total)[3] = "Comparison"

#write.csv(tab_total, "disturbance_diversity_stats.csv", row.names = FALSE)
