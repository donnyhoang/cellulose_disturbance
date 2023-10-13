library(BiodiversityR) #this loads vegan too
library(ggplot2)
#library(rgl) #for 3D plot
library(plot3D)
library(scatterplot3d)
library(usedist)
library(stringr)
library(purrr)
library(ggpubr)
library(tidyverse)
library(cowplot)

#https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html
#http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization
#https://chrischizinski.github.io/rstats/adonis/

metadata <- read.csv("disturbance_total_metadata_no2CDay5.csv", header=TRUE)
otu <- read.csv("disturbance_total_otu_genus_no2CDay5.csv", header=TRUE, row.names=1, check.names = FALSE)

#restructure factors for easier graphing+visualization
metadata$Frequency <- as.character(metadata$Frequency)
metadata$Frequency <- factor(metadata$Frequency, levels=c("7","5","3","2","1"))
metadata$Substrate <- factor(metadata$Substrate, levels=c("Cellulose","Glucose"))

attach(metadata)

#ordinate in 3 dimensions bc it cannot find solution with 2
ordination.model <- metaMDS(otu, distance='bray', k=3)

#extract coordinates so I can plot with ggplot2
nmds_coords <- ordination.model$points

#attach to metadata
data <- cbind(metadata, nmds_coords)
#write.csv(data, "disturbance_ordination_coords_no2CDay5.csv", row.names=FALSE)

####2d nmds plot#######
plot <- ggplot() +
  geom_point(data=data,
             aes(x=MDS1, y=MDS2, color=Substrate, bg=Substrate, shape=Frequency),
             size=4) +
  theme_classic() +
  xlab("NMDS1") +
  ylab("NMDS2") +
  coord_fixed(ratio=1) +
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  scale_fill_manual(values=c("#7fbf7b","#e9a3c9")) +
  scale_shape_manual(values=c(25,24,23,22,21))
plot

#color by freq. Coloring by Freq easier to read imo, keep this one
plot2 <- ggplot() +
  geom_point(data=data,
             aes(x=MDS1, y=MDS2, color=Frequency, bg=Frequency, shape=Substrate),
             size=3) +
  theme_classic(base_size=15) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  #coord_fixed(ratio=1) +
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854")) +
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854")) +
  scale_shape_manual(values=c(21,22)) +
  guides(fill=guide_legend(ncol=2)) +
  ggtitle("Biplot of Bray-Curtis Distance Matrix")
plot2

######## stats ###########
#https://chrischizinski.github.io/rstats/adonis/

dist <- vegdist(otu) #create distance matrix
attach(metadata) #attach metadata


#anosim
#can i reject null hypothesis that there are no differences between groups?
anoSub <- anosim(dist, Substrate) #similarity analysis
summary(anoSub) # see summary
plot(anoSub)
anoFre <- anosim(dist, Frequency)
anoRep <-anosim(dist, Replicate)
anoTime <- anosim(dist, Time)
anoRun<- anosim(dist, Run)

summary(anoSub)
summary(anoFre)
summary(anoRep)
summary(anoTime)
summary(anoRun)


#adonis
#null hypothesis: groups do not different in spread or position in multivariate space
adSubstrate <- adonis2(dist ~ Substrate, data=metadata, permutations=99, method="bray")
adFrequency <- adonis2(dist ~ Frequency, data=metadata,  permutations=99, method="bray")
adReplicate <- adonis2(dist ~ Replicate, data=metadata,  permutations=99, method="bray")
adTime <- adonis2(dist ~ Time, data=metadata,  permutations=99, method="bray")
adRun <- adonis2(dist ~ Run, data=metadata,  permutations=99, method="bray")
adSubstrate
adonis2(dist ~ Substrate, data=metadata, permutations=99, method="bray")

adSubFre <- adonis2(dist ~ Substrate * Frequency, data=metadata,  permutations=99, method="bray")
adTiFre <- adonis2(dist ~ Time * Frequency, data=metadata,  permutations=99, method="bray")



Variables <- c("Substrate", "Frequency", "Replicate", "Time", "Illumina Run")
ANOSIM_R <- c(0.912,0.083,-0.009,0.008,0.041)
ANOSIM_p <- c(0.001,0.001,0.993,0.165,0.001)
PERMANOVA <- c(0.624,0.115,0.003,0.11,0.023)
PERMANOVA_p <- c(0.010,0.010,1,0.040,0.010)

stat_table <- data.frame(Variables, ANOSIM_R, ANOSIM_p, PERMANOVA, PERMANOVA_p)
colnames(stat_table) <- c("Variable", "ANOSIM", "ANOSIM \n p-value", "PERMANOVA", "PERMANOVA \n p-value")

p_stat_table <- ggtexttable(stat_table, rows = NULL,
                            theme = ttheme("classic"))
p_stat_table


p1 <- plot2 + annotation_custom(ggplotGrob(p_stat_table),xmin = 0.25, ymin= 0.75, xmax= 1)
p1

######### centroid calculations #############
dist <- vegdist(otu) #create distance matrix
attach(metadata) #attach metadata


#distance between centroids (how much to disturbance frequencies vary between cellulose and glucose?) ###

start <- c("1C","1C","1C","1C","2C","2C","2C","3C","3C","5C")
end <- c("2C","3C","5C","7C","3C","5C","7C","5C","7C","7C")

dbc_cellulose <- map2_dfr(start, end, ~ {
  idx1 <- which(metadata$Substrate_Frequency == .x)
  idx2 <- which(metadata$Substrate_Frequency == .y)
  tibble(
    start_centroid = .x,
    end_centroid = .y,
    distance = dist_between_centroids(dist, idx1, idx2)
  )
})
dbc_cellulose$Substrate <- "Cellulose"



start2 <- c("1G","1G","1G","1G","2G","2G","2G","3G","3G","5G")
end2 <- c("2G","3G","5G","7G","3G","5G","7G","5G","7G","7G")

dbc_glucose <- map2_dfr(start2, end2, ~ {
  idx1 <- which(metadata$Substrate_Frequency == .x)
  idx2 <- which(metadata$Substrate_Frequency == .y)
  tibble(
    start_centroid = .x,
    end_centroid = .y,
    distance = dist_between_centroids(dist, idx1, idx2)
  )
})
dbc_glucose$Substrate <- "Glucose"


dbc <- rbind(dbc_cellulose,dbc_glucose)

dbc_mean <- dbc %>%
  group_by(Substrate) %>%
  summarise(mean = mean(distance))


###distance to centroid (how much variation is there within disturbance frequency groupings, between substrates?) #######

dtc <- dist_to_centroids(dist, Substrate_Frequency)

###you calculated all samples to all centroids. keep only the ones that match their respective grouping.
#str_split to create variable to filter with
dtc[c("Substrate_Frequency", "Replicate","Day")] <- str_split_fixed(dtc$Item, "_", 3)
dtc <- dtc[which(dtc[,2] == dtc[,4]), ]
dtc[c("Frequency", "Substrate")] <- str_split_fixed(dtc$Substrate_Frequency, "", 2)
dtc <- dtc %>%
  mutate(across("Substrate", str_replace, "C", "Cellulose"))
dtc <- dtc %>%
  mutate(across("Substrate", str_replace, "G", "Glucose"))

dtc$Frequency <- factor(dtc$Frequency, levels = c("7","5","3","2","1"))


#plot centroid calculations


stats_data <- function(y) {
  return(data.frame(
    y=0.85,
    label = paste('n =', length(y), '\n')
  ))
}


p_dbc <- ggplot(dbc, aes(x = Substrate, y = distance, color = Substrate, fill = Substrate)) +
  geom_boxplot(lwd = 1, outlier.shape = NA) +
  geom_jitter(size = 2, width = 0.15) +
  theme_classic(base_size = 15) +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) +
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  stat_compare_means(method = "t.test", label = "p.signif", label.x = c(1.5))+
  stat_summary(fun.data = stats_data,
               geom = "text",
               position = position_dodge(width = 0.75)) +
  labs(y="Distance") +
  ggtitle("Distance Between Centroids")
p_dbc



stats_data2 <- function(y) {
  return(data.frame(
    y=0.57,
    label = paste('n =', length(y), '\n')
  ))
}


p_dtc <- ggplot(dtc, aes(x = Frequency, y = CentroidDistance, color = Substrate, fill = Substrate)) +
  geom_boxplot(lwd = 1, outlier.size = 2) +
  #geom_jitter(width = 0.15) +
  theme_classic(base_size = 15) +
  #theme(legend.position="none") +
  scale_fill_manual(values=c("#e1f8df","#fff1f9")) +
  scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
  #facet_grid(~Substrate, scales = "free") +
  scale_x_discrete(labels = c("1/7","1/5","1/3","1/2","1/1"))+
  stat_compare_means(method = "t.test", label = "p.signif", label.x = c(1.5))+
  stat_summary(fun.data = stats_data2,
               geom = "text",
               position = position_dodge(width = 0.75)) +
  labs(y="Distance", x = "Disturbance Frequency") +
  ggtitle("Distance to Centroids")
p_dtc





p1 <- plot2 
p2 <- p_stat_table
p3 <- p_dbc +theme(legend.position="none") 
p4 <- p_dtc

fig3 <- ggdraw () +
  draw_plot(p1, x = 0, y = .4, width = .6, height = .6) +
  draw_plot(p2, x = 0.7, y = .8, width = .2, height = .2) +
  draw_plot(p3, x = 0.6, y = 0.4, width = 0.4, height = 0.4) +
  draw_plot(p4, x = 0, y = 0, width = 1, height = 0.4) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0.6, 0.6, 0), y = c(1, 1, 0.8, 0.4))
fig3

ggsave("Figure_3.tiff", device = "tiff", dpi = 700)


####### ######## 3D plot ########################
mycolors <- c("#7fbf7b","#e9a3c9")
data$color <- mycolors[as.numeric(data$Substrate)]

##subset data to add to graph
data1 <- subset(data, data$Frequency=="1")
data2 <- subset(data, data$Frequency=="2")
data3 <- subset(data, data$Frequency=="3")
data5 <- subset(data, data$Frequency=="5")
data7 <- subset(data, data$Frequency=="7")

##plot, add

open3d()
pch3d(x=data1$MDS1, y=data1$MDS2, z=data1$MDS3, cex=0.7, bg = data1$color, pch=21)
pch3d(x=data2$MDS1, y=data2$MDS2, z=data2$MDS3, cex=0.7, bg = data2$color, pch=22)
pch3d(x=data3$MDS1, y=data3$MDS2, z=data3$MDS3, cex=0.7, bg = data3$color, pch=23)
pch3d(x=data5$MDS1, y=data5$MDS2, z=data5$MDS3, cex=0.7, bg = data5$color, pch=24)
pch3d(x=data7$MDS1, y=data7$MDS2, z=data7$MDS3, cex=0.7, bg = data7$color, pch=25)
axes3d()

title3d(xlab="NMDS1", ylab="NMDS2", zlab="NMDS3")


legend3d( "topright", title="Frequency", c("7", "5", "3","2","1"), 
          pch = c(25,24,23,22,21))

k <- sort(unique(data$Substrate))
legend3d("right", legend=k, pch=18, col=mycolors, title="Substrate")


##3d plot, color by frequency


mycolors2 <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854")
data$colors <- mycolors2[as.numeric(data$Frequency)]

datac <- subset(data, data$Substrate=="Cellulose")
datag <- subset(data, data$Substrate=="Glucose")


open3d()
par3d(windowRect = c(100, 100, 612, 612))
pch3d(x=datac$MDS1, y=datac$MDS2, z=datac$MDS3, cex=0.4, bg = datac$color, pch=21)
pch3d(x=datag$MDS1, y=datag$MDS2, z=datag$MDS3, cex=0.7, bg = datag$color, pch=22)
axes3d()

title3d(xlab="NMDS1", ylab="NMDS2", zlab="NMDS3")

l <- sort(unique(data$Frequency))
legend3d("right", legend=l, pch=18, col=mycolors2, title="Frequency")
legend3d("topright", title="Substrate", c("Cellulose", "Glucose"), 
          pch = c(21,22))

rgl.postscript("plot.pdf",fmt="pdf")

####### 3D plot #############
source('~/hubiC/Documents/R/function/addgrids3d.r')


colors <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854")
data$colors <- colors[as.numeric(data$Frequency)]

shapes <- c(19, 15)
data$shapes <- shapes[as.numeric(data$Substrate)]

scatterplot3d(data[,8:10], pch = data$shapes, color = data$colors, box = TRUE, angle = 60, xlab="NMDS1", ylab="NMDS2", zlab="NMDS3")
legend("topleft", legend = levels(data$Frequency), col = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854"), pch = 16, title = "Frequency")
legend("bottomright", legend = levels(data$Substrate), pch = c(19.15), xpd = TRUE, title = "Substrate")


tiff('Supplemental_Figure_4.tiff', units="in", width=8, height=9, res=1200, compression = 'lzw')

scatterplot3d(data[,8:10], pch = data$shapes, color = data$colors, box = TRUE, angle = 60, xlab="NMDS1", ylab="NMDS2", zlab="NMDS3")
legend("topleft", legend = levels(data$Frequency), col = c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854"), pch = 16, title = "Frequency")
legend("bottomright", legend = levels(data$Substrate), pch = c(19.15), xpd = TRUE, title = "Substrate")



dev.off()
       
