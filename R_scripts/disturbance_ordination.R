library(BiodiversityR) #this loads vegan too
library(ggplot2)
library(rgl) #for 3D plot

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
             size=4) +
  theme_classic(base_size=22) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  coord_fixed(ratio=1) +
  scale_color_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854")) +
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854")) +
  scale_shape_manual(values=c(21,22))
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

adonis2(dist ~ Substrate, data=metadata, permutations=99, method="bray")

adSubFre <- adonis2(dist ~ Substrate * Frequency, data=metadata,  permutations=99, method="bray")
adTiFre <- adonis2(dist ~ Time * Frequency, data=metadata,  permutations=99, method="bray")


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
pch3d(x=datac$MDS1, y=datac$MDS2, z=datac$MDS3, cex=0.4, bg = datac$color, pch=21)
pch3d(x=datag$MDS1, y=datag$MDS2, z=datag$MDS3, cex=0.7, bg = datag$color, pch=22)
axes3d()

title3d(xlab="NMDS1", ylab="NMDS2", zlab="NMDS3", main = "NMDS plot of community composition")

l <- sort(unique(data$Frequency))
legend3d("right", legend=l, pch=18, col=mycolors2, title="Frequency")
legend3d("topright", title="Substrate", c("Cellulose", "Glucose"), 
          pch = c(21,22))








