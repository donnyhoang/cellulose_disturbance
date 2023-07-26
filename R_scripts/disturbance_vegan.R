library(vegan)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
#references
#https://rpubs.com/an-bui/vegan-cheat-sheet
#https://sites.google.com/site/mb3gustame/home
#https://www.flutterbys.com.au/stats/tut/tut13.2.html

metadata <- read.csv("disturbance_total_metadata_no2CDay5.csv", header=TRUE)
otu <- read.csv("disturbance_total_otu_genus_no2CDay5.csv", header=TRUE, check.names = FALSE)

metadata$Frequency <- as.character(metadata$Frequency)
metadata$Frequency <- factor(metadata$Frequency, levels=c("7","5","3","2","1"))
metadata$Substrate <- factor(metadata$Substrate, levels=c("Cellulose","Glucose"))

#attach otu to metadata
data <- cbind(metadata,otu)
#drop columns
data <-data[-c(8)]
data$Frequency <- as.character(data$Frequency)

### ALPHA DIVERSITY ###

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

###### Abundance Bar graph ######

metadata <- read.csv("disturbance_total_metadata_top15_no2CDay5.csv", header=TRUE)
otu <- read.csv("disturbance_total_otu_genus_top15_no2CDay5.csv", header=TRUE, check.names = FALSE)



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


datac <- subset(data, data$Substrate=="Cellulose")
datag <- subset(data, data$Substrate=="Glucose")

#create list of new label names and label function for graphs
#freq.labs <- c("Frequency 1/7","Frequency 1/5","Frequency 1/3","Frequency 1/2","Frequency 1/1")
#names(freq.labs) <- c("7","5","3","2","1")

#namesg <- datag$Time_Replicate
#namesc <- datac$Time_Replicate

#Can't get xaxis write in R, bc I'm essentially implying a continuous xaxis when I'm plotting discretely.
#Also, number of replicates within timepoints samples are different

#Make dummy labels, then MANUALLY EDIT XAXIS OUTSIDE OF R 
#Probably did this in Paint, please don't clown on me
names <- c("1"," "," "," "," ","2"," "," "," "," ","3"," "," "," "," ","4"," "," "," "," ",
           "5"," "," "," "," ","6"," "," "," "," ","7"," "," "," "," ")

## Bar graphs (these ended up in supplemental)
pg <-ggplot(datag) + 
  geom_bar(aes(y=value,x=Sample_Name,fill=variable,color=variable),stat="identity") + 
  theme_classic(base_size=15) + #make it pretty
  theme(#axis.text.x=element_blank(),
    legend.position="bottom") + 
  #theme(axis.text.x=element_text(angle=45, vjust=0.5, hjust=1)) +
  ylab("Relative Abundance (%)") + 
  xlab("Samples arranged by Time (Day Sampled) in replicates") +
  labs(fill='Genus')+ #change the legend title
  scale_fill_manual(values=getPalette(colorCount)) +
  scale_color_manual(values=getPalette(colorCount)) +
  scale_x_discrete(labels=names) +
  facet_grid( . ~ Frequency, scales = "free_x") + 
  ggtitle("Glucose Samples faceted by Disturbance Frequency (1/n Days)") +
  guides(color="none", fill=guide_legend("Genus"))

pc <-ggplot(datac) + 
  geom_bar(aes(y=value,x=Sample_Name,fill=variable,color=variable),stat="identity") + 
  theme_classic(base_size=15) + #make it pretty
  theme(#axis.text.x=element_blank(),
    legend.position="bottom") + 
  #theme(axis.text.x=element_text(angle=45, vjust=0.5)) +
  ylab("Relative Abundance (%)") + 
  xlab("Samples arranged by Time (Day Sampled)in replicates") +
  labs(fill='Genus')+ #change the legend title
  scale_fill_manual(values=getPalette(colorCount)) +
  scale_color_manual(values=getPalette(colorCount)) +
  scale_x_discrete(labels=names) +
  facet_grid( . ~ Frequency , scales = "free_x")  + 
  ggtitle("Cellulose Samples faceted by Disturbance Frequency (1/n Days)") +
  guides(color="none", fill=guide_legend("Genus"))


ggarrange(pc + rremove("xlab"), pg, common.legend = TRUE, ncol=1, nrow=2,legend = "bottom")


###### line graph (use this for main figure since there's less going on)

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


datag2 <- subset(datag2, variable %in% tax_list)
datac2 <- subset(datac2, variable %in% tax_list)

pgline <- ggplot(datag2,aes(y=mean, x= Time, color=variable)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
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




pcline <- ggplot(datac2,aes(y=mean, x= Time, color=variable)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.2) +
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





#For looking at diversity, facetted by days. Ended up not using, can ignore.
shan_plotter <- function(...){
  return(...%>%
           ggboxplot(x="Frequency",y="Shannon", color="Substrate", fill="Substrate", lwd=1) +
           theme_classic(base_size=12)+
           theme(aspect.ratio=1) +
           scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
           scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
           labs(y="Shannon Index", x="Disturbance Frequency (1/n days)") + 
           ggtitle("Shannon Diversity") +
           stat_compare_means(method="t.test", aes(group=Substrate), label="p.signif"))
}

isimp_plotter <- function(...){
  return(...%>%
           ggboxplot(x="Frequency",y="InverseSimpson", color="Substrate", fill="Substrate", lwd=1) +
           theme_classic(base_size=12)+
           theme(aspect.ratio=1) +
           scale_fill_manual(values=c("#e1f8df","#fff1f9")) +
           scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
           labs(y="Inverse Simpson Index", x="Disturbance Frequency (1/n days)") + 
           ggtitle("Inverse Simpson") +
           stat_compare_means(method="t.test", aes(group=Substrate), label="p.signif"))
}

pil_plotter <- function(...){
  return(...%>%
           ggboxplot(x="Frequency",y="Pielou", color="Substrate", fill="Substrate", lwd=1) +
           theme_classic(base_size=12)+
           theme(aspect.ratio=1) +
           scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
           scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
           labs(y="Pielou's Evenness", x="Disturbance Frequency (1/n days)") + 
           ggtitle("Pielou's Index") +
           stat_compare_means(method="t.test", aes(group=Substrate), label="p.signif"))
}

richness_plotter <- function(...){
  return(...%>%
           ggboxplot(x="Frequency",y="Richness", color="Substrate", fill="Substrate", lwd=1) +
           theme_classic(base_size=12)+
           theme(aspect.ratio=1) +
           scale_fill_manual(values=c("#e1f8df","#fff1f9")) + 
           scale_color_manual(values=c("#7fbf7b","#e9a3c9")) +
           labs(y="Count of Unique ASVs", x="Disturbance Frequency (1/n days)") + 
           ggtitle("Observed Richness") +
           stat_compare_means(method="t.test", aes(group=Substrate), label="p.signif"))
}



### Diversity by day ####



day1 <- subset(metadata, metadata$Time=="1")
ps1 <- shan_plotter(day1) + labs(title = "Day 1")
pi1 <- isimp_plotter(day1) + labs(title = "Day 1")
pp1 <- pil_plotter(day1) + labs(title = "Day 1")
pr1 <- richness_plotter(day1) + labs(title = "Day 1")
#ggarrange(ps1 + rremove("xlab"), pi1 + rremove("xlab"), pp1,pr1, common.legend = TRUE, ncol=2, nrow=2,legend = "right")

day2 <- subset(metadata, metadata$Time=="2")
ps2 <- shan_plotter(day2) + labs(title = "Day 2")
pi2 <- isimp_plotter(day2) + labs(title = "Day 2")
pp2 <- pil_plotter(day2) + labs(title = "Day 2")
pr2 <- richness_plotter(day2) + labs(title = "Day 2")
#ggarrange(ps2 + rremove("xlab"), pi2 + rremove("xlab"), pp2,pr2, common.legend = TRUE, ncol=2, nrow=2,legend = "right")

day3 <- subset(metadata, metadata$Time=="3")
ps3 <- shan_plotter(day3) + labs(title = "Day 3")
pi3 <- isimp_plotter(day3)  + labs(title = "Day 3")
pp3 <- pil_plotter(day3)  + labs(title = "Day 3")
pr3 <- richness_plotter(day3)  + labs(title = "Day 3")
#ggarrange(ps3 + rremove("xlab"), pi3 + rremove("xlab"), pp3,pr3, common.legend = TRUE, ncol=2, nrow=2,legend = "right")

day4 <- subset(metadata, metadata$Time=="4")
ps4 <- shan_plotter(day4)  + labs(title = "Day 4")
pi4 <- isimp_plotter(day4)  + labs(title = "Day 4")
pp4 <- pil_plotter(day4)  + labs(title = "Day 4")
pr4 <- richness_plotter(day4)  + labs(title = "Day 4")
#ggarrange(ps4 + rremove("xlab"), pi4 + rremove("xlab"), pp4,pr4, common.legend = TRUE, ncol=2, nrow=2,legend = "right")

day5 <- subset(metadata, metadata$Time=="5")
ps5 <- shan_plotter(day5)  + labs(title = "Day 5")
pi5 <- isimp_plotter(day5)  + labs(title = "Day 5")
pp5 <- pil_plotter(day5)  + labs(title = "Day 5")
pr5 <- richness_plotter(day5)  + labs(title = "Day 5")
#ggarrange(ps5 + rremove("xlab"), pi5 + rremove("xlab"), pp5,pr5, common.legend = TRUE, ncol=2, nrow=2,legend = "right")

day6 <- subset(metadata, metadata$Time=="6")
ps6 <- shan_plotter(day6)  + labs(title = "Day 6")
pi6 <- isimp_plotter(day6)  + labs(title = "Day 6")
pp6 <- pil_plotter(day6)  + labs(title = "Day 6")
pr6 <- richness_plotter(day6)  + labs(title = "Day 6")
#ggarrange(ps6 + rremove("xlab"), pi6 + rremove("xlab"), pp6,pr6, common.legend = TRUE, ncol=2, nrow=2,legend = "right")

day7 <- subset(metadata, metadata$Time=="7")
ps7 <- shan_plotter(day7)  + labs(title = "Day 7")
pi7 <- isimp_plotter(day7)  + labs(title = "Day 7")
pp7 <- pil_plotter(day7)  + labs(title = "Day 7")
pr7 <- richness_plotter(day7)  + labs(title = "Day 7")
#ggarrange(ps7 + rremove("xlab"), pi7 + rremove("xlab"), pp7,pr7, common.legend = TRUE, ncol=2, nrow=2,legend = "right")

ggarrange(ps1+ rremove("xlab"),
          ps2+ rremove("xlab") + rremove("ylab"),
          ps3+ rremove("xlab") + rremove("ylab"),
          ps4+ rremove("xlab"),
          ps5+ rremove("ylab"),
          ps6+ rremove("ylab"),
          ps7, 
          common.legend = TRUE, ncol=3, nrow=3,legend = "bottom")

ggarrange(pi1+ rremove("xlab"),
          pi2+ rremove("xlab") + rremove("ylab"),
          pi3+ rremove("xlab") + rremove("ylab"),
          pi4+ rremove("xlab"),
          pi5+ rremove("ylab"),
          pi6+ rremove("ylab"),
          pi7, 
          common.legend = TRUE, ncol=3, nrow=3,legend = "bottom")

ggarrange(pp1+ rremove("xlab"),
          pp2+ rremove("xlab") + rremove("ylab"),
          pp3+ rremove("xlab") + rremove("ylab"),
          pp4+ rremove("xlab"),
          pp5+ rremove("ylab"),
          pp6+ rremove("ylab"),
          pp7, 
          common.legend = TRUE, ncol=3, nrow=3,legend = "bottom")

ggarrange(pr1+ rremove("xlab"),
          pr2+ rremove("xlab") + rremove("ylab"),
          pr3+ rremove("xlab") + rremove("ylab"),
          pr4+ rremove("xlab"),
          pr5+ rremove("ylab"),
          pr6+ rremove("ylab"),
          pr7, 
          common.legend = TRUE, ncol=3, nrow=3,legend = "bottom")

