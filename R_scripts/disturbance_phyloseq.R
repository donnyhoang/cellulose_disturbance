library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(QsRutils)
library(dplyr)
#references
#http://www.hiercourse.com/docs/microbial/02_alphaDiversity.html
#https://rpubs.com/an-bui/vegan-cheat-sheet

#load data, have two different outputs from DADA2 so read them in individually
tax1 <- read.csv("disturbance_run1_taxonomy.csv", header=TRUE, row.names=1)
counts1 <- read.csv("disturbance_run1_counts.csv", header=TRUE, row.names = 1, check.names = FALSE)
metadata1 <- read.csv("disturbance_metadata_run1.csv", header=TRUE, row.names= 1)
#change things to character (not numeric)
metadata1$Frequency <- as.character(metadata1$Frequency)
#metadata$Time <- as.character(metadata$Time)

tax2 <- read.csv("disturbance_run2_taxonomy.csv", header=TRUE, row.names=1)
counts2 <- read.csv("disturbance_run2_counts.csv", header=TRUE, row.names = 1, check.names = FALSE)
metadata2 <- read.csv("disturbance_metadata_run2.csv", header=TRUE, row.names = 1)
#change things to character (not numeric)
metadata2$Frequency <- as.character(metadata2$Frequency)
#metadata$Time <- as.character(metadata$Time)



#get sequencing summary, not part of phyloseq wrangling, but just get summary table of sequencing stats while files are loaded
seq_summary1 <- as.data.frame(colSums(counts1 !=0))
seq_summary2 <- as.data.frame(colSums(counts2 !=0))
colnames(seq_summary1)[1] <- "count"
colnames(seq_summary2)[1] <- "count"
seq_summary <- rbind(seq_summary1, seq_summary2)
seq_summary %>%
  summarise(mean=mean(count),
            sd=sd(count))

reads1 <- read.csv("reads_tracked_run1.csv", header=TRUE, row.names =1)
reads2 <- read.csv("reads_tracked_run2.csv", header=TRUE, row.names =1)
reads <- rbind(reads1, reads2)
rows_remove <- c("sum", "sum1")
reads <- reads[!(row.names(reads) %in% rows_remove),]
reads %>%
  summarise(mean=mean(nonchim),
            sd=sd(nonchim))


#create phyloseq object
count.ps1 <- otu_table(as.matrix(counts1),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps1 <- tax_table(as.matrix(tax1)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps1 <- sample_data(metadata1) #define our metadata as sample information
ps1 <- phyloseq(count.ps1, tax.ps1, meta.ps1) #put them all together into a phyloseq object

#create second phyloseq ob
count.ps2 <- otu_table(as.matrix(counts2),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps2 <- tax_table(as.matrix(tax2)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps2 <- sample_data(metadata2) #define our metadata as sample information
ps2 <- phyloseq(count.ps2, tax.ps2, meta.ps2) #put them all together into a phyloseq object



#subset samples, taking only experimental samples and anything over 100- reads
ps1 <- subset_samples(ps1, Frequency %in% c("1","2","3","5","7"))
ps1 <- subset_samples(ps1, sample_sums(ps1) >=1000)


##take top 20 abundant, transform to proportion
#actually, take 14 since it reduces into something comfortable for Set3 color palette later
top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:14] 
ps1 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps1 <- prune_taxa(top20, ps1)

#write otu_table were asv is replaced with genus 
#https://github.com/joey711/phyloseq/issues/1559
ps1 <- tax_glom(ps1, "Genus")
otu_genus <- otu_table(ps1)
otu_genus <- t(otu_genus)
df_genus1 <- as.data.frame(otu_genus)
colnames(df_genus1) <- as.data.frame(tax_table(ps1))$Genus


#############repeat for second set
#subset 
ps2 <- subset_samples(ps2, Frequency %in% c("1","2","3","5","7"))
ps2 <- subset_samples(ps2, sample_sums(ps2) >=1000)


##take top 20 abundant, transform to proportion
top20 <- names(sort(taxa_sums(ps2), decreasing=TRUE))[1:14] 
ps2 <- transform_sample_counts(ps2, function(OTU) OTU/sum(OTU))
ps2 <- prune_taxa(top20, ps2)

#write otu_table were asv is replaced with genus 
#https://github.com/joey711/phyloseq/issues/1559
ps2 <- tax_glom(ps2, "Genus")
otu_genus <- otu_table(ps2)
otu_genus <- t(otu_genus)
df_genus2 <- as.data.frame(otu_genus)
colnames(df_genus2) <- as.data.frame(tax_table(ps2))$Genus



#########merge data 
newmeta1 <- vegansam(ps1)
newmeta2 <- vegansam(ps2)
metadata <- rbind(newmeta1, newmeta2)
otu <- bind_rows(df_genus1,df_genus2)
otu[is.na(otu)] <- 0

#write data, including separate files for datasets 1 and 2
write.csv(newmeta1, "disturbance_run1_metadata_top15.csv", row.names=TRUE)
write.csv(df_genus1, "disturbance_run1_otu_genus_top15.csv", row.names=TRUE)

write.csv(newmeta2, "disturbance_run2_metadata_top15.csv", row.names=TRUE)
write.csv(df_genus2, "disturbance_run2_otu_genus_top15.csv", row.names=TRUE)

write.csv(metadata, "disturbance_total_metadata_top15.csv", row.names=TRUE)
write.csv(otu, "disturbance_total_otu_genus_top15.csv", row.names=TRUE)



#rarefaction, see how many reads we need to reach saturation for sampling unique ASVs
#https://ucdavis-bioinformatics-training.github.io/2021-May-Microbial-Community-Analysis/data_analysis/mca_part2
nice_colors = c("#999999", "#E69F00", "#56B4E9","#e98756","#c08160","#5800e6", "#CDDC49", "#C475D3", 
                "#E94B30", "#233F57", "#FEE659", "#A1CFDD", "#F4755E", "#D6F6F7","#EB6D58", "#6898BF")


#ps1 <- subset_samples(ps1, sample_sums(ps1) <=2000)
#play around with numbers, but 1000 seems good for the first dataset
otu1 <- as.data.frame(otu_table(ps1))
otu1 <- as.data.frame(t(otu1))
rare1 <- rarecurve(otu1,
                   step =50,
                   lwd=2,
                   col=nice_colors,
                   ylab="ASVs", 
                   main = "Rarefaction curve for all samples (Illumina Run 1)",
                   label=FALSE)

#second dataset has so many reads, i don't need to filter any of these samples
ps2 <- subset_samples(ps2, sample_sums(ps2) <=2000)
otu2 <- as.data.frame(otu_table(ps2))
otu2 <- as.data.frame(t(otu2))
rare2 <- rarecurve(otu2,
                   step =50,
                   lwd=2,
                   col=nice_colors,
                   ylab="ASVs", 
                   main = "Rarefaction curve for samples under 2000 reads (Illumina Run 2)",
                   label=FALSE)


