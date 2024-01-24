library(data.table)
library(stringr)

data <- read.csv("run1_samples.txt", header = FALSE)

fn1 <- as.data.frame(data[data$V1 %like% "R1",])
colnames(fn1)[1] <- "filename1"
fn2 <- as.data.frame(data[data$V1 %like% "R2",])
colnames(fn2)[1] <- "filename2"

fn1[c("Sample", "unused")] <- str_split_fixed(fn1$filename1, "_", 2)
fn2[c("Sample", "unused")] <- str_split_fixed(fn2$filename2, "_", 2)

total_run1 <- merge(fn1, fn2, by = "Sample")
total_run1 <- total_run1[,c(1,2,4)]


data2 <- read.csv("run2_samples.txt", header = FALSE)

fn12 <- as.data.frame(data2[data2$V1 %like% "R1",])
colnames(fn12)[1] <- "filename1"
fn22 <- as.data.frame(data2[data2$V1 %like% "R2",])
colnames(fn22)[1] <- "filename2"

fn12[c("Sample", "unused")] <- str_split_fixed(fn12$filename1, "_", 2)
fn22[c("Sample", "unused")] <- str_split_fixed(fn22$filename2, "_", 2)

total_run2 <- merge(fn12, fn22, by = "Sample")
total_run2 <- total_run2[,c(1,2,4)]

total_run <- rbind(total_run1, total_run2)

write.csv(total_run, "ncbi_metadata.csv", row.names = FALSE)
