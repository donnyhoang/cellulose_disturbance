library(dada2); packageVersion("dada2")
#reference this https://benjjneb.github.io/dada2/tutorial_1_8.html


#change path as needed for the two different Illumina runs. We'll combine the datasets later since they may have different error rates
path <- "disturbance_run2/"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE)) #import forward reads
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE)) #import reverse reads
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1) #extract sample names


plotQualityProfile(fnFs[1:2]) #plot quality score of forward reads
plotQualityProfile(fnRs[1:2]) #plot quality score of reverse reads

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz")) #create object holding filtered F reads.
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz")) #create object holding filtered R reads
names(filtFs) <- sample.names #name
names(filtRs) <- sample.names #name

#This step takes awhile
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, #specify input+output
                     truncLen=c(240,160), #where we are cutting off our sequences
                     trimLeft = c(17,21), #play with these, to see if they improve sequence quality (# reads merged, kept, etc)
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, #default parameters
                     compress=TRUE, multithread=FALSE) #my multitread is set to FALSE bc I run on Windows.

head(out)

errF <- learnErrors(filtFs, multithread=FALSE) #my multitread is set to FALSE bc I run on Windows.
errR <- learnErrors(filtRs, multithread=FALSE) #my multitread is set to FALSE bc I run on Windows.
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=FALSE) #merge paired reads
head(mergers[[1]]) # inspect the merger data.frame from the first sample

seqtab <- makeSequenceTable(mergers) #create sequence table that is imput for chimera removing function
dim(seqtab) #Check number of ASV's before chimera removal

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=FALSE) #my multitread is set to FALSE bc I run on Windows.

dim(seqtab.nochim)#check number of ASV's after chimera removal

sum(seqtab.nochim)/sum(seqtab) 


getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#This step takes a few minutes.
taxa <- assignTaxonomy(seqtab.nochim, "tax/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE) #my multitread is set to FALSE bc I run on Windows.
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

seqs <- colnames(seqtab.nochim)
headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  headers[i] <- paste(">ASV", i, sep="_")
}

#fasta file
fasta <- c(rbind(headers, seqs)) #make fasta file
write(fasta, "disturbance_run2.fa") #write out as .fasta file

#count table
count_tab <- t(seqtab.nochim) #make count table
row.names(count_tab) <- sub(">", "", headers) #replace rownames with something that is easy to follow
write.csv(count_tab, "disturbance_run2_counts.csv", row.names=TRUE) #write out as .csv file


# tax table:
tax <- taxa #make taxonomy table
row.names(tax) <- sub(">", "", headers) #replace rownames with something that is easy to follow
write.csv(tax, "disturbance_run2_taxonomy.csv", row.names=TRUE) #write out as .csv file
write.csv(out, "reads_filtered_run2.csv", row.names=TRUE) #write out as .csv file
write.csv(track, "reads_tracked_run2.csv", row.names=TRUE) #write out as .csv file
