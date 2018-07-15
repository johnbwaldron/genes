#Rterm.exe --vanilla
#install.packages(c("nnls", "Hmisc", "mgcv", "plyr"))
# Not necessary to load the libraries because they're loaded in the CANOES.R functions
#library("nnls")
#library("Hmisc")
#library("mgcv")
#library("plyr")
#install.packages("tidyverse")
#library("tidyverse")

#Read in the gc content per probe/capture oligos (do you call them probes when their in solution and tagged with beads?
gc <- read.table("https://raw.githubusercontent.com/johnbwaldron/genes/master/lungcancerwes81616.gc.txt")$V2
canoes.reads <- read.table("C:/Users/4wald/Documents/203samples_canoes.reads.txt")


#for whatever reason coming out of bedtools, our canoes.reads.txt column 4 is full of "." 
# so I need to delete those (which is done in unix in the existining pipeline). 
canoes.reads <- canoes.reads[,-4]

#rename columns 5 and onward with sample names:
# bam name only in CANOES run order found here: "https://raw.githubusercontent.com/johnbwaldron/genes/master/bamNameOnly.tsv"
sample.names <- read.table("https://raw.githubusercontent.com/johnbwaldron/genes/master/bamNameOnly.tsv", row.names="V1")
sample.names <- rownames(sample.names)
#sample.names <- paste("WES_GELCC", c(105:110), sep="") #just made up
names(canoes.reads) <- c("chromosome", "start", "end", sample.names)

#create a vector of consecutive target/probe ids
target <- seq(1, nrow(canoes.reads))

#cbind it all into one dataframe
canoes.reads <- cbind(target, gc, canoes.reads)

#to load the functions for CANOES
source("http://www.columbia.edu/~ys2411/canoes/CANOES.R") # also: source("https://raw.githubusercontent.com/johnbwaldron/genes/master/CANOES.R")

# create a vector to hold the results for each sample
xcnv.list <- vector('list', length(sample.names))

# call CNVs in each sample
for (i in 1:length(sample.names)){
    xcnv.list[[i]] <- CallCNVs(sample.names[i], canoes.reads)
  }

# combine the results into one data frame
xcnvs <- do.call('rbind', xcnv.list)

####################################################FOR FIGURING OUT GENOTYPING########################################################
#write.table(xcnvs, file="C:/Users/Public/xcnvs.csv", row.names=FALSE, quote=TRUE, sep=",")
#write.table(canoes.reads, file="C:/Users/Public/canoesReads.csv", row.names=FALSE, quote=TRUE, sep=",")
xcnvs <- read.table(file="C:/Users/Public/xcnvs.csv",  header = TRUE, sep=",")
canoes.reads <- read.table(file="C:/Users/Public/canoesReads.csv",  header = TRUE, sep=",")
counts <- subset(canoes.reads, !chromosome %in% c("chrX", "chrY", "X", "Y"))
source("https://raw.githubusercontent.com/johnbwaldron/genes/master/CANOES.R")
i<-2
genotypesBySample <- GenotypeCNVs(xcnvs, colnames(counts)[-c(1:5)][i], counts)
########################################################################################################################################

##########################Rather than reformatting canoes.reads each time, just load it below###########################
#xcnvs <- read.table(file="C:/Users/Public/xcnvs.csv",  header = TRUE, sep=",")
#canoes.reads <- read.table(file="C:/Users/Public/canoesReads.csv",  header = TRUE, sep=",")
# genotype all the CNV calls made above in samples
    # rather than append all CNVs for each sample to single file, I think I'll create a file for each genotyped cnv. 
    # with that goal in mind, there is no need to create the list, it makes sense to me to directly write the file with wach sample
    #dir.create("C:/Users/Public/GELCCgenotypes/")

for (i in 1:length(sample.names)){
    genotypesBySample <- GenotypeCNVs(xcnvs, sample.names[i], canoes.reads)

    write.table(genotypesBySample, file=paste("C:/Users/Public/GELCCgenotypes/", sample.names[i], "Genotypes.csv", sep=""), row.names=FALSE, quote=TRUE, sep=",")
    cat("completed ", i, " of ", length(sample.names),"\n")
    }

# GenotypeCNVs
#     Genotype CNVs in sample of interest
# Arguments:
#   xcnv
#     data frame with the following columns, and one row for each
#     CNV to genotype
#      INTERVAL: CNV coordinates in the form chr:start-stop
#      TARGETS: target numbers of CNV in the form start..stop
#               these should correspond to the target numbers in counts
#   sample.name:
#     sample to genotype CNVs in (should correspond to a column in counts)
#   counts: 
#     count matrix, first five columns should be 
#       target: consecutive numbers for targets (integer)
#       chromosome: chromosome number (integer-valued) 
#         (support for sex chromosomes to come)
#       start: start position of probe (integer)
#       end: end position of probe (integer)
#       gc: gc content (real between 0 and 1)
#       subsequent columns should include counts for each probe for samples
#   p:
#     average rate of occurrence of CNVs (real) default is 1e-08
#   D:
#     expected distance between targets in a CNV (integer) default is 70,000
#   Tnum:
#     expected number of targets in a CNV (integer) default is 6
#   numrefs
#     maximum number of reference samples to use (integer) default is 30
#     the weighted variance calculations will take a long time if too 
#     many reference samples are used
#   emission.probs and distances are for internal use only
# Returns: 
#   data frame with the following columns and one row for each genotyped CNV:
#      INTERVAL: CNV coordinates in the form chr:start-stop
#      NQDEL: a Phred-scaled quality score that sample.name has no deletion 
#             in the interval
#      SQDEL: a Phred-scaled quality score that sample.name has a deletion 
#             in the interval
#      NQDUP and SQDUP: same, but for a duplication
