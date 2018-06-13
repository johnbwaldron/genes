#Rterm.exe --vanilla
#install.packages(c("nnls", "Hmisc", "mgcv", "plyr"))
# Not necessary to load the libraries because their loaded in the CANOES.R functions
#library("nnls")
#library("Hmisc")
#library("mgcv")
#library("plyr")

#Read in the gc content per probe/capture oligos (do you call them probes when their in solution and tagged with beads?
gc <- read.table("https://raw.githubusercontent.com/johnbwaldron/genes/master/lungcancerwes81616.gc.txt")$V2
canoes.reads <- read.table("https://raw.githubusercontent.com/johnbwaldron/genes/master/WES_GELCC105.canoes.reads.txt")

#for whatever reason coming out of bedtools, our canoes.reads.txt column 4 is full of "." 
# so I need to delete those (which is done in unix in the existining pipeline). 
canoes.reads <- canoes.reads[,-4]

#rename columns 5 and onward with sample names:
sample.names <- paste("WES_GELCC", c(105:110), sep="") #just made up
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

write.table(xcnvs, file="C:/Users/Public/xcnvs.csv", row.names=FALSE, quote=TRUE, sep=",")
