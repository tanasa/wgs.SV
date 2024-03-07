library("reshape2")
library("ggplot2")
library("ggthemes")

library("plyr")
library("dplyr")
library("tidyr")
library("magrittr")

library("cowplot")
library("gridExtra")

library("superheat")
library("corrplot")

library("data.table")
library("splitstackshape")

######################################################################################
######################################################################################

# read.file.annovar <- function(file)
# {
#   NAMES <- read.table(file, nrow = 1, stringsAsFactors = FALSE, sep = "\t", quote="", fill=FALSE)
#   DATA <- read.table(file, skip = 1, stringsAsFactors = FALSE, sep = "\t", quote="", fill=FALSE)
#   DATA <- DATA[, 1:58]
#   names(DATA) <- NAMES 
#   return(DATA)
# }

#######################################################################################
#######################################################################################  
############################ reading the files whose name end in "FUSIONS"

files.names <- list.files(pattern="*with.FUSIONS.txt") 
# It reads the text files.

list.of.files <- lapply(files.names, read.delim, header=TRUE, stringsAsFactors=F)
# list.of.files <- lapply(files.names, read.table)
# list.of.files <- lapply(files.names, read.delim, header=TRUE, stringsAsFactors=F)
# it reads in all the files and places them in a list, with each element of the list being a data.frame.

names(list.of.files) <- files.names

length.list.of.files <- length(list.of.files)

#######################################################################################
### in order to decompose the name of the file into a SIMPLER NAMES :
### having the particle after vcf.... ie. vcf.NAME....
### the names of the SAMPLES will start with SPCG- ...
#######################################################################################

for (i in 1:length(list.of.files))
{
     short.name <- strsplit(names(list.of.files)[i],".", fixed=TRUE)[[1]][3]
     ## print(paste("The gene is", short.name))
     names(list.of.files)[[i]] <- short.name
}

########################################################################################

### here setting up the name of the SAMPLE :

NAME  <- unique(names(list.of.files))

########################################################################################
########################################################################################

### here combining ALL the FUSIONS :

ALL_fusions <- do.call(rbind.data.frame, list.of.files)

### rownames(ALL)
### colnames(ALL)

df.name.fusions <- data.frame(TYPE="FUSION")

if (dim(ALL_fusions)[1] != 0)
{
   ALL.and.sample.fusions <- cbind(df.name.fusions, NAME, ALL_fusions)

   TOTAL_FUSIONS = dim(ALL.and.sample.fusions)[1]  

   write.table(ALL.and.sample.fusions, file=paste(NAME, ".list.of.fusions.all.info", sep=""),
                                       sep="\t", 
                                       quote = FALSE, 
                                       row.names=FALSE)
}

if (dim(ALL_fusions)[1] == 0)
{
   ALL.and.sample.fusions <- NULL

   TOTAL_FUSIONS = 0

   write.table(ALL.and.sample.fusions, file=paste(NAME, ".list.of.fusions.all.info", sep=""),
                                       sep="\t", 
                                       quote = FALSE, 
                                       row.names=FALSE)
}

########################################################################################
########################################################################################
########################################################################################
############################ doing the same with TRUNCATIONS :
#######################################################################################
#######################################################################################  
############################ reading the files whose name end in "TRUNCATIONS"

files.names <- list.files(pattern="*with.TRUNCATIONS.txt") 
# It reads the text files.

list.of.files <- lapply(files.names, read.delim, header=TRUE, stringsAsFactors=F)
# list.of.files <- lapply(files.names, read.table)
# list.of.files <- lapply(files.names, read.delim, header=TRUE, stringsAsFactors=F)
# it reads in all the files and places them in a list, with each element of the list being a data.frame.

names(list.of.files) <- files.names

length.list.of.files <- length(list.of.files)

#######################################################################################
### in order to decompose the name of the file into a SIMPLER NAMES :
### having the particle after vcf.... ie. vcf.NAME....
### the names of the SAMPLES will start with SPCG- ...
#######################################################################################

for (i in 1:length(list.of.files))
{
     short.name <- strsplit(names(list.of.files)[i],".", fixed=TRUE)[[1]][3]
     ## print(paste("The gene is", short.name))
     names(list.of.files)[[i]] <- short.name
}

########################################################################################

### here setting up the name of the SAMPLE :

NAME  <- unique(names(list.of.files))

########################################################################################
########################################################################################

### here combining ALL the TRUNCATIONS :

ALL_truncations <- do.call(rbind.data.frame, list.of.files)

### rownames(ALL)
### colnames(ALL)
df.name.truncations <- data.frame(TYPE="TRUNCATION")

if (dim(ALL_truncations)[1] != 0)
{
   ALL.and.sample.truncations <- cbind(df.name.truncations, NAME, ALL_truncations)

   TOTAL_TRUNCATIONS = dim(ALL.and.sample.truncations)[1]

   write.table(ALL.and.sample.truncations, file=paste(NAME, ".list.of.truncations.all.info", sep=""),
                                                   sep="\t", 
                                                   quote = FALSE, 
                                                   row.names=FALSE)
}

if (dim(ALL_truncations)[1] == 0)
{
   ALL.and.sample.truncations <- NULL

   TOTAL_TRUNCATIONS = 0

   write.table(ALL.and.sample.truncations, file=paste(NAME, ".list.of.truncations.all.info", sep=""),
                                                   sep="\t", 
                                                   quote = FALSE, 
                                                   row.names=FALSE)
}

################################################################################################
################################################################################################
################################################################################################
############################ COMBINING THE DATAFRAMES with FUSIONS and TRUNCATIONS :
### ALL.and.sample.fusions
### ALL.and.sample.truncations
################################################################################################
################################################################################################
################################################################################################  

if ( (dim(ALL_fusions)[1] == 0) && (dim(ALL_truncations)[1] == 0) )
{

    ALL.and.sample.fusions.and.truncations <- NULL

    TOTAL_FUSIONS_TRUNCATIONS = 0

    write.table(ALL.and.sample.fusions.and.truncations, 
                           file=paste(NAME, ".list.of.fusions.and.truncations.all.info", sep=""),
                                      sep="\t", 
                                      quote = FALSE, 
                                      row.names=FALSE)
}

if ( (dim(ALL_fusions)[1] != 0) || (dim(ALL_truncations)[1] != 0) )
{

    ALL.and.sample.fusions.and.truncations <- rbind(ALL.and.sample.fusions, ALL.and.sample.truncations) 

    TOTAL_FUSIONS_TRUNCATIONS = dim(ALL.and.sample.fusions.and.truncations)[1]

    write.table(ALL.and.sample.fusions.and.truncations, 
                           file=paste(NAME, ".list.of.fusions.and.truncations.all.info", sep=""),
                                      sep="\t", 
                                      quote = FALSE, 
                                      row.names=FALSE)
}

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
##################### here combining the files that contain FILTERED SV, function of AD and AF :

files.names <- list.files(pattern="*filter.AD.and.AF.df") 
# It reads the text files.

list.of.files <- lapply(files.names, read.delim, header=TRUE, stringsAsFactors=F)
# list.of.files <- lapply(files.names, read.table)
# list.of.files <- lapply(files.names, read.delim, header=TRUE, stringsAsFactors=F)
# it reads in all the files and places them in a list, with each element of the list being a data.frame.

names(list.of.files) <- files.names

length.list.of.files <- length(list.of.files)

#######################################################################################
### in order to decompose the name of the file into a SIMPLER NAMES :
### having the particle after vcf.... ie. vcf.NAME....
### the names of the SAMPLES will start with SPCG- ...
#######################################################################################

for (i in 1:length(list.of.files))
{
     short.name <- strsplit(names(list.of.files)[i],".", fixed=TRUE)[[1]][3]
     ## print(paste("The gene is", short.name))
     names(list.of.files)[[i]] <- short.name
}

########################################################################################

### here setting up the name of the SAMPLE :

NAME  <- unique(names(list.of.files))

########################################################################################
########################################################################################

### here combining ALL the FILTERED SV :

ALL_SV <- do.call(rbind.data.frame, list.of.files)

### rownames(ALL)
### colnames(ALL)
df.name.SV <- data.frame(TYPE="SV")

if (dim(ALL_SV)[1] != 0)
{
   ALL.and.sample.SV <- cbind(df.name.SV, NAME, ALL_SV)

   write.table(ALL.and.sample.SV, file=paste(NAME, ".list.of.filtered.SV.all.info", sep=""),
                                             sep="\t", 
                                             quote = FALSE, 
                                             row.names=FALSE)
}

if (dim(ALL_SV)[1] == 0)
{
   ALL.and.sample.SV <- NULL

   write.table(ALL.and.sample.SV, file=paste(NAME, ".list.of.filtered.SV.all.info", sep=""),
                                             sep="\t", 
                                             quote = FALSE, 
                                             row.names=FALSE)
}

########################################################################################
########################################################################################
########################################################################################
########################################################################################
###
### here setting up a new DATAFRAME that counts the number of DEL, DUP, INS, INV, TRA

DEL <- sum(grepl("DEL", ALL.and.sample.SV$SV_ID))
DUP <- sum(grepl("DUP", ALL.and.sample.SV$SV_ID))
INS <- sum(grepl("INS", ALL.and.sample.SV$SV_ID))
INV <- sum(grepl("INV", ALL.and.sample.SV$SV_ID))
TRA <- sum(grepl("TRA", ALL.and.sample.SV$SV_ID))

TOTAL_SV_FILTERED = DEL + DUP + INS + INV + TRA

# TOTAL_FUSIONS = dim(ALL.and.sample.fusions)[1]
# TOTAL_TRUNCATIONS = dim(ALL.and.sample.truncations)[1]
# TOTAL_FUSIONS_TRUNCATIONS = dim(ALL.and.sample.fusions.and.truncations)[1]

summary.df <- data.frame("NAME" = NAME,
                         "DEL" = DEL, 
                         "DUP" = DUP, 
                         "INS" = INS,
                         "INV" = INV, 
                         "TRA" = TRA,
                         "TOTAL_SV_FILTERED" = TOTAL_SV_FILTERED,
                         "TOTAL_FUSIONS" = TOTAL_FUSIONS, 
                         "TOTAL_TRUNCATIONS" = TOTAL_TRUNCATIONS, 
                         "TOTAL_FUSIONS_TRUNCATIONS" = TOTAL_FUSIONS_TRUNCATIONS,     
                         stringsAsFactors=FALSE)

write.table(summary.df, file=paste(NAME, ".list.of.filtered.SV.all.info.SUMMARY", sep=""),
                                          sep="\t", 
                                          quote = FALSE, 
                                          row.names=FALSE)  
 
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################                                
