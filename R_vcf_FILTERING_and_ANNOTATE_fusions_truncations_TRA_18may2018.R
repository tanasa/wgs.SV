## devtools::install_github('seandavi/VCFWrenchR')
## install_github("d-cameron/StructuralVariantAnnotation")

# CRAN packages
library(devtools)
library(stringr)
library(dplyr)
library(ggplot2)
library(VCFWrenchR)

# bioconductor packages

library(VariantAnnotation)
library(VariantTools)

library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(org.Hs.eg.db)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(StructuralVariantAnnotation) 

#########################################################
#########################################################
#########################################################

args <- commandArgs(TRUE)

FILE_VCF <- args[1]  
       
name <- basename(FILE_VCF)

vcf <- readVcf(FILE_VCF, genome="hg38")

### samples(header(vcf))[1] ### 
### samples(header(vcf))[2] ### 

TUMOR  <-  samples(header(vcf))[1] 
NORMAL <-  samples(header(vcf))[2]

##########################################################
##########################################################
##########################################################

vcf.df <- as.data.frame(vcf)

vcf.df2 <- apply(vcf.df, 2, as.character)

write.table(vcf.df2, 
            file=paste(name, ".df", sep=""), 
            sep="\t", row.names=F)

#########################################################
#########################################################
#########################################################

# info(vcf)
# info(header(vcf))

# geno(vcf)
# geno(header(vcf))

# rowRanges(vcf)
# fixed(vcf)

########################

# fixed(vcf)
# ref(vcf)
# alt(vcf)
# qual(vcf)
# filt(vcf)

#########################

# colData(vcf)
# rowRanges(vcf)
# genome(vcf)
# seqlevels(vcf)

#########################################################
#########################################################
#########################################################

DEPTH_threshold <- 8
FRACTION_threshold <- 0.2 
 
##########################################################################
##########################################################################
##########################################################################

PASS_IMPRECISE_filters = function(x) { DV_germline <- ( geno(vcf)$DV[,NORMAL] < 1 )
                                       RV_germline <- ( geno(vcf)$RV[,NORMAL] < 1 )

                                       DV_tumor <- geno(vcf)$DV[,TUMOR]
                                       RV_tumor <- geno(vcf)$RV[,TUMOR]  

                                       DR_tumor <- geno(vcf)$DR[,TUMOR]
                                       RR_tumor <- geno(vcf)$RR[,TUMOR] 

                                       AD <-  (DV_tumor > DEPTH_threshold)      
                                       AF <- ((DV_tumor / (DV_tumor + DR_tumor)) > FRACTION_threshold)
                                       
                                       DV_germline &
                                       RV_germline & 
                                       AD &
                                       AF & 
                                       (filt(vcf) == "PASS") & 
                                       (info(vcf)$IMPRECISE) 
                                      }

vcf_PASS_IMPRECISE_FILTERS <- vcf[PASS_IMPRECISE_filters(vcf)]

writeVcf(vcf_PASS_IMPRECISE_FILTERS,
         file=paste(name, ".filter.PASS.IMPRECISE.vcf", sep=""))

##### transforming into a DATAFRAME : 
##### head(as.data.frame(vcf_PASS_IMPRECISE_FILTERS))

vcf_PASS_IMPRECISE_FILTERS.df <- as.data.frame(vcf_PASS_IMPRECISE_FILTERS)

##### here computing AD and AF :
##### AD = DV_TUMOR
##### AF = DV_TUMOR / (DV_TUMOR + DR_TUMOR)

vcf_PASS_IMPRECISE_FILTERS.df$AD <- vcf_PASS_IMPRECISE_FILTERS.df[,46]
vcf_PASS_IMPRECISE_FILTERS.df$AF <- vcf_PASS_IMPRECISE_FILTERS.df[,46] / (vcf_PASS_IMPRECISE_FILTERS.df[,46] + vcf_PASS_IMPRECISE_FILTERS.df[,44])

## vcf_PASS_IMPRECISE_FILTERS.df2 <- apply(vcf_PASS_IMPRECISE_FILTERS.df, 2, as.character)
vcf_PASS_IMPRECISE_FILTERS.df2 <- as.matrix(vcf_PASS_IMPRECISE_FILTERS.df)

write.table(vcf_PASS_IMPRECISE_FILTERS.df2, 
            file=paste(name, ".filter.PASS.IMPRECISE.df", sep=""), 
            sep="\t", row.names=F)

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

PASS_PRECISE_filters = function(x) {   DV_germline <- ( geno(vcf)$DV[,NORMAL] < 1 )
                                       RV_germline <- ( geno(vcf)$RV[,NORMAL] < 1 )

                                       DV_tumor <- geno(vcf)$DV[,TUMOR]
                                       RV_tumor <- geno(vcf)$RV[,TUMOR]  

                                       DR_tumor <- geno(vcf)$DR[,TUMOR]
                                       RR_tumor <- geno(vcf)$RR[,TUMOR] 

                                       AD <- ((DV_tumor + RV_tumor) > DEPTH_threshold)      
                                       AF <- ((RV_tumor / (RV_tumor + RR_tumor)) > FRACTION_threshold) 
                                       
                                       DV_germline &
                                       RV_germline & 
                                       AD &
                                       AF & 
                                       (filt(vcf) == "PASS") & 
                                       (info(vcf)$PRECISE) 
                                    }

vcf_PASS_PRECISE_FILTERS <- vcf[PASS_PRECISE_filters(vcf)]

writeVcf(vcf_PASS_PRECISE_FILTERS, 
         file=paste(name, ".filter.PASS.PRECISE.vcf", sep=""))

##### transforming into a DATAFRAME : 
##### head(as.data.frame(vcf_PASS_PRECISE_FILTERS))

vcf_PASS_PRECISE_FILTERS.df <- as.data.frame(vcf_PASS_PRECISE_FILTERS)

##### here computing AD and AF :
##### AD = DV_TUMOR + RV_TUMOR
##### AF = RV_TUMOR / (RV_TUMOR + RR_TUMOR)

vcf_PASS_PRECISE_FILTERS.df$AD <- vcf_PASS_PRECISE_FILTERS.df[,46] + vcf_PASS_PRECISE_FILTERS.df[,50] 
vcf_PASS_PRECISE_FILTERS.df$AF <- vcf_PASS_PRECISE_FILTERS.df[,50] / (vcf_PASS_PRECISE_FILTERS.df[,50] + vcf_PASS_PRECISE_FILTERS.df[,48])

## vcf_PASS_PRECISE_FILTERS.df2 <- apply(vcf_PASS_PRECISE_FILTERS.df, 2, as.character)
   vcf_PASS_PRECISE_FILTERS.df2 <- as.matrix(vcf_PASS_PRECISE_FILTERS.df)

write.table(vcf_PASS_PRECISE_FILTERS.df2, 
            file=paste(name, ".filter.PASS.PRECISE.df", sep=""), 
            sep="\t", row.names=F)

##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################

vcf_PASS_COMBINED_FILTERS <- rbind(vcf_PASS_IMPRECISE_FILTERS, vcf_PASS_PRECISE_FILTERS) 

writeVcf(vcf_PASS_COMBINED_FILTERS, 
         file=paste(name, ".filter.PASS.PRECISE.and.IMPRECISE.filter.AD.and.AF.vcf", sep="")) 

##### head(as.data.frame(vcf_PASS_PRECISE_FILTERS))
##### vcf_PASS_COMBINED_FILTERS.df <- as.data.frame(vcf_PASS_COMBINED_FILTERS)
##### or we can do :

vcf_PASS_COMBINED_FILTERS.df <- rbind(vcf_PASS_IMPRECISE_FILTERS.df, vcf_PASS_PRECISE_FILTERS.df) 

vcf_PASS_COMBINED_FILTERS.df$SV_ID <- rownames(vcf_PASS_COMBINED_FILTERS.df)

### vcf_PASS_COMBINED_FILTERS.df2 <- apply(vcf_PASS_COMBINED_FILTERS.df, 2, as.character)
vcf_PASS_COMBINED_FILTERS.df2 <- as.matrix(vcf_PASS_COMBINED_FILTERS.df)

write.table(vcf_PASS_COMBINED_FILTERS.df2, 
            file=paste(name, ".filter.PASS.PRECISE.and.IMPRECISE.filter.AD.and.AF.df", sep=""), 
            sep="\t", row.names=F)

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
##### here working with the dataframe : vcf_PASS_COMBINED_FILTERS 
# breakpointRanges
# extracts the structural variants as a BreakendGRanges
# convert to breakend GRanges

gr <- breakpointRanges(vcf_PASS_COMBINED_FILTERS)

# retain only primary chromosomes
####### seqlevelsStyle(gr) <- "UCSC"
####### gr <- gr[seqnames(gr) %in% paste0("chr", c(1:22, "X", "Y"))]
####### seqlevels(gr) <- paste0("chr", c(1:22, "X", "Y"))

# gr <- breakpointRanges(vcf)

###################################

# remove breakends that now don't have a partner (eg: chr1 -> chrMT)

gr <- gr[gr$partner %in% names(gr)]

## annotate breakends with gene names and gene orientation
## making a change : adding refGene 
## hg38 <-  makeTxDbFromUCSC(genome="hg38", tablename="refGene")
## gns <- genes(hg38)

gns <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

hits <- as.data.frame(findOverlaps(gr, gns, ignore.strand=TRUE))

#############################################################################################
#############################################################################################

file_exit = paste(name, ".if.breakpoints.intersecting.genes", sep="")

if (dim(hits)[1] == 0) { 
                          message <- "NO breakpoints intersecting the GENES"
                          write.table(message, file_exit, sep="\t", row.names=F)
                          stop("no breakpoints intersecting the genes", call=FALSE)
                       
} else {
                         message <- "there are breakpoints intersecting the GENES"
                         write.table(message, file_exit, sep="\t", row.names=F)
                         ## continue      
} 

############################################################################################
############################################################################################

## here it is a condition saying whether to continue the SCRIPT or to STOP
## depending on the value of GNS

############################################################################################
############################################################################################
############################################################################################

hits$SYMBOL <- biomaRt::select(org.Hs.eg.db, gns[hits$subjectHits]$gene_id, "SYMBOL")$SYMBOL

hits$gene_strand <- as.character(strand(gns[hits$subjectHits]))

hits <- hits %>%
  group_by(queryHits) %>%
  summarise(SYMBOL=paste(SYMBOL, collapse=","), gene_strand=paste0(gene_strand, collapse=""))

gr$SYMBOL <- ""
gr$geneStrand <- ""
gr$SYMBOL[hits$queryHits] <- hits$SYMBOL
gr$geneStrand[hits$queryHits] <- hits$gene_strand

#################################################################################################
################################################################################################# 
#################################################################################################

gr$couldBeThreePrimeStart <- str_detect(gr$geneStrand, stringr::fixed(as.character(strand(gr))))

gr$couldBeFivePrimeEnd <- str_detect(gr$geneStrand, stringr::fixed(ifelse(as.character(strand(gr))=="+", "-", "+")))

################################################################################################# 
#################################################################################################

### here printing the annotated dataframe with 3'START and 5'END :

gr.df <- as.data.frame(gr)

write.table(gr.df, 
            file=paste(name, ".filter.PASS.PRECISE.and.IMPRECISE.filter.AD.and.AF.with.35.info", sep=""), 
            sep="\t", row.names=F)

################### TRUNCATIONS ########################################################################
########################################################################################################
########################################################################################################
########################################################################################################
################## possibly now to use the dataframe GR or GR.DF in order to infer the TRUNCATIONS :

gr.use.truncations <- gr

gr.truncations <- gr.use.truncations[
                     (gr.use.truncations$couldBeThreePrimeStart & (! partner(gr.use.truncations)$couldBeFivePrimeEnd)) |
                     ( (! gr.use.truncations$couldBeFivePrimeEnd) & partner(gr.use.truncations)$couldBeThreePrimeStart), 
                  ]

gr.truncations.df <- as.data.frame(gr.truncations)

write.table(gr.truncations.df, 
            file=paste(name, ".filter.PASS.PRECISE.and.IMPRECISE.filter.AD.and.AF.with.TRUNCATIONS.Rform", sep=""), 
            sep="\t", row.names=F)

#######################################################################################################
########################################### here it is a dataframe where we write the TRUNCATIONS : 

TRUNCATIONS.df <- data.frame( ID = character(), 
                              Partner1 = character(), 
                              Partner2 = character(), 
                              stringsAsFactors=FALSE )


truncation.df <- data.frame( ID = character(), 
                             Partner1 = character(), 
                             Partner2 = character(), 
                             stringsAsFactors=FALSE )

###################################################################################################
########################################################### here below extracting the TRUNCATIONS :

ids <- unique(gr.truncations.df$vcfId)

for(id in ids)
{
   temp <- gr.truncations.df[gr.truncations.df$vcfId == id, ]

   if(nrow(temp) == 2) 
   {
         ### here it is the possibility of an orientation :  
 
         if (  (  temp$couldBeThreePrimeStart[1] == "TRUE") & 
               (  temp$couldBeFivePrimeEnd[2] == "FALSE") ) 
           {
               truncation.df <- data.frame( ID = unique(temp$vcfId[1], temp$vcfId[2]), 
                                            Partner1 = temp$SYMBOL[1] , 
                                            Partner2 = temp$SYMBOL[2] , 
                                            stringsAsFactors=FALSE 
                                          )
               
               TRUNCATIONS.df <- rbind(TRUNCATIONS.df, truncation.df)
           } 
 
         ### here it is the possibility of another orientation :
          
          if ( ( temp$couldBeThreePrimeStart[2] == "TRUE") & 
               ( temp$couldBeFivePrimeEnd[1] == "FALSE") ) 
           {
               truncation.df <- data.frame( ID = unique(temp$vcfId[1], temp$vcfId[2]), 
                                            Partner1 = temp$SYMBOL[2] , 
                                            Partner2 = temp$SYMBOL[1] , 
                                            stringsAsFactors=FALSE 
                                          )
               
               TRUNCATIONS.df <- rbind(TRUNCATIONS.df, truncation.df)
           } 
   }
}

############################################################################
### THE OUTPUT is TRUNCATIONS.df
###         head(TRUNCATIONS.df)
############################################################################
############################################################################
### now if we want to integrate the TRUNCATIONS.df with 
### the dataframe vcf_PASS_COMBINED_FILTERS.df:
############################################################################
############################################################################

TRUNCATIONS.df.with.INFO <- merge(TRUNCATIONS.df, 
                              as.data.frame(vcf_PASS_COMBINED_FILTERS.df), 
                              by.x = "ID", 
                              by.y = "SV_ID", 
                              all.x = TRUE)


write.table(as.matrix(TRUNCATIONS.df.with.INFO), 
            file=paste(name, ".filter.PASS.PRECISE.and.IMPRECISE.filter.AD.and.AF.with.TRUNCATIONS.txt", sep=""), 
            sep="\t", row.names=F)

#########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
####################################################################################################################
####################################################################################################################
################### now to use the INITIAL dataframe in order to infer the FUSIONS :

gr.use.fusions <- gr

# requires that THE BREAKPOINT to be BETWEEN DIFFERENT GENES

gr.use.fusions <- gr.use.fusions[gr.use.fusions$SYMBOL != partner(gr.use.fusions)$SYMBOL,]
gr.use.fusions <- gr.use.fusions[gr.use.fusions$SYMBOL != "" & partner(gr.use.fusions)$SYMBOL != "",]

####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

gr.fusions <- gr.use.fusions[(gr.use.fusions$couldBeThreePrimeStart & partner(gr.use.fusions)$couldBeFivePrimeEnd) |
                             (gr.use.fusions$couldBeFivePrimeEnd & partner(gr.use.fusions)$couldBeThreePrimeStart), 
                ]
gr.fusions.df <- as.data.frame(gr.fusions)

write.table(gr.fusions.df, 
            file=paste(name, ".filter.PASS.PRECISE.and.IMPRECISE.filter.AD.and.AF.with.FUSIONS.Rform", sep=""), 
            sep="\t", row.names=F)

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
### here it is a dataframe where we write the FUSIONS : 

FUSIONS.df <- data.frame( ID = character(), 
                          Partner1 = character(), 
                          Partner2 = character(), 
                          stringsAsFactors=FALSE )

fusion.df <- data.frame( ID = character(), 
                          Partner1 = character(), 
                          Partner2 = character(), 
                          stringsAsFactors=FALSE )

############################### here below extracting the FUSED GENES :

ids <- unique(gr.fusions.df$vcfId)

for(id in ids)
{

 temp <- gr.fusions.df[gr.fusions.df$vcfId == id, ]

 if(nrow(temp) == 2) 
 {

         ### here it is the possibility of an orientation :  
 
         if (  (  temp$couldBeThreePrimeStart[1] == "TRUE") & 
               ( temp$couldBeFivePrimeEnd[2] == "TRUE")  ) 
           {
               fusion.df <- data.frame( ID = unique(temp$vcfId[1], temp$vcfId[2]), 
                                        Partner1 = temp$SYMBOL[1] , 
                                        Partner2 = temp$SYMBOL[2] , 
                                        stringsAsFactors=FALSE 
                                       )
               
               FUSIONS.df <- rbind(FUSIONS.df, fusion.df)
           } 
 
          ### here it is the possibility of another orientation :
          
          if ( ( temp$couldBeThreePrimeStart[2] == "TRUE") & 
               ( temp$couldBeFivePrimeEnd[1] == "TRUE") ) 
           {
               fusion.df <- data.frame( ID = unique(temp$vcfId[1], temp$vcfId[2]), 
                                        Partner1 = temp$SYMBOL[2] , 
                                        Partner2 = temp$SYMBOL[1] , 
                                        stringsAsFactors=FALSE 
                                       )
               
               FUSIONS.df <- rbind(FUSIONS.df, fusion.df)
           } 
   }
}

############################################################################
### THE OUTPUT is FUSIONS.df
###         head(FUSIONS.df)
############################################################################
############################################################################
############################################################################
############################################################################

### now if we want to integrate the FUSIONS.df with 
### the dataframe vcf_PASS_COMBINED_FILTERS.df:

############################################################################
############################################################################

FUSIONS.df.with.INFO <- merge(FUSIONS.df, 
                              as.data.frame(vcf_PASS_COMBINED_FILTERS.df), 
                              by.x = "ID", 
                              by.y = "SV_ID", 
                              all.x = TRUE)


write.table(as.matrix(FUSIONS.df.with.INFO), 
            file=paste(name, ".filter.PASS.PRECISE.and.IMPRECISE.filter.AD.and.AF.with.FUSIONS.txt", sep=""), 
            sep="\t", row.names=F)

############################################################################
############################################################################
############################################################################
