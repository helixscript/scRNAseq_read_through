library(dplyr)
library(ShortRead)
library(ggplot2)
source('lib.R')

alignmentMinSeqId <- 97
refGenome <- '/home/ubuntu/AAVengeR/data/referenceGenomes/blat/hg38.2bit'
TUs   <- '/home/ubuntu/AAVengeR/data/genomeAnnotations/hg38.TUs.rds'
exons <- '/home/ubuntu/AAVengeR/data/genomeAnnotations/hg38.exons.rds'

# Read in 10x allowed cell bar code list.
w <- readLines('cellRangerBarcodeWhiteList.txt')

# Collate the results from run_FASTQ_hmm.R
r <- bind_rows(lapply(list.files('output', pattern = '^result$', recursive = TRUE, full.names = TRUE), function(x){
       o <- readr::read_tsv(x)
       o$sample <- unlist(strsplit(x, '/'))[2]
       select(o, sample, targetName, tlen, fullEval, fullScore, targetStart, targetEnd, envStart, envEnd, hmmStart, hmmEnd, R1_seq, R2_seq) %>% distinct()
      }))

readr::write_tsv(r, 'allHMM_hits.tsv')


# Select the strongest HMM hit for each read.
# Break ties by selecting hits with HMM alignments closer to the ends of the LTR HMM.
r <- group_by(r, targetName) %>%
     dplyr::slice_max(fullScore, with_ties = TRUE) %>% 
     dplyr::slice_max(hmmEnd, with_ties = FALSE) %>%
     ungroup()


# Extract cell and transcript bar codes.
r$cellBarcode <- substr(r$R1_seq, 1, 16)
r$transcriptBarcode <- substr(r$R1_seq, 17, 26)


# Correct cell bar codes that are 1 off from whitelist.
a <- r[r$cellBarcode %in% w,]
b <- r[! r$cellBarcode %in% w,]

for(x in unique(b$cellBarcode)){
  m <- stringdist::stringdist(x, w)
  i <- which(m == 1)
  
  if(length(i) == 1){
    b[which(b$cellBarcode == x),]$cellBarcode <- w[i]
  } else {
    b[which(b$cellBarcode == x),]$cellBarcode <- NA
  }
}

b <- b[! is.na(b$cellBarcode),]
r <- bind_rows(a, b)


# Apply a mild filter to candidate reads.
r <- subset(r, fullScore >= 20 & (r$targetEnd - r$targetStart + 1) >= 20)
r$TSO_match <- vcountPattern('TTTCTTATATGGG', DNAStringSet(r$R1_seq), max.mismatch = 3) == 1



# Run R1 blat
#--------------------------------------------------------------------------------------------------
o <- DNAStringSet(substr(r$R1_seq, r$targetEnd+1, nchar(r$R1_seq)))
names(o) <-  r$targetName
o <- o[width(o) >= 15]
writeXStringSet(o, 'R1.ff')

system(paste('blat -stepSize=5 -repMatch=5000 -minScore=0 -minIdentity=0 -out=psl -noHead', 
             refGenome, 'R1.ff', 'R1.psl'))

R1_blat <- parseBLAToutput('R1.psl') %>% 
  dplyr::filter(alignmentPercentID >= alignmentMinSeqId & 
                  tNumInsert <= 1 & qNumInsert <= 1 & tBaseInsert <= 1 & qBaseInsert <= 1 & qStart <= 5)

R1_blat <- left_join(R1_blat, distinct(select(r, targetName, sample, cellBarcode, targetStart, targetEnd, hmmStart, hmmEnd)), by = c('qName' = 'targetName'))

invisible(file.remove(c('R1.ff', 'R1.psl')))

# Select matches within 3 of the best match.
R1_blat <- bind_rows(lapply(split(R1_blat, R1_blat$qName), function(x){
              o <- subset(x, matches >= (max(x$matches) - 3))
              k <- subset(r, targetName == x$qName[1]) %>% dplyr::slice_max(fullScore, with_ties = FALSE)
              o$hmmStart <- k$hmmStart
              o$hmmEnd   <- k$hmmEnd
              o
           }))


# Run R2 blat
#--------------------------------------------------------------------------------------------------
o <- DNAStringSet(r$R2_seq)
names(o) <-  r$targetName
writeXStringSet(o, 'R2.ff')

# Trim off NTs if they appear to be from short fragments where we read back into the LTR. 
# Here we use cutadapt and the reverse compliment of the end of the LTR for trimming.

system('cutadapt -e 0.15 -a TGCTAGAGATTT R2.ff > R2.trimmed.ff')
o <- readDNAStringSet('R2.trimmed.ff')
o <- o[width(o) >= 15]
writeXStringSet(o, 'R2.ff')

system(paste('blat -stepSize=5 -repMatch=5000 -minScore=0 -minIdentity=0 -out=psl -noHead', 
             refGenome, 'R2.ff', 'R2.psl'))

R2_blat <- parseBLAToutput('R2.psl') %>% 
  dplyr::filter(alignmentPercentID >= alignmentMinSeqId & 
                  tNumInsert <= 1 & qNumInsert <= 1 & tBaseInsert <= 1 & qBaseInsert <= 1 & qStart <= 5)

R2_blat <- left_join(R2_blat, distinct(select(r, targetName, sample, cellBarcode)), by = c('qName' = 'targetName'))

invisible(file.remove(c('R2.trimmed.ff', 'R2.ff', 'R2.psl')))

# Select matches within 3 of the best match.
R2_blat <- bind_rows(lapply(split(R2_blat, R2_blat$qName), function(x){
             subset(x, matches >= (max(x$matches) - 3))
           }))

R1_blat <- distinct(select(R1_blat, matches, strand, qName, qStart, qEnd, tName, tStart, tEnd, sample, cellBarcode, targetStart, targetEnd, hmmStart, hmmEnd))
R2_blat <- distinct(select(R2_blat, matches, strand, qName, qStart, qEnd, tName, tStart, tEnd, sample, cellBarcode))

# Use the R2 alignments to resolve instances when R1 aligns to multiple locations.
# Use GenomicRanges::reduce() with a large margin to determine if reads align near one another.

reconcile_alignments <- function(R1, R2){
  R2$strand <- ifelse(R2$strand == '+', '-', '+')  # Flip strands to allow merging of sites with reduce().
    
  o <- bind_rows(lapply(1:nrow(R1), function(x){
           a <- GenomicRanges::makeGRangesFromDataFrame(R1[x,], seqnames.field = 'tName', start.field = 'tStart', end.field = 'tEnd')
           b <- GenomicRanges::makeGRangesFromDataFrame(R2, seqnames.field = 'tName', start.field = 'tStart', end.field = 'tEnd')
           o <- c(a, b)
           d <- GenomicRanges::reduce(o, min.gapwidth = 1000, with.revmap = TRUE)  # Attempt to merge ranges even if they are up to 1000 NT apart.

           z <- NULL
           if(length(d) < length(o)){
             i <- unlist(lapply(d$revmap, function(x) 1 %in% x))  # Which revmap has sequence 1 which cones from this R1 alignment (R1[x,])
             k <- unlist(d[which(i)]$revmap)                      # Pull the revmap with sequence 1.
             z <- R2[k[k != 1] - 1,]                              # Pull the associated R2 sequence
             z$strand <- ifelse(z$strand == '+', '-', '+')        # Restore the original strand.
           }
      
           for(n in c('R2_matches', 'R2_strand', 'R2_tName', 'R2_tStart', 'R2_tEnd')){
                 R1[x, n] = ifelse(! is.null(z), z[1, sub('R2_', '', n)] , NA)
           }
           
           R1[x,]
         }))
    
  if(nrow(o) > 1 & ! all(is.na(o$R2_matches))) o <- o[! is.na(o$R2_matches),]
  
  o
}


z <- bind_rows(lapply(unique(c(R1_blat$qName, R2_blat$qName)), function(x){
       message(x)
       R1 <- subset(R1_blat, qName == x)
       R2 <- subset(R2_blat, qName == x)
  
       if(nrow(R1) == 0){
         # There is no R1 for this read ID -- return nothing.
         return(tibble())
       } else if (nrow(R1) > 0 & nrow(R2) == 0){
         # There is no R2 for this read ID, only return R1.
         return(R1)
       } else {
         reconcile_alignments(R1, R2)
       }
    }))


# Select uniquely mapped reads.
u <- z[! duplicated(z$qName),]


# Remove known artifacts that arise from reading out of U5 into vector or vector-plasmid.
# These positions were identified from earlier versions of this analysis from positions that 
# spanned multiple subjects.
a1 <- subset(u, tName == 'chr3' & tStart >= 17402188-1000 & tStart <= 17402188+1000)$qName
a2 <- subset(u, tName == 'chr6' & tStart >= 131054616-1000 & tStart <= 131054616+1000)$qName
u <- subset(u, ! qName %in% c(a1, a2))

# Add nearest gene annotations.
u$posid <- paste0(u$tName, u$strand, ifelse(u$strand == '+', u$tStart, u$tEnd))
nearestGenes <- nearestGene(u$posid, readRDS(TUs), readRDS(exons))
nearestGenes$posid <- paste0(nearestGenes$chromosome, nearestGenes$strand, nearestGenes$position)
u <- left_join(u, distinct(select(nearestGenes, posid, nearestGene)), by = 'posid')


# Identify genes with candidate hits nearby that have HMM hits near the ends of the HMM.
s <- group_by(subset(u, hmmEnd >= 95), nearestGene, sample, posid) %>%
     summarise(nCellBarcodes = n_distinct(cellBarcode), cellBarCodes = paste0(unique(cellBarcode), collapse = ';')) %>%
     ungroup() %>%
     arrange(desc(nCellBarcodes)) %>% filter(nearestGene != 'OSBP')
openxlsx::write.xlsx(s, 's.xlsx')
