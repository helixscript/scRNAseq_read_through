library(dplyr)
library(Biostrings)
source('lib.R')

alignmentMinSeqId <- 97
refGenome <- '/home/ubuntu/AAVengeR/data/referenceGenomes/hg38.2bit'

r <- readr::read_tsv('filtered_HMM_hits.tsv')


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

R1_blat <- left_join(R1_blat, distinct(select(r, targetName, sample, cellBarcode)), by = c('qName' = 'targetName'))

R1_blat$tStart <- as.integer(R1_blat$tStart)
R1_blat$tEnd <- as.integer(R1_blat$tEnd)

invisible(file.remove(c('R1.ff', 'R1.psl')))


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

R2_blat$tStart <- as.integer(R2_blat$tStart)
R2_blat$tEnd   <- as.integer(R2_blat$tEnd)

invisible(file.remove(c('R2.trimmed.ff', 'R2.ff', 'R2.psl')))

saveRDS(list('R1' = R1_blat, 'R2' = R2_blat), 'HMM_hit_genome_alignments.rds')
