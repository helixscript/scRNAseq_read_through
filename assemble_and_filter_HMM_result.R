library(dplyr)
library(ShortRead)
source('lib.R')

w <- readLines('cellRangerBarcodeWhiteList.txt')

r <- bind_rows(lapply(list.files('output', pattern = '^result$', recursive = TRUE, full.names = TRUE), function(x){
       o <- readr::read_tsv(x)
       o$sample <- unlist(strsplit(x, '/'))[2]
       select(o, sample, targetName, tlen, fullEval, fullScore, targetStart, targetEnd, envStart, envEnd, hmmStart, hmmEnd, R1_seq, R2_seq) %>% distinct()
      }))

readr::write_tsv(r, 'allHMM_hits.tsv')


# Apply a mild filter to candidate reads.
r$TSO_match <- vcountPattern('TTTCTTATATGGG', subseq(DNAStringSet(r$R1_seq), 25, 40), max.mismatch = 3) == 1
r2 <- subset(r, hmmEnd >= 98 & TSO_match == TRUE)


# Extract cell and transcript bar codes.
r2$cellBarcode <- substr(r2$R1_seq, 1, 16)
r2$transcriptBarcode <- substr(r2$R1_seq, 17, 26)


# Correct cell bar codes that are 1 off from whitelist.
a <- r2[r2$cellBarcode %in% w,]
b <- r2[! r2$cellBarcode %in% w,]

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
r2 <- bind_rows(a, b)

readr::write_tsv(r2, 'filtered_HMM_hits.tsv')
