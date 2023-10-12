library(dplyr)
library(GenomicRanges)

o <- readRDS('HMM_hit_genome_alignments.rds')
R1_blat <- o$R1
R2_blat <- o$R2
rm(o)

# Group blat results by read id and identify the hit with the highest match score.
# Include other hits that are within {maxMatchDistFromLeader} of the highest score.
# This will often return a single result which can be considered the proper alignments. 
# For each hit, look for correponding R2 hit within 10000 NT.
#--------------------------------------------------------------------------------------------------

maxMatchDistFromLeader <- 3

r <- distinct(bind_rows(lapply(split(R1_blat, R1_blat$qName), function(x){
  m <- dplyr::filter(x, matches >= max(x$matches) - maxMatchDistFromLeader)
  
  m <- bind_rows(lapply(split(m, 1:nrow(m)), function(x2){
    if(x2$strand == '+'){
      pos1 <- x2$tStart
      pos2 <- pos1 + 10000
      k <- subset(R2_blat, sample == x2$sample & tName == x2$tName & strand != x2$strand & tEnd >= pos1 & tEnd <= pos2)
    } else {
      pos2 <- x2$tEnd
      pos1 <- pos2 - 10000
      k <- subset(R2_blat, sample == x2$sample & tName == x2$tName & strand != x2$strand & tStart >= pos1 & tStart <= pos2)
    }
    
    x2$r_supported <- FALSE
    if(nrow(k) > 0) x2$r_supported <- TRUE
    
    x2
  }))
  
  # Eliminate non-supported hits if only one hit is supported.
  if(nrow(m) > 1 & sum(m$r_supported == 1)) m <- m[m$r_supported,] 
  
  posids <- paste0(gtools::mixedsort(paste0(m$tName, m$strand, ifelse(m$strand == '+', m$tStart, m$tEnd))), collapse = ';')
  
  tibble(sample = m$sample, targetName = m$qName[1], 
         cellBarCode = paste0(unique(m$cellBarcode), collapse = ';'), 
         hits = posids, r_supported = any(m$r_supported))
})))

o <- bind_rows(lapply(split(r, paste(r$sample, r$hits)), function(x){
  tibble(sample = x$sample[1], 
         reads = n_distinct(x$targetName),
         hits = x$hits[1], 
         R2_supported = x$r_supported[1], 
         cellBarCodes = n_distinct(x$cellBarCode))
}))

# Remove known artifacts.
o <- o[! o$hits %in% c('chr3+17402188','chr6-131054616'),]

openxlsx::write.xlsx(r, 'readHits.xlsx')
openxlsx::write.xlsx(o, 'readHitsSummary.xlsx')



# # Correct instances where positions ids vary by a couple of nucleotides.
# #--------------------------------------------------------------------------------------------------
# single <- r[! grepl(';', r$hits),]
# multiple <- r[grepl(';', r$hits),]
# 
# d <- makeGRangesFromDataFrame(
#   data.frame(
#     seqnames = unlist(lapply(strsplit(single$hits, '[\\+\\-]'), '[', 1)),
#     start    = as.integer(unlist(lapply(strsplit(single$hits, '[\\+\\-]'), '[', 2))) - 5,
#     end      = as.integer(unlist(lapply(strsplit(single$hits, '[\\+\\-]'), '[', 2))) + 5,
#     strand   = stringr::str_extract(single$hits, '[\\+\\-]'),
#     sample   = single$sample, 
#     posid    = single$hits), 
#   keep.extra.columns = TRUE)
# 
# mce <- function(x) names(sort(table(x), decreasing = TRUE))[1]
# 
# d2 <- Reduce('append', lapply(split(d, d$sample), function(x){
#   k <- reduce(x, with.revmap = TRUE)
#   
#   Reduce('append', lapply(k$revmap, function(x2){
#     z <- x[x2]
#     z$posid2 <- mce(z$posid)
#     z
#   }))
# }))
# 
# single <- left_join(single, distinct(tibble(hits = d2$posid, hits2 = d2$posid2)), by = 'hits')
# single$hits <- single$hits2
# single$hits2 <- NULL
# 
# r <- bind_rows(single, multiple)
