
parseBLAToutput <- function(f){
  if(! file.exists(f) | file.info(f)$size == 0) return(tibble::tibble())
  b <- readr::read_delim(f, delim = '\t', col_names = FALSE, col_types = readr::cols())
  
  names(b) <- c('matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand',
                'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts')
  
  b$tEnd <- b$tEnd - 1
  
  b$queryPercentID       <- (b$matches/b$qSize)*100
  b$tAlignmentWidth      <- (b$tEnd - b$tStart) + 1
  b$queryWidth           <- (b$qEnd - b$qStart) + 1
  b$alignmentPercentID   <- (b$matches/b$tAlignmentWidth)*100
  b$percentQueryCoverage <- (b$queryWidth/b$qSize)*100
  b$qStarts              <- as.character(b$qStarts)
  b$tStarts              <- as.character(b$tStarts)
  b
}



nearestGene <- function(posids, genes, exons, CPUs = 20){
  library(dplyr)
  library(parallel)
  
  o <- base::strsplit(posids, '[\\+\\-]')
  d <- tibble(chromosome = unlist(lapply(o, '[', 1)),
              position = unlist(lapply(o, '[', 2)),
              strand = stringr::str_extract(posids, '[\\+\\-]'))
  d$n <- ntile(1:nrow(d), CPUs)
  d$position <- as.integer(d$position)
  
  cluster <- makeCluster(CPUs)
  clusterExport(cluster, c('genes', 'exons'), envir = environment())
  
  r <- bind_rows(parLapply(cluster, split(d, d$n), function(x){
    #r <- bind_rows(lapply(split(d, d$n), function(x){
    library(dplyr)
    library(GenomicRanges)
    
    x$exon <- NA
    x$nearestGene <- NA
    x$nearestGeneStrand <- NA
    x$nearestGeneDist <- NA
    x$inGene <- FALSE
    x$inExon <- FALSE
    x$beforeNearestGene <- NA
    
    r <- GenomicRanges::makeGRangesFromDataFrame(x, 
                                                 seqnames.field = 'chromosome',
                                                 start.field = 'position',
                                                 end.field = 'position')
    
    o <- suppressWarnings(GenomicRanges::findOverlaps(r, exons, select='all', ignore.strand=TRUE, type='any'))
    
    if(length(queryHits(o)) > 0){
      for(i in unique(queryHits(o))){
        x[i,]$exon <- paste0(unique(exons[subjectHits(o[queryHits(o) == i]),]$name2), collapse = ', ')
        x[i,]$nearestGeneStrand <- paste0(unique(as.character(strand(exons[subjectHits(o[queryHits(o) == i])]))), collapse = ',')
      }
      
      if(any(! is.na(x$exon))){
        x[! is.na(x$exon),]$nearestGene <- x[! is.na(x$exon),]$exon
        x[! is.na(x$exon),]$nearestGeneDist <- 0
      }
    }
    
    x1 <- tibble()
    x2 <- tibble()
    
    if(any(! is.na(x$exon))) x1 <- x[! is.na(x$exon),]  # Exon hits
    if(any(is.na(x$exon)))   x2 <- x[is.na(x$exon),]    # Non-exon hits
    
    if(nrow(x2) > 0){
      r <- GenomicRanges::makeGRangesFromDataFrame(x2, 
                                                   seqnames.field = 'chromosome',
                                                   start.field = 'position',
                                                   end.field = 'position')
      
      o <- suppressWarnings(GenomicRanges::distanceToNearest(r, genes, select='all', ignore.strand=TRUE))
      
      if(length(queryHits(o)) > 0){
        for(i in unique(queryHits(o))){
          x2[i,]$nearestGene <- paste0(unique(genes[subjectHits(o[queryHits(o) == i]),]$name2), collapse = ', ')
          x2[i,]$nearestGeneStrand <- paste0(unique(as.character(strand(genes[subjectHits(o[queryHits(o) == i])]))), collapse = ',')
          x2[i,]$nearestGeneDist <- unique(mcols(o[queryHits(o) == i])$distance) # Always single?
          x2[i,]$beforeNearestGene <- ifelse(x2[i,]$position < min(start(genes[subjectHits(o[queryHits(o) == i])])), TRUE, FALSE)
        }
      }
      
      x2 <- x2[! is.na(x2$nearestGene),]
      
      if(any(x2$nearestGeneDist == 0)) x2[which(x2$nearestGeneDist == 0),]$beforeNearestGene <- NA
    }
    
    if(nrow(x1) > 0  & nrow(x2) > 0)  x <- bind_rows(x1, x2)
    if(nrow(x1) == 0 & nrow(x2) > 0)  x <- x2
    if(nrow(x1) > 0  & nrow(x2) == 0) x <- x1
    
    x
  }))
  
  stopCluster(cluster)
  
  if(any(r$nearestGeneDist == 0)) r[which(r$nearestGeneDist == 0),]$inGene <- TRUE
  if(any(! is.na(r$exon))) r[which(! is.na(r$exon)),]$inExon <- TRUE
  if('n' %in% names(r)) r$n <- NULL
  if('exon' %in% names(r)) r$exon <- NULL
  r
}



alignChunk <- function(f){
  library(dplyr)
  library(ShortRead)
  source(file.path(analysisDir, 'lib.R'))
  message(f)
  
  t2 <- tempfile(tmpdir = file.path(analysisDir, 'tmp'))
  
  # Read in sequence data chunk, quality trim to Q20, export as fasta for blat.
  # gDNA are expected to start at position 40 after cell bar code + transcript code + TSO.
  # reads <- trimTails(readFastq(f), 2, '5', 5)
  
  reads <- readFastq(f)
  
  if(length(reads) == 0) return(tibble())
  
  reads@id <- BStringSet(sub('\\s.+$', '', reads@id))
  
  barCodes <- tibble(readID = as.character(reads@id), 
                     cellBarCode = as.character(narrow(reads, 1, 16)@sread),
                     transcriptBarCode = as.character(narrow(reads, 17, 26)@sread))  
  
  reads <- narrow(reads, 41, width(reads))
  
  o <- reads@sread
  names(o) <- reads@id
  
  # Process all reads so that reads in earlier analysis are not lost.
  rm(reads)
  gc()
  
  writeXStringSet(o, paste0(t2, '.fasta'))
  
  system(paste('blat -stepSize=5 -repMatch=5000 -minScore=0 -minIdentity=0 -out=psl -noHead', 
               paste0(t, '.2bit'), paste0(t2, '.fasta'), paste0(t2, '.psl')))
  
  b <- parseBLAToutput(paste0(t2, '.psl'))
  
  if(nrow(b) > 0){
    b <- dplyr::filter(b, alignmentPercentID >= 97 & matches >= 20 & 
                         tNumInsert <= 1 & qNumInsert <= 1 & tBaseInsert <= 1 & qBaseInsert <= 1) %>%
      dplyr::mutate(sample = h$sample[1], hit = h$hits[1])
    
    if(nrow(b) > 0){
      
      b <- left_join(b, barCodes, by = c('qName' = 'readID'))
      b <- distinct(select(b, transcriptBarCode, cellBarCode, tName, strand, qStart, qEnd, tStart, tEnd, qName))
      b$tName <- h$chromosome
      b$tStart <- b$tStart + (h$position - 1000)
      b$tEnd <- b$tEnd + (h$position - 1000)
      b <- left_join(b, tibble(qName = names(o), R1_readSeq = as.character(o)), by = 'qName')
      b$sample <- h$sample
      b$hit <- h$hits
    }
  }
  
  rm(o, barCodes)
  gc()
  
  invisible(file.remove(list.files(file.path(analysisDir, 'tmp'), 
                                   recursive = TRUE, 
                                   pattern = rev(unlist(strsplit(t2, '/')))[1], 
                                   full.names = TRUE)))
  return(b)
}
