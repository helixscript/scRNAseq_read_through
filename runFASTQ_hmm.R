#!/usr/bin/Rscript
library(dplyr)
library(optparse)
library(Biostrings)
library(ShortRead)
library(parallel)

CPUs <- 15
hmmFile <- 'ELPS_CD19_BBzeta.hmm'
minHMMscore <- 5


tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')) }

if(! dir.exists('output')) dir.create('output')

write(date(), 'log', append = FALSE)

dirs <- list.dirs('data', recursive = TRUE, full.names = TRUE)

for(x in dirs){
  if(x == 'data') next
  
  sampleName <- unlist(strsplit(x, '/'))[2]
  dir.create(file.path('output', sampleName))
  
  for(x2 in list.files(x, pattern = 'R1', full.names = TRUE)){

    partName <- stringr::str_extract(x2, 'part\\d+')
    dir.create(file.path('output', sampleName, partName))
    
    scriptPath <- file.path('output', sampleName, partName, 'script.sh')
    
    write(c('#!/usr/bin/sh', 
            paste0('./runLTRhmm.R', 
                   ' --R1 ', x2, 
                   ' --R2 ', sub('R1_', 'R2_', x2),
                   ' --outputDir ', file.path('output', sampleName, partName), 
                   ' --CPUs ', CPUs, 
                   ' --hmmFile ', hmmFile, 
                   ' --minHMMscore ', minHMMscore)), 
          file = scriptPath)
    
    system(paste0('chmod 755 ./', scriptPath))
    system(scriptPath, wait = TRUE)
  }
}