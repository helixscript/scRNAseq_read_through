# runFASTQ_HMM.R 
This script locates the FASTQ data chunks in the data directory and creates system calls to runLTRhmm.R which runs an HMM created 
from the last 100 NT of the expected LTR against the reads in each data chunk. Significant hits are written to the output folder 
and latter collated. This script only analyzes R1 though both R1 and R2 read sequences are stored in the output so that you do not 
have to search through the massive dataset afterwards. This script takes a couple of days to run using 25 CPUs.  

# assemble_align_resolveMultiHits.R
This script collates all the HMM hits from runFASTQ_HMM.R and filters results, aligns reads to the human genome, and attempts to 
resolve multi-hit R1 alignments by considering R2 reads that align nearby in the opposite direction.
