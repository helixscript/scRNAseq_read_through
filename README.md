# runFASTQ_HMM.R 
This script locates the FASTQ data chunks in the data directory and creates system calls to runLTRhmm.R which runs an HMM created 
from the last 100 NT of the expected LTR against the reads in each data chunk. Significant hits are written to the output folder 
and latter collated. This script only analyzes R1 though both R1 and R2 read sequences are stored in the output so that you do not 
have to search through the massive dataset afterwards. This script takes a couple of days to run using 25 CPUs.  

# assemble_and filter_HMM_result.R
This script collates all the HMM hits from runFASTQ_HMM.R and filters results. HMM hits are retained if HMM alignments include up to 
at least position 98 in the HMM which is immediately before the terminal CA. R1 reads must include at least a weak match to the expected 
TSO oligo sequence in the expected read position. Cell barcodes are corrected to match the 10x expected barcode list if they deviate by 
1 mismatch. 

# align_HMM_hits.R
This script isolates the NTs following HMM hits and aligns them to the human genome including the corresponding R2 read sequences.

# process_alignments.R
This script analyzes the blat alignments and looks for instances where R2 aligned near R1 in order to add additional support to the call. 
