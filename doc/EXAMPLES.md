## Running BLASR

  Typing 'blasr -h' or 'blasr -help' on the command line will give you a
  list of options.  At the least, provide a fasta, fastq, or bas.h5 file,
  and a genome.

### Some typical use cases

Align reads from reads.bas.h5 to ecoli_K12 genome, and output in SAM format.

    blasr reads.bas.h5  ecoli_K12.fasta --sam --out alignments.sam

Same as above, but with soft clipping

    blasr reads.bas.h5  ecoli_K12.fasta --sam --clipping soft --out alignments.sam

Use multiple threads

    blasr reads.bas.h5  ecoli_K12.fasta --sam --clipping soft --out alignments.sam --nproc 16

Include a larger minimal match, for faster but less sensitive alignments

    blasr reads.bas.h5  ecoli_K12.fasta --sam --clipping soft --minMatch 15

Produce alignments in a pairwise human readable format

    blasr reads.bas.h5  ecoli_K12.fasta -m 0

Use a precomputed suffix array for faster startup

    sawriter hg19.fasta.sa hg19.fasta #First precompute the suffix array
    blasr reads.bas.h5 hg19.fasta --sa hg19.fasta.sa

Use a precomputed BWT-FM index for smaller runtime memory footprint, but slower alignments.

    sa2bwt hg19.fasta hg19.fasta.sa hg19.fasta.bwt
    blasr reads.bas.h5 hg19.fasta --bwt hg19.fasta.bwt

