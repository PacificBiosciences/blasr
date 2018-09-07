Test input query is not pb bam, output is sam/bam, alignment has more than 65535 cigar ops.

Set up
  $ mkdir -p $OUTDIR

  $ Q=$DATDIR/long-cigar/query.fasta
  $ T=$DATDIR/long-cigar/ref.fasta

Test m4 output
  $ O=$OUTDIR/long-cigar.sam
  $ $BLASR_EXE $Q $T --minMatch 15 --maxMatch 25 --advanceHalf --advanceExactMatches 10  --bestn 1 --nproc 1 --noSplitSubreads --out $O 1>/dev/null 2>/dev/null && echo $?
  0
  $ grep 'CG:' $O |wc -l
  1
