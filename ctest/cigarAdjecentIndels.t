Set up
  $ . $TESTDIR/setup.sh

Without -allowAdjacentIndels, adjacent indels should not exist in SAM/BAM CIGAR strings
  $ $EXEC $DATDIR/test_dataset/nofilter.subreadset.xml $DATDIR/ecoli_reference.fasta -bam -out $OUTDIR/noAdjacentIndels.bam -concordant -bestn 1 && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0

  $ $SAMTOOLS view $OUTDIR/noAdjacentIndels.bam |cut -f 6 > $TMP1

  $ grep 'ID' $TMP1 |wc -l
  0

  $ grep 'DI' $TMP1 |wc -l
  0

With -allowAdjacentIndels
  $ $EXEC $DATDIR/test_dataset/nofilter.subreadset.xml $DATDIR/ecoli_reference.fasta -bam -out $OUTDIR/allowAdjacentIndels.bam -concordant -bestn 1 -allowAdjacentIndels && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0
