Set up
  $ . $TESTDIR/setup.sh

Without -allowAdjacentIndels, adjacent indels should not exist in SAM/BAM CIGAR strings
  $ $EXEC $DATDIR/test_dataset/nofilter.subreadset.xml $DATDIR/ecoli_reference.fasta -bam -out $OUTDIR/noAdjacentIndels.bam -concordant -bestn 1 && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0

  $ $SAMTOOLS view $OUTDIR/noAdjacentIndels.bam |cut -f 6 |sed 's/[0-9]//g' > $TMP1 && echo $?
  0

  $ sed -n 's/ID/ID/g' $TMP1
  $ sed -n 's/DI/DI/g' $TMP1

With -allowAdjacentIndels, adjacent indels may exist in SAM/BAM CIGAR strings
  $ $EXEC $DATDIR/test_dataset/nofilter.subreadset.xml $DATDIR/ecoli_reference.fasta -bam -out $OUTDIR/allowAdjacentIndels.bam -concordant -bestn 1 -allowAdjacentIndels && echo $?
  [INFO]* (glob)
  [INFO]* (glob)
  0

  $ $SAMTOOLS view $OUTDIR/allowAdjacentIndels.bam |cut -f 6 |sed 's/[0-9]//g' > $TMP2 && echo $?
  0

  $ sed -n 's/ID/ID/g' $TMP1
  $ sed -n 's/DI/DI/g' $TMP1

