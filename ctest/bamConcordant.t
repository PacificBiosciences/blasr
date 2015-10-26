Set up
  $ . $TESTDIR/setup.sh

Test using bam as input, use -concordant
  $ $EXEC $DATDIR/test_bam/tiny_bam.fofn $DATDIR/bamConcordantRef.fasta -bam -concordant -bestn 1 -out $OUTDIR/bamConcordant.bam
  [INFO]* (glob)
  [INFO]* (glob)

Check whether sam out and bam out have identical alignments, not checking qvs
  $ $SAMTOOLS view $OUTDIR/bamConcordant.bam |cut -f 4 
  1
  1
  8??? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)
  86?? (glob)

  $ $EXEC /pbi/dept/secondary/siv/testdata/SA3-RS/lambda/2372215/0007_tiny/Analysis_Results/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.subreads.bam $DATDIR/lambda_ref.fasta -m 4 -concordant -bestn 1 -holeNumbers 17417 -out $OUTDIR/tmp.m4 -V 2 > $OUTDIR/bamConcordant.log
  [INFO]* (glob)
  [INFO]* (glob)

  $ grep "Concordant template" $OUTDIR/bamConcordant.log
  Concordant template subread index: 8, 17417/14708_16595
