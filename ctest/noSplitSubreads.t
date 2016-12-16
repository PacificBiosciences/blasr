Set up
  $ . $TESTDIR/setup.sh

Test blasr with --noSplitSubreads 
#  $ rm -rf $OUTDIR/lambda_bax_noSplitSubreads.m4 
#  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta --noSplitSubreads -m 4 --out lambda_bax_noSplitSubreads_tmp.m4 --nproc 15
#  [INFO]* (glob)
#  [INFO]* (glob)
#  $ sort lambda_bax_noSplitSubreads_tmp.m4 > $OUTDIR/lambda_bax_noSplitSubreads.m4
#  $ diff $OUTDIR/lambda_bax_noSplitSubreads.m4 $STDDIR/lambda_bax_noSplitSubreads.m4
# This test takes a long time, use a subset instad. 

  $ rm -rf $OUTDIR/lambda_bax_noSplitSubreads_subset.m4 
  $ $EXEC $DATDIR/lambda_bax.fofn $DATDIR/lambda_ref.fasta --noSplitSubreads -m 4 --out $OUTDIR/lambda_bax_noSplitSubreads_tmp_subset.m4 --nproc 15 --holeNumbers 1--1000 --sa $DATDIR/lambda_ref.sa
  [INFO]* (glob)
  [INFO]* (glob)
  $ sort $OUTDIR/lambda_bax_noSplitSubreads_tmp_subset.m4 > $OUTDIR/lambda_bax_noSplitSubreads_subset.m4
  $ diff $OUTDIR/lambda_bax_noSplitSubreads_subset.m4 $STDDIR/lambda_bax_noSplitSubreads_subset.m4

# Test key command of unrolled resequencing, check bam header and alignments in output
  $ outbam=$OUTDIR/unrolled-4mer.bam
  $ outsam=$OUTDIR/unrolled-4mer.sam
  $ query=$DATDIR/unrolled/m54006_151021_185942.subreadset.xml
  $ ref=$DATDIR/unrolled/All4mer_V2_11_V2_13_V2_15_V2_44_circular_72x_l50256.fasta
  $ stdsam=$STDDIR/unrolled-4mer.sam
  $ rm -rf $outbam $outsam
  $ $EXEC $query $ref --out $outbam --noSplitSubreads --fastMaxInterval --bam
  [INFO]* (glob)
  [INFO]* (glob)
  $ $SAMTOOLS view -h $outbam -o $outsam
  $ diff $outsam $stdsam
  $ grep '@RG' $outsam
  @RG\tID:e6043908* (glob)
  $ grep 'RG:Z:e6043908' $outsam |wc -l
  4


  $ query=$DATDIR/unrolled/m54006_151021_185942.subreads.bam
  $ outbam=$OUTDIR/unrolled-4mer-bam-in.bam
  $ outsam=$OUTDIR/unrolled-4mer-bam-in.sam
  $ rm -rf $outbam $outsam
  $ $EXEC $query $ref --out $outbam --noSplitSubreads --fastMaxInterval --bam
  [INFO]* (glob)
  [INFO]* (glob)
  $ $SAMTOOLS view -h $outbam -o $outsam
  $ grep -v '^@PG' $outsam > $TMP1 && grep -v '^@PG' $stdsam > $TMP2 && diff $TMP1 $TMP2
  $ grep '@RG' $outsam
  @RG\tID:e6043908* (glob)
  $ grep 'RG:Z:e6043908' $outsam |wc -l
  4


