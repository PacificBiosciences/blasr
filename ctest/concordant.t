Set up
  $ . $TESTDIR/setup.sh

Test -concordant
  $ rm -rf $OUTDIR/concordant_subset.sam
  $ $EXEC $DATDIR/ecoli_lp.fofn $DATDIR/ecoli_reference.fasta -concordant -sam -out $OUTDIR/concordant_subset.sam -nproc 12 -holeNumbers 1-10000 -sa $DATDIR/ecoli_reference.sa
  [INFO]* (glob)
  [INFO]* (glob)
  $ sed -n 6,110864p $OUTDIR/concordant_subset.sam > $OUTDIR/tmp1 
  $ sort $OUTDIR/tmp1 > $OUTDIR/tmp11
  $ sed -n 6,110864p $STDDIR/$UPDATEDATE/concordant_subset.sam > $OUTDIR/tmp2
  $ sort $OUTDIR/tmp2 > $OUTDIR/tmp22
  $ diff $OUTDIR/tmp11 $OUTDIR/tmp22
  $ rm -rf $OUTDIR/tmp1 $OUTDIR/tmp2 $OUTDIR/tmp11 $OUTDIR/tmp22
#2014_05_28  --> changelist 135254, use MAX_BAND_SIZE to contrain GuidedAlign
#2014_08_21  --> changelist 138516, added YS, YE, ZM tags. 
#2014_08_28  --> changelist 139176, update SAM MD5 
#2014_09_12  --> changelist 140410, changed the default value of '-concordantTemplate' from 'longestsubread' to 'typicalsubread'
#2014_09_17  --> changelist 140573, changed SDPFragment LessThan to make sure blasr compiled with gcc 4.4 and 4.8 can produce identical results. 
#2014_10_16  --> changelist 141378, changed the default value of '-concordantTemplate' from 'typicalsubread' to 'mediansubread'
#2015_03_01  --> changelist 146599, reads from the same movie should have unique readGroupId
#2015_03_28  --> changelist 148101, 148080 updated read group id, 148100 updated TLEN
#2015_04_09  --> changelist 148796, updated read group id
#2015_04_25  --> changelist 149721, update CIGAR string, replace M with X=.
#2015_04_25  --> changelist ?, force refine all concordant alignments

Test -concordant FMR1 case (the 'typical subread' is selected as template for concordant mapping)
  $ FOFN=$DATDIR/FMR1_concordant.fofn
  $ REF=$DATDIR/FMR1_130CGG.fasta
  $ $EXEC $FOFN $REF -concordant -out $OUTDIR/FMR1_zmw_37927.m4 -m 4 -holeNumbers 37927
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/FMR1_zmw_37927.m4 $STDDIR/$UPDATEDATE/FMR1_zmw_37927.m4
