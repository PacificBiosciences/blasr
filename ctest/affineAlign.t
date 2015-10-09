Set up
  $ . $TESTDIR/setup.sh

Test affineAlign
  $ rm -rf $OUTDIR/affineAlign.m0
  $ $EXEC $DATDIR/affineAlign.fofn $DATDIR/substr_with_ins.fasta -m 0 -out $OUTDIR/affineAlign.m0  -affineAlign  -holeNumbers 493 -insertion 100 -deletion 100
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/affineAlign.m0 $STDDIR/affineAlign_2014_06_10.m0

  $ rm -rf $OUTDIR/ecoli_affine.m0
  $ $EXEC $DATDIR/ecoli_affine.fasta $DATDIR/ecoli_reference.fasta -m 0 -out $OUTDIR/ecoli_affine.m0 -affineAlign -insertion 100 -deletion 100
  [INFO]* (glob)
  [INFO]* (glob)
  $ diff $OUTDIR/ecoli_affine.m0 $STDDIR/ecoli_affine_2014_06_10.m0

# Note that MapQV for -affineAlign has been fixed in 2014 04 18, bug 24363 
