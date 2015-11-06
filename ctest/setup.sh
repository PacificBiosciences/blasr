# Set up directories
CURDIR=$TESTDIR
REMOTEDIR=/pbi/dept/secondary/siv/testdata/BlasrTestData/ctest
DATDIR=$REMOTEDIR/data
OUTDIR=$CURDIR/out
STDDIR=$REMOTEDIR/stdout

# Set up the executable: blasr.
EXEC=${BLASR_PATH}/blasr

# Define tmporary files
TMP1=$OUTDIR/$$.tmp.out
TMP2=$OUTDIR/$$.tmp.stdout

# Make OUTDIR
mkdir -p $OUTDIR

#FIXME: make samtools independent of absolute build path.
SAMTOOLS=/mnt/secondary/Smrtpipe/builds/Internal_Mainline_Nightly_LastSuccessfulBuild/analysis/bin/samtools

#Update date
UPDATEDATE=2015_11_05

# 2014_08_21 --> change 138516: added YS, YE, ZM tags
# 2014_08_28 --> change 139176: Update SAM MD5 
# 2015_03_28 --> change 148101: 148080 update read group id, 148100 update TLEN. 
# 2015_04_09 --> change 148796: update read group id
# 2015_04_25 --> change 149721, update CIGAR string, replace M with X=
# 2015_04_26 --> change 149749, add opiton -cigarUseSeqMatch (default: false). If -cigarUseSeqMatch is turned on, CIGAR strings use '=' and 'X' to represent sequence match and mismatch instead of 'M'.
# 2015_11_05 --> change 166177, update CIGAR string, DO NOT allow adjacent indels unless -allowAdjacentIndels is ON.
