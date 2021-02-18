blasr_cram_script = find_program('cram')

# third-party tools
blasr_samtools = find_program([
  'samtools',
  '/mnt/secondary/Smrtpipe/builds/Internal_Mainline_Nightly_LastSuccessfulBuild/analysis/bin/samtools'])

blasr_fast_test_list = [
  ['ecoli', 'FAST'],
  ['fastMaxInterval', 'FAST'],
  ['aggressiveIntervalCut', 'FAST'],
  ['multipart', 'FAST'],
  ['affineAlign', 'FAST'],
  ['bamOut', 'FAST'],
  ['ccsH5', 'FAST'],
  ['filtercriteria', 'FAST'],
  ['m0-5', 'FAST'],
  ['fofn', 'FAST'],
  ['alignScore', 'FAST'],
  ['hitpolicy', 'FAST'],
  ['noSplitSubreads', 'FAST'],
  ['bamIn', 'FAST'],
  ['open_fail', 'FAST'],
  ['verbose', 'FAST'],
  ['deterministic', 'FAST'],
  ['pgc-naive', 'FAST'],
  ['pgc-fasta', 'FAST'],
  ['pgc-concordant', 'FAST'],
  ['pgc-concordant-naive', 'FAST'],
#  ['concordant', 'INTERMEDIATE'],
  ['bug25766', 'INTERMEDIATE'],
  ['holeNumbers', 'INTERMEDIATE'],
#  ['bug25328', 'SLOW'],
#  ['pgc-big', 'SLOW'],
#  ['long-cigar', 'SLOW'],
]

blasr_test_remotedir = '/pbi/dept/secondary/siv/testdata/BlasrTestData/ctest'
blasr_test_testdir = meson.current_build_dir()
blasr_test_outdir = join_paths(blasr_test_testdir, 'out')

foreach i : blasr_fast_test_list
  test(
    'blasr cram test - ' + i[0],
    blasr_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.build_root(), 'blasr-cram-' + i[0] + '.xml'),
      '--verbose'] +
      files(i[0] + '.t'),
    env : [
      'BLASR_EXE=' + blasr_main.full_path(),
      'SAMTOOLS_EXE=' + blasr_samtools.path(),

      'REMOTEDIR=' + blasr_test_remotedir,
      'DATDIR=' + join_paths(blasr_test_remotedir, 'data'),
      'STDDIR=' + join_paths(blasr_test_remotedir, 'stdout'),
      'TMP1=' + join_paths(blasr_test_outdir, '$$.tmp.out'),
      'TMP2=' + join_paths(blasr_test_outdir, '$$.tmp.stdout'),
      'TESTDIR=' + blasr_test_testdir,
      'OUTDIR=' + blasr_test_outdir,
      'UPDATEDATE=2015_11_09'],
    is_parallel : false,
    timeout : 600)
endforeach
