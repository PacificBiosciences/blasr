% PLS2FASTA(1) 3ca7fe8
%
% July 2015

# NAME

pls2fasta - convert plx.h5/bax.h5/fofn files to fasta or fastq files

# SYNOPSIS

**pls2fasta** *in.bax.h5* *out.fasta* [**options**]

# DESCRIPTION

Although fasta files are provided with every run, they are not trimmed nor
split into subreads. This program takes additional annotation
information, such as the subread coordinates and high quality regions,
and uses them to create fasta sequences that are substrings of all
bases called. Most of the time, you will want to trim low quality
reads, so you should specify **-trimByRegion**.

# OPTIONS

*in.bax.h5*
:   Input plx.h5/bax.h5/fofn file.

*out.fasta*
:   Output fasta/fastq file.

**-trimByRegion**
:   Trim away low quality regions.

**-maskByRegion**
:   Mask low quality regions with 'N'.

**-regionTable** *value*
:   Optional HDF file with a /PulseData/Regions dataset.

**-minSubreadLength** *value*
:   Do not write subreads less than the specified length.

**-noSplitSubreads**
:   Do not split reads on adapter sequences.

**-holeNumber**
:   Only print this hole number (or list of numbers).

**-fastq**
:   Print in FASTQ format with quality.

**-ccs**
:   Print de novo circular consensus (ccs) sequences

**-lineLength** *value*
:   Specify fasta/fastq line length

**-minReadScore** *value*
:   Minimum read score to print a read.
	The score is a number between 0 and 1000 and represents the expected
	accuracy percentage \* 10. A typical value would be between 750 and 800.
	This does not apply to ccs reads.

**-best**
:   If a ccs sequence exists, print this. Otherwise, print the longest subread.
	This does not support fastq.

# SEE ALSO

**blasr**(1)
**loadPulses**(1)
**samFilter**(1)
**samtoh5**(1)
**samtom4**(1)
**sawriter**(1)
**sdpMatcher**(1)
**toAfg**(1)
