% SAMTOH5(1) 3ca7fe8
%
% July 2015

# NAME

samtoh5 - convert a SAM file to cmp.h5 format

# SYNOPSIS

**samtoh5** *in.sam* *reference.fasta* *out.cmp.h5* [*options*]

# OPTIONS

*in.sam*
:   Input SAM file.

*reference.fasta*
:   Reference used to generate reads.

*out.cmp.h5*
:   Output cmp.h5 file.

**-smrtTitle**
:   Use this option when converting alignments generated from reads
	produced by the **pls2fasta**(1) from bas.h5 files by parsing read
	coordinates from the SMRT read title. The title is in the
	format */name/hole/coordinates*, where *coordinates* are in the
	format `\d+_\d+` and represent the interval of the read that was
	aligned.

**-readType** *value*
:   Set the read type: 'standard', 'strobe', 'CCS', or 'cDNA'

**-verbosity** *value*
:   Set desired verbosity.

**-useShortRefName**
:   Use abbreviated reference names obtained from *file.sam* instead
	of using full names from *reference.fasta*.

**-copyQVs**
:   Copy all QVs available in the SAM file into the cmp.h5 file.
	This includes things like InsertionQV and DeletionTag.

# NOTES

Because SAM has optional tags that have different meanings in different
programs, careful usage is required in order to have proper output. The
"xs" tag in bwa-sw is used to show the suboptimal score, but in PacBio SAM
(**blasr**(1)), it is defined as the start in the query sequence of the
alignment. When **-smrtTitle** is specified, the xs tag is ignored, but when
it is not specified, the coordinates given by the xs and xe tags are used to
define the interval of a read that is aligned. The CIGAR string is relative to
this interval.

# SEE ALSO

**blasr**(1)
**loadPulses**(1)
**pls2fasta**(1)
**samFilter**(1)
**samtom4**(1)
**sawriter**(1)
**sdpMatcher**(1)
**toAfg**(1)
