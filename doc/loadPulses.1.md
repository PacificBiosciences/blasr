% LOADPULSES(1) 3ca7fe8
%
% July 2015

# NAME

loadPulses - load sequencing pulse information

# SYNOPSIS

**loadPulses** *basFileName* *cmpFileName* [*options*]

# DESCRIPTION

**loadPulses** loads sequencing pulse information such as inter-pulse distance
or quality information into the cmp.h5 file.
This allows one to analyze kinetic and quality information by alignment column.

# OPTIONS

*basFileName*
:   The input {bas,pls}.h5 or input.fofn.

*cmpFileName*
:   The cmp.h5 file to load pulse information into.

**-metrics** *value*
:   A comma separated list of metrics (with no spaces).
	Valid options are:

	WhenStarted, QualityValue, InsertionQV, MergeQV,
	DeletionQV, DeletionTag, SubstitutionTag, SubstitutionQV, PreBaseFrames,
	IPD, StartFrame, PulseWidth, WidthInFrames, Light, pkmid, ClassifierQV,
	PulseIndex

	Default options are:
	QualityValue, ClassifierQV, StartFrame, PulseWidth,
	WidthInFrames, pkmid, IPD

**-failOnMissingData**
:   Exit if any data fields are missing from the bas.h5 or pls.h5
	input that are required to load a metric. Defualt is a warning.

**-byread**
:   Load pulse information by read rather than buffering metrics.

**-bymetric**
:   Load pulse information by metric rather than by read. This uses
	more memory than **-byread**, but can be faster.

**-maxElements** *value*
:   Set a limit on the size of pls/bas file to buffer in with
	**-bymetric** (default value: maximum int). Use **-byread** if the limit
	is exceeded.

**-maxMemory** *value*
:   Set a limit (in GB) on the memory to buffer data with **-bymetric**
	(default value: 4 GB). Use **-byread** if the limit is exceeded.

**-metaNElements** *value*
:   Set number of elements in meta data cache for reading
	bas/bax/pls.h5 file.

**-rawNElements** *value*
:   Set number of elements in raw data cache for reading the bas/bax/pls.h5
	file.

**-rawChunkSize** *value*
:   Set chunk size of raw data cache for reading the bas/bax/pls.h5 file.

# SEE ALSO
**blasr**(1)
**pls2fasta**(1)
**samFilter**(1)
**samtoh5**(1)
**samtom4**(1)
**sawriter**(1)
**sdpMatcher**(1)
**toAfg**(1)
