% SAWRITER(1) 3ca7fe8
%
%  July 2015

# NAME

sawriter - generate suffix arrays for nucleotide sequences

# SYNOPSIS

**sawriter** *saOut* *fastaIn* [*fastaIn2* *fastaIn3* ...] [**-blt** *p*]
	[**-4bit**] [**-larsson**|**-manmy**|**-kar**|**-mafe**|**-welter**]

**sawriter** *fastaIn* (writes to fastaIn.sa)

# OPTIONS

**-blt** *p*
:   Build a lookup table on prefixes of length *p*. This speeds
	up lookups considerably (more than the LCP table), but misses matches
	less than *p* when searching.

**-4bit**
:   Read in (one) fasta file as a compressed sequence file.

**-welterweight** *N*
:   use a difference cover of size \fIN\fR for building the suffix array.
	Valid values are 7,32,64,111, and 2281.

## Methods

**-larsson**
:   (default)  Uses the method of Larsson and Sadakane to build the array.

**-mamy**
:   Uses the method of MAnber and MYers to build the array (slower than larsson,
	and produces the same result. This is mainly for double checking
	the correctness of larsson).

**-kark**
:   Use Karkkainen DS3 method for building the suffix array.
	This will probably be slower than larsson, but takes only an extra
	N/(sqrt 3) extra space.

**-mafe**
:   (disabled for now!) Use the lightweight construction algorithm from
	Manzini and Ferragina

**-welter**
:   Use lightweight (sort of light) suffix array construction.
	This is a bit more slow than normal larsson.

# SEE ALSO

**blasr**(1)
**loadPulses**(1)
**pls2fasta**(1)
**samFilter**(1)
**samtoh5**(1)
**samtom4**(1)
**sdpMatcher**(1)
**toAfg**(1)
