% SDPMATCHER(1) 3ca7fe8
%
%  July 2015

# NAME

sdpMatcher - nucleotide sequence matching using sparse dynamic programming

# SYNOPSIS
**sdpMatcher** *query* *target* *k* [**-indelRate** *delta*]
	[**-showalign**] [**-printsw**] [**-noRefine**] [**-indel** *i*]
	[**-local**] [**-match** *m*] [**-sdpIndel** *i*]
	
# DESCRIPTION
	
Performs sparse dynamic programming (SDP) between pairs of sequences as they 
are given in two FASTA files, called *query* and *target* for
convenience. *k* is the size of the k-mer used for the SDP algorithm.

# SEE ALSO

**blasr**(1)
**loadPulses**(1)
**pls2fasta**(1)
**samFilter**(1)
**samtoh5**(1)
**samtom4**(1)
**sawriter**(1)
**toAfg**(1)
