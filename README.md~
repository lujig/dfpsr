# dfpsr
dedisperse and fold search mode psrfits data

dfpsr.py: 
	Dedisperse and fold the search mode psrfits file with pulsar name (or PSRCAT ephemris file, or DM and period). The folded file is saved as ld format, which is a new format to save pulsar data and information. The data in ld file always have 4 dimensions (nchan, nsub, nbin, npol).

	dfpsr.py [-f FREQRANGE] [-d DM] [-p PERIOD] [-n PSR_NAME] [-e PAR_FILE] [-b NBIN] filename [filename ...]

ld.py:
	Provide some functions to access ld format data. With these functions, one can read data and information of ld file, or write data and information in a ld file.

compress.py:
	Compress the ld format file with given nchan, nsub, nbin, and save resutls in a new ld file.

	compress.py [-f NCHAN] [-F] [-t NSUB] [-T] [-b NBIN] [-B] [-r FREQU_RANGE] [-s SUBINT_RANGE] filename

para.py:
	View the information of ld format file.

	para.py [-c PARAMETER_NAME_LIST] filename

plot.py:
	Plot the time-domain or frequency-domain image or pulse profile of a ld file.

	plot [-f] [-t] [-p] [-b PHASE_RANGE] [-r FREQ_RANGE] [-s SUBINT_RANGE] filename

Independence: 

	Software: PSRCAT and TEMPO2

	Python module: numpy, matplotlib, astropy (or pyfits) and Tkinter (or gtk)
