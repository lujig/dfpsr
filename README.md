# dfpsr
dedisperse and fold search mode psrfits data

dfpsr.py: 

	Dedisperse and fold the search mode psrfits file with pulsar name (or PSRCAT ephemris file, or DM and period). The folded file is saved as ld format, which is a new format to save pulsar data and information. The data in ld file always have 4 dimensions (nchan, nsub, nbin, npol).

ld.py:

	Provide some functions to access ld format data. With these functions, one can read data and information of ld file, or write data and information in a ld file.

compress.py:

	Compress the ld format file with given nchan, nsub, nbin, and save resutls in a new ld file.

para.py:

	View the information of ld format file.

plot.py:

	Plot the time-domain or frequency-domain image or pulse profile of a ld file.

Independence: 

	Software: PSRCAT and TEMPO2

	Python module: numpy, matplotlib, astropy (or pyfits) and Tkinter (or gtk)
