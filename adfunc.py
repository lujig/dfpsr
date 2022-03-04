import numpy as np
import numpy.fft as fft
import scipy.optimize as so
#
def shift(y,x):
	ffts=y*np.exp(x*1j)
	fftr=fft.irfft(ffts)
	return fftr
#
def dmdet(fftdata,dmconst,dm0,dmw,polynum,prec=0):
	length=100
	dm=np.linspace(dm0-dmw,dm0+dmw,length)
	value=np.zeros(length)
	for i in np.arange(length):
		disp=dm[i]*dmconst
		value[i]=np.max(shift(fftdata,disp).sum(0))
		value[i]=(shift(fftdata,disp).sum(0)**2).sum()
	polyfunc=np.polyfit(dm,value,polynum)
	fitvalue=np.polyval(polyfunc,dm)
	roots=np.roots(np.polyder(polyfunc))
	roots=np.real(roots[np.isreal(roots)])
	if len(roots):
		dmmax=roots[np.argmin(np.abs(roots-dm[np.argmax(value)]))]
		if dmmax<dm[-1] and dmmax>dm[0] and np.polyval(np.polyder(polyfunc,2),dmmax)<0:
			res=value-fitvalue
			error=np.std(res)
			errfunc=np.append(polyfunc[:-1],polyfunc[-1]-np.polyval(polyfunc,dmmax)+error)
			roots=np.roots(errfunc)
			roots=np.real(roots[np.isreal(roots)])
			dmerr=np.mean(np.abs([roots[np.argmin(np.abs(roots-dmmax+dmw/10))],roots[np.argmin(np.abs(roots-dmmax-dmw/10))]]-dmmax))
			error1=np.std(res[1:]-res[:-1])
			if prec>0:
				if dmerr<prec or error1>error: return dmmax,dmerr
				else: return dmdet(fftdata,dmconst,dmmax,dmerr*20,polynum,prec=prec)
			else: return dmmax,dmerr,dm,value,fitvalue
	return 0,0


