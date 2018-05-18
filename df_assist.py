#!/usr/bin/env python
import numpy as np
import os,ld
import numpy.fft as fft
import argparse as ap
try:
	import astropy.io.fits as ps
except:
	import pyfits as ps
#
parser=ap.ArgumentParser(prog='df_assist',description='Assistance program.')
parser.add_argument("filename",nargs='+',help="name of file or filelist")
parser.add_argument("-s","--start",dest='start_chan',default=1,type=int,help="start channel number")
parser.add_argument("-n","--nchan",dest='nchan',default=0,type=int,help="total channel number")
args=(parser.parse_args())
#
name=args.filename[0]
d=ld.ld(name+'.ld')
info=np.load(name+'_df_assist_temp.npy').item()
nchan_new=info['nchan_new']
zchan=info['zchan']
totalbin=info['totalbin']
temp_multi=info['temp_multi']
nbin0=info['nbin0']
nbin01=info['nbin01']
nbin_old=info['nbin_old']
disp=info['disp']
judge0=info['judge0']
judge1=info['judge1']
phasenum=info['phasenum']
phaseres=info['phaseres']
binint=info['binint']
def shift(y,x):
	fftp=fft.rfft(y)
	ffts=fftp*np.exp(x*1j*np.arange(len(fftp)))
	fftr=fft.irfft(ffts)
	return fftr
for k in np.arange(args.start_chan,args.start_chan+args.nchan):
	if k in zchan:
		d.write_chan(np.zeros(totalbin/temp_multi),k)
		continue
	data=np.zeros(nbin01)
	data[:nbin_old]=d.__read_chan0__(k,ndata_chan0=nbin_old)
	ddata=np.float64(shift(data,disp[k])[:nbin0])[judge0:judge1]
	tpdata=np.zeros(totalbin+1)
	tpdata[phasenum]=ddata*(1-phaseres)
	tpdata[phasenum+1]+=ddata*phaseres
	tpdata=tpdata[:-1].reshape(-1,temp_multi).sum(1)/binint
	d.write_chan(tpdata,k)

