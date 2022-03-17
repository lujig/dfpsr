#!/usr/bin/env python
import numpy as np
import time_eph as te
import os,sys,time,ld
import scipy.signal as ss
import argparse as ap
#
version='JigLu_20220317'
#
parser=ap.ArgumentParser(prog='update_tdiff',description='Updating the time difference between local clock and gps.',epilog='Ver '+version)
parser.add_argument('-v','--version', action='version', version=version)
parser.add_argument('--verbose', action="store_true",default=False,help="print detailed information")
parser.add_argument("-w","--rewrite",action="store_true",default=False,help="update the whole clock difference file.")
args=(parser.parse_args())
#
dirname=os.path.split(os.path.realpath(__file__))[0]
outfile=dirname+'/conventions/local2gps.txt'
if args.rewrite:
	d0=open(outfile,'w+')
	d0.close()
	endunix=0
else:
	d0=np.loadtxt(outfile)
	endunix=d0[-1,0]
	end=te.time(endunix,0,scale='unix',unit=1).unix2local().mjd[0]
	endyr,endmo,endday,_,_=te.mjd2datetime(end)
#
str0='FAST_CLOCKDIFF-'
l0=len(str0)
deltat=50
flist=os.listdir(dirname+'/clock')
t0=[]
dt0=[]
for i in flist:
	if i[:l0]!=str0: continue
	date=i[l0:-4]
	if not args.rewrite:
		if int(date[:4])<endyr: continue
		if len(date)==6:
			if int(date[:4])==endyr and int(date[-2:])<endmo: continue
	tmp=np.genfromtxt(dirname+'/clock/'+i,skip_header=21)
	if not args.rewrite:
		if tmp[-1,0]-endunix<=50: continue
	t0.extend(np.int64(tmp[:,0]))
	dt0.extend(np.float64(tmp[:,1]))
#
if len(t0)<20:	print('All present clock data have been included in the clock-difference file.'); exit()
dt0=np.array(dt0)[np.argsort(t0)]
t0=np.sort(t0)
jt0=(t0>(endunix-deltat))&(dt0>0)
t1=t0[jt0]
dt1=dt0[jt0]
if len(t1)<20:	print('All present clock data have been included in the clock-difference file.'); exit()
lt=len(t1)
t1a=t1[:int(lt//10*10)].reshape(-1,10).mean(1)
dt1a=dt1[:int(lt//10*10)].reshape(-1,10)
jj=(np.max(dt1a,1)-np.min(dt1a,1))<1e-5
#dt1b=np.polyval(np.polyfit(t1a[jj],dt1a[jj].mean(1),3),t1)
dt1b=np.interp(t1,t1a[jj],dt1a[jj].mean(1))
jt1=np.abs(dt1b-dt1)<1e-5
t2=t1[jt1]
dt2=dt1[jt1]
if len(t2)<20:	print('All present clock data have been included in the clock-difference file.'); exit()
if endunix==0: sttunix=t2[0]
else: sttunix=endunix
t3=np.arange(sttunix,t2[-1],10)
if len(t3)<20:	print('All present clock data have been included in the clock-difference file.'); exit()
dt3=np.interp(t3,t2,dt2)
b,a=ss.butter(15,0.15,'lowpass')
dt3a=ss.filtfilt(b,a,dt3)
t4=np.arange(sttunix+deltat,t3[-1],deltat)
dt4=dt3a[int(deltat//10)::int(deltat//10)]
f=open(outfile,'a+')
for i in range(len(t4)):
	f.write(str(int(t4[i]))+'     '+'{:.12f}'.format(dt4[i])+'\n')
f.close()

