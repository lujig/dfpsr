#!/usr/bin/env python
import urllib.request as ur
import time_eph as rd
import numpy as np
import scipy.signal as ss
import os
#
dirname=os.path.split(os.path.realpath(__file__))[0]
#
# utc2tai, leap second
a=ur.urlopen('https://hpiers.obspm.fr/eop-pc/earthor/utc/UTC-offsets_tab.html').readlines()
b=list(map(lambda x:x.strip().decode(),a))
lenb=len(b)
for i in range(lenb):
	num=lenb-1-i
	t=b[num]
	if t:
		if not t[0:4].isdigit(): 
			b.pop(num)
		elif int(t[0:4])<1972:
			b.pop(num)
		else:
			b[num]=t.replace('Jan.','1').replace('Jul.','7').replace('1s','').strip()+'1\n'
			if b[num].split()[:2]==['1972','1']:
				b.pop(num)
	else:
		b.pop(num)
f=open(dirname+'/conventions/'+'leap.txt','w')
f.writelines(b)
f.close()
#
# utc2ut1
a=ur.urlopen('https://hpiers.obspm.fr/eoppc/eop/eopc04/eopc04.62-now').readlines()
b=list(map(lambda x:x.strip().decode(),a))
lenb=len(b)
for i in range(lenb):
	num=lenb-1-i
	t=b[num]
	if t:
		if not t[0:4].isdigit(): 
			b.pop(num)
		else:
			b[num]=t+'\n'
	else:
		b.pop(num)
f=open(dirname+'/conventions/'+'eopc.txt','w')
f.writelines(b)
f.close()
#
# tai2ut1
mjd,deltat=np.float64(list(map(lambda x:x.split(),b)))[:,[3,6]].T
taimjd=(rd.time(mjd,np.zeros_like(mjd),scale='utc')).tai().mjd
f=np.loadtxt(dirname+'/conventions/'+'leap.txt')
leap_time0=np.array(list(map(lambda x:rd.datetime2mjd([x[0],x[1],x[2],0,0,0]).mjd,f))).reshape(-1)
leap_time=np.array(list(map(lambda x:f[leap_time0<=x,3].sum(),mjd)))
deltat+=leap_time-10
f=open(dirname+'/conventions/'+'tai2ut1.txt','w')
f.writelines(list(map(lambda x: str(x[0])+' '+str(x[1])+' '+str(x[2])+'\n',zip(mjd,taimjd,deltat))))
f.close()
#
# tai2tt
a=ur.urlopen('ftp://ftp2.bipm.org/pub/tai/ttbipm/TTBIPM.2019').readlines()
b=list(map(lambda x:x.decode(),a))
lenb=len(b)
for i in range(lenb):
	num=lenb-1-i
	t=b[num]
	if t:
		if not t[0:5].isdigit(): 
			b.pop(num)
	else:
		b.pop(num)
f=open(dirname+'/conventions/'+'tai2tt.txt','w')
f.writelines(b)
f.close()
#
# gps2utc
a=ur.urlopen('ftp://ftp2.bipm.org/pub/tai/other-products/utcgnss/utc-gnss').readlines()[20:]
b=list(map(lambda x:x.decode().strip(),a))
lenb=len(b)
for i in range(lenb):
	num=lenb-1-i
	t=b[num]
	if t:
		if not t[0:5].isdigit(): 
			b.pop(num)
		else:
			b[num]=t[:23]+'\n'
	else:
		b.pop(num)
f=open(dirname+'/conventions/'+'gps2utc.txt','w')
f.writelines(b)
f.close()
#
# local2gps
a=ur.urlopen('https://crafts.bao.ac.cn/pub/fast/time/fast-hm-gps-diff.txt').readlines()
b=np.float64(list(map(lambda x:x.decode().strip().split(),a)))
lenb=len(b)
t0,dt0=b.T
t1=np.arange(t0[0],t0[-1]+1000,1000)
dt1=np.interp(t1,t0,dt0)
b0=ss.firwin(101,0.005)
dt2=ss.lfilter(b0,1,dt1) 
f=open(dirname+'/conventions/'+'local2gps.txt','w')
for i in range(lenb-100):
	f.write(str(int(t1[i+50]))+'     '+'{:.12f}'.format(dt2[i+100])+'\n')
f.close()

