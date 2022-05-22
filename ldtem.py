#!/usr/bin/env python
import numpy as np
import numpy.polynomial.chebyshev as nc
import argparse as ap
import numpy.fft as fft
import os,ld,time,sys
import scipy.optimize as so
import time_eph as te
import psr_read as pr
import psr_model as pm
import warnings as wn
import adfunc as af
import matplotlib.pyplot as plt
wn.filterwarnings('ignore')
#
version='JigLu_20201202'
parser=ap.ArgumentParser(prog='ldtem',description='Generate the profile template with multi-profiles.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",nargs='+',help="input ld file or files")
parser.add_argument('-T','--tscrunch',action='store_true',default=False,dest='tscrunch',help='time scrunch to one subint for each file')
parser.add_argument('-r','--frequency_range',default=0,dest='freqrange',help='calculate in the frequency range (FREQ0,FREQ1)')
parser.add_argument('-s','--subint_range',default=0,dest='subint_range',help='calculate in the subint range (SUBINT0,SUBINT1)')
parser.add_argument('-z',"--zap",dest="zap_file",default=0,help="file recording zap channels")
parser.add_argument('-o',"--output",dest="output",default="template",help="outputfile name")
parser.add_argument('-d',action='store_true',default=False,dest='dm_corr',help='the progress will not correcting the DM deviation before calculating rotating phase')
parser.add_argument('-l',"--linear_number",dest="lnumber",type=np.int8,default=0,help="the number of frequency-domain points for linear fitting")
parser.add_argument('-b',"--nbin",dest="nbin",type=np.int16,default=256,help="the bin number of output profile")
parser.add_argument('-a',action='store_true',default=False,dest='auto',help='do not discard the low-rms data')
parser.add_argument('-m',action='store_true',default=False,dest='mean',help='use mean profile as template')
args=(parser.parse_args())
command=['ldtem.py']
#
if args.zap_file:
	command.append('-z')
	if not os.path.isfile(args.zap_file):
		parser.error('The zap channel file is invalid.')
	zchan=np.loadtxt(args.zap_file,dtype=np.int32)
	if np.min(zchan)<0:
		parser.error('The zapped channel number is overrange.')
else:
	zchan=np.int32([])
#
if args.subint_range:
	command.append('-s '+args.subint_range)
	sub_start,sub_end=np.int32(args.subint_range.split(','))
	if sub_end>0 and sub_start>sub_end:
		parser.error("Starting subint number larger than ending subint number.")
	elif sub_start<0:
		parser.error("Input subint is overrange.")
#
if args.tscrunch:
	command.append('-T')
#
if args.freqrange:
	command.append('-r '+args.freqrange)
	freq_s,freq_e=np.float64(args.freqrange.split(','))
	if freq_s>freq_e:
		parser.error("Starting frequency larger than ending frequency.")
else:
	freq_s,freq_e=0,1e9
#
filelist=args.filename
filenum=len(filelist)
nsub_new=[]
psrname=''
dm=0
def ld_check(fname,filetype='Ld file',notfirst=True):
	global freq_s,freq_e,psrname,dm
	if not os.path.isfile(fname):
		parser.error(filetype+' name '+fname+' '+'is invalid.')
	try:
		f=ld.ld(filelist[i])
		finfo=f.read_info()
	except:
		parser.error(filetype+' '+fname+' is invalid.')
	if notfirst:
		tmpname=pr.psr(finfo['psr_par']).name
		if psrname!=tmpname:
			parser.error('The pulsar recorded in '+fname+' is different from the template.')
	else:
		psr=pr.psr(finfo['psr_par'])
		psrname=psr.name
		dm=psr.dm
	#
	if 'compressed' in finfo.keys():
		nchan=int(finfo['nchan_new'])
		nperiod=int(finfo['nsub_new'])
	else:
		nchan=int(finfo['nchan'])
		nperiod=int(finfo['nsub'])
	#
	if args.zap_file:
		if np.max(zchan)>=nchan or np.min(zchan)<0: parser.error('The zapped channel number is overrange.')
	#
	if args.subint_range:
		if sub_end<0: sub_end_tmp=nperiod+sub_end
		else: sub_end_tmp=sub_end
		if sub_start>sub_end_tmp: parser.error("Starting subint number larger than ending subint number.")
		elif sub_end_tmp>nperiod or sub_end_tmp<0: parser.error("Input subint is overrange.")
	#
	if args.tscrunch: nsub_new.append(1)
	elif args.subint_range: nsub_new.append(sub_end_tmp-sub_start)
	else: nsub_new.append(nperiod)
	#
	freq_start,freq_end=np.float64(finfo['freq_start']),np.float64(finfo['freq_end'])
	freq=(freq_start+freq_end)/2.0
	channel_width=(freq_end-freq_start)/nchan
	if args.freqrange:
		if freq_s<freq_start or freq_e>freq_end: parser.error("Input frequency is overrange.")
#
sys.stdout=open(os.devnull,'w')
for i in np.arange(filenum):
	ld_check(filelist[i],notfirst=i)
sys.stdout=sys.__stdout__
#
nsub_new=np.array(nsub_new)
#
name=args.output
name_tmp='    '+name
if name_tmp[-3:]=='.ld':
	name=name[:-3]
if os.path.isfile(name+'.ld'):
	parser.error('The name of output file already existed. Please provide a new name.')
#
def shift(y,x):
	fftdata=fft.rfft(y,axis=1)
	tmp=int(args.nbin/2+1)
	ns,nb=fftdata.shape
	if nb>tmp: fftdata=fftdata[:,:tmp]
	elif nb<tmp: fftdata=np.concatenate((fftdata,np.zeros([ns,tmp-nb])),axis=1)
	if x is not int(0):
		fftdata=fftdata*np.exp(x.repeat(tmp).reshape(-1,tmp)*np.arange(tmp)*1j)
	fftr=fft.irfft(fftdata)
	return fftr
#
if args.dm_corr: command.append('-d ')
dm_zone=np.max([0.1,dm])
dm_zone=np.min([0.5,dm_zone])
#
def dmcor(data,freq,rchan,period,output=1):
	data=data[rchan]
	freq=freq[rchan]
	fftdata=fft.rfft(data,axis=1)
	tmp=np.shape(fftdata)[-1]
	const=(1/freq**2*4148.808/period*np.pi*2.0).repeat(tmp).reshape(-1,tmp)*np.arange(tmp)
	ddm,ddmerr=af.dmdet(fftdata,const,0,dm_zone,9,prec=1e-4)
	if output==1:
		const=(1/freq**2*4148.808/period*np.pi*2.0)*ddm
		data=shift(data[rchan],const)
		return data,ddm,ddmerr
	else:
		return ddm,ddmerr
#
if args.lnumber:
	command.append('-l '+str(args.lnumber))
#
if args.auto:
	command.append('-a ')
#
if args.mean:
	command.append('-m ')
#
def lin(x,k):
	return k*x
#
command=' '.join(command)
#
result=np.zeros([nsub_new.sum(),args.nbin])
cumsub=np.append(0,np.cumsum(nsub_new)[:-1])
for k in np.arange(filenum):
	d=ld.ld(filelist[k])
	info=d.read_info()
	#
	if 'compressed' in info.keys():
		nchan=int(info['nchan_new'])
		nbin=int(info['nbin_new'])
		nperiod=int(info['nsub_new'])
		npol=int(info['npol_new'])
	else:
		nchan=int(info['nchan'])
		nbin=int(info['nbin'])
		nperiod=int(info['nsub'])
		npol=int(info['npol'])
	#
	if args.zap_file:
		if 'zchan' in info.keys():
			zchan_tmp=np.array(list(set(map(int,info['zchan'].split(','))).union(zchan)))
	elif 'zchan' in info.keys():
		zchan_tmp=np.int32(info['zchan'].split(','))
	else:
		zchan_tmp=np.int32([])
	#
	if args.subint_range:
		if sub_end<0:
			sub_e=nperiod+sub_end
		else:
			sub_e=sub_end
		sub_s=sub_start
	else:
		sub_s,sub_e=0,nperiod
	#
	freq_start,freq_end=np.float64(info['freq_start']),np.float64(info['freq_end'])
	freq=(freq_start+freq_end)/2.0
	channel_width=(freq_end-freq_start)/nchan
	chanstart,chanend=np.int16(np.round((np.array([max(freq_s,freq_start),min(freq_e,freq_end)])-freq)/channel_width+0.5*nchan))
	#
	nchan_new=chanend-chanstart
	rchan=np.array(list(set(range(chanstart,chanend))-set(list(zchan_tmp))))-chanstart
	if not args.dm_corr:
		data1=d.period_scrunch(sub_s,sub_e)[chanstart:chanend,:,0]
		freq_real=np.linspace(freq_start,freq_end,nchan+1)[:-1]
		if 'best_dm' in info.keys():
			ddm=np.float64(info['best_dm'][0])-np.float64(info['dm'])
		else:
			ddm,ddmerr=dmcor(data1,freq_real,rchan,np.float64(info['period']),output=0)
		const=(1/freq_real**2*4148.808/np.float64(info['period'])*np.pi*2.0)*ddm
	else:
		const=0
	if not args.tscrunch:
		for s in np.arange(nsub_new[k]):
			data=d.read_period(s+sub_s)[chanstart:chanend,:,0]
			data=shift(data,const)
			result[cumsub[k]+s]=data[rchan].mean(0)
	else:
		data=shift(data1,const)[rchan]
		result[cumsub[k]]=data.mean(0)
#
result-=result.mean(1).reshape(-1,1)
result/=np.sqrt((result**2).sum(1)).reshape(-1,1)
jj=(result.sum(1)**2)>=0
result=result[jj]
fres=fft.rfft(result,axis=1)
if not args.auto:
	fabs=np.abs(fres)
	fabss=fabs.sum(0)
	ntmp=int(np.round(fabss.sum()/fabss.max()))
	ntmp=max(ntmp,20)
	tmp=fabs[:,1:ntmp]
	l=(tmp-tmp.mean(0))/tmp.std(0)
	jj=np.abs(l).max(1)<2
	result=result[jj]
	fres=fres[jj]
	tmp=result.max(1)-result.min(1)
	l=(tmp-tmp.mean(0))/tmp.std(0)
#
fabs=np.abs(fres)
if args.lnumber:
	command.append('-l '+str(args.lnumber))
	lnumber=args.lnumber
	if lnumber>100 or lnumber<=2:
		parser.error('The number of frequency-domain points for linear fitting is invalid.')
else:
	fabss=fabs.sum(0)
	lnumber=int(np.round(fabss.sum()/fabss.max()))
	lnumber=max(4,lnumber)
#
knumber=lnumber*2
nsub,nbin=result.shape
dt=np.zeros(nsub)
angerr=np.zeros([nsub,knumber-1])
nb=int(nbin/2+1)
errbinnum=np.min([int(nb/6),20])
for i in np.arange(nsub):
	f=fres[-1]
	tmpnum=np.argmax(fft.irfft(fres[i]*f.conj()))
	d0=np.append(result[i,tmpnum:],result[i,:tmpnum])
	f0=fft.rfft(d0)
	df=f0/f
	ang=np.angle(df)
	sr0=np.std(f0.real[-errbinnum:])
	si0=np.std(f0.imag[-errbinnum:])
	sr=np.std(f.real[-errbinnum:])
	si=np.std(f.imag[-errbinnum:])
	err=np.sqrt((sr**2+si**2)/fabs[-1]**2+(sr0**2+si0**2)/fabs[i]**2)
	popt,pcov=so.curve_fit(lin,np.arange(1,lnumber),ang[1:lnumber],p0=[0.0],sigma=err[1:lnumber])
	dt[i]=popt[0]/(2*np.pi)-tmpnum/((nb-1)*2)
	angerr[i]=err[1:knumber]
#
result[:-1]=shift(result[:-1],-dt[:-1]*2*np.pi)
#
if args.mean:
	prof=result[jj].mean(0)
else:
	fres=fft.rfft(result,axis=1)
	fang=np.angle(fres)[:,1:knumber]
	fang0=fang[0]*1
	fang-=fang0
	fang[fang>np.pi]-=np.pi*2
	fang[fang<-np.pi]+=np.pi*2
	sum0=fang.mean(0).sum()
	def angle(para):
		ang0=para[:(knumber-2)]
		ang0=np.append(ang0,sum0-ang0.sum())
		k=para[(knumber-2):]
		angc=ang0+np.arange(1,knumber)*k.reshape(-1,1)-fang
		return (angc/angerr).reshape(-1)
	#
	p0=np.zeros(knumber-2+nsub)
	res=so.leastsq(angle,x0=p0,full_output=True)
	popt=res[0]
	if not args.auto:
		k=popt[(knumber-2):]*1
		k-=k.mean()
		k/=k.std()
		jj=np.abs(k)<1
	else:
		jj=np.ones(nsub,dtype=np.bool)
	ang1=popt[:(knumber-2)]
	ang1=np.append(ang1,sum0-ang1.sum())
	weight=1/(np.std(fres[:,-errbinnum:].real,1)**2+np.std(fres[:,-errbinnum:].real,1)**2)[jj]
	abs0=(np.abs(fres[jj])*weight.reshape(-1,1)).sum(0)/weight.sum()
	prof0=abs0[1:knumber]*np.exp(1j*(ang1+fang0))
	prof=fft.irfft(np.concatenate(([0],prof0,np.zeros(nb-knumber))))
#
info={'mode':'template','nchan_new':1, 'nsub_new':1, 'nbin_new':nbin, 'npol_new':1, 'file_time':time.strftime('%Y-%m-%dT%H:%M:%S',time.gmtime()), 'pol_type':'I','compressed':True,'history':command,'psr_name':psrname}
do=ld.ld(name)
do.write_shape([1,1,nbin,1])
do.write_chan(prof,0)
do.write_info(info)
#