#!/usr/bin/env python
import numpy as np
import numpy.polynomial.chebyshev as nc
import argparse as ap
import numpy.fft as fft
import os,ld,time
import scipy.optimize as so
import time_eph as te
import psr_read as pr
import psr_model as pm
import warnings as wn
import adfunc as af
wn.filterwarnings('ignore')
#
version='JigLu_20201202'
parser=ap.ArgumentParser(prog='ldtoa',description='Get the relative pulse rotating phase and DM of the ld file.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",help="input ld file")
parser.add_argument('-t',dest='template',help="template ld file")
parser.add_argument('-f',dest='fscrunch',default=0,type=np.int16,help="obtain results after frequency scrunch by this factor")
parser.add_argument('-T',action='store_true',default=False,dest='tscrunch',help='time scrunch to one subint to obtain result')
parser.add_argument('-r','--frequency_range',default=0,dest='freqrange',help='calculate in the frequency range (FREQ0,FREQ1)')
parser.add_argument('-s','--subint_range',default=0,dest='subint_range',help='calculate in the subint range (SUBINT0,SUBINT1)')
parser.add_argument("-z","--zap",dest="zap_file",default=0,help="file recording zap channels")
parser.add_argument("-o","--output",dest="output",default="toa",help="outputfile name")
parser.add_argument('-d',action='store_true',default=False,dest='dm_corr',help='correcting the DM deviation before calculating rotating phase')
parser.add_argument("-n","--linear_number",dest="lnumber",type=np.int8,default=10,help="the number of frequency-domain points for linear fitting")
args=(parser.parse_args())
command=['ldtoa.py']
#
if not os.path.isfile(args.filename):
	parser.error('A valid ld file name is required.')
d=ld.ld(args.filename)
info=d.read_info()
#
if not os.path.isfile(args.template):
	parser.error('A valid ld file name is required.')
d0=ld.ld(args.template)
info0=d0.read_info()
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
if 'compressed' in info0.keys():
	nchan0=int(info0['nchan_new'])
	nsub0=int(info0['nsub_new'])
	nbin0=int(info0['nbin_new'])
else:
	nchan0=int(info0['nchan'])
	nsub0=int(info0['nsub'])
	nbin0=int(info0['nbin'])
#
if args.zap_file:
	command.append('-z')
	if not os.path.isfile(args.zap_file):
		parser.error('The zap channel file is invalid.')
	zchan=np.loadtxt(args.zap_file,dtype=np.int32)
	if np.max(zchan)>=nchan or np.min(zchan)<0:
		parser.error('The zapped channel number is overrange.')
	if 'zchan' in info.keys():
		info['zchan']=str(list(set(map(int,info['zchan'].split(','))).update(zchan)))[1:-1]
	else:
		info['zchan']=str(list(zchan))[1:-1]
	#if nchan_new!=nchan:
	#	info.pop('zchan')
else:
	if 'zchan' in info.keys():
		info.pop('zchan')
	zchan=[]
#
if args.subint_range:
	command.append('-s '+args.subint_range)
	sub_start,sub_end=np.int32(args.subint_range.split(','))
	if sub_start>sub_end:
		parser.error("Starting subint number larger than ending subint number.")
	elif sub_start<0 or sub_end>nperiod:
		parser.error("Input subint is overrange.")
else:
	sub_start,sub_end=0,nperiod
#
if args.tscrunch:
	command.append('-T')
	nsub_new=1
else:
	nsub_new=sub_end-sub_start
#
if args.lnumber:
	lnumber=args.lnumber
	if lnumber>100 or lnumber<=2:
		parser.error('The number of frequency-domain points for linear fitting is invalid.')
else:
	lnumber=10
#
freq0,freq1=np.float64(info['freq_start']),np.float64(info['freq_end'])
freq00,freq10=np.float64(info0['freq_start']),np.float64(info0['freq_end'])
if freq00>=freq1 or freq10<=freq0:
	parser.error('The template has different frequency band from the data.')
#poltype,poltype0=info['pol_type'],info0['pol_type']
#if poltype!=poltype0:
#	parser.error('The polarization type of data and template are not consistent.')
freq=(freq00+freq10)/2.0
bandwidth=freq10-freq00
channel_width=bandwidth/nchan0
if args.freqrange:
	command.append('-r '+args.freqrange)
	freq_start,freq_end=np.float64(args.freqrange.split(','))
	chanstart,chanend=np.int16(np.round((np.array([freq_start,freq_end])-freq)/channel_width+0.5*nchan))
	if freq_start>freq_end:
		parser.error("Starting frequency larger than ending frequency.")
	elif freq_start<max(freq00,freq0) or freq_end>min(freq10,freq1):
		parser.error("Input frequency is overrange.")
else:
	freq_start,freq_end=max(freq00,freq0),min(freq10,freq1)
	chanstart,chanend=np.int16(np.round((np.array([freq_start,freq_end])-freq)/channel_width+0.5*nchan0))
#
if args.fscrunch:
	command.append('-f '+str(args.fscrunch))
	if (chanend-chanstart)<args.fscrunch:
		parser.error('The selected template frequency range cannot be scrunched by the given factor.')
	nchan_new=np.int32((chanend-chanstart)//args.fscrunch)
	if (freq_end-freq_start)/nchan_new<(freq1-freq0)/nchan:
		parser.error('The data are so sparse in frequency domain for the given scrunch factor.')
else:
	nchan_new=chanend-chanstart
#
command=' '.join(command)
#
name=args.output
if os.path.isfile(name):
	parser.error('The name of output file already existed. Please provide a new name.')
if len(name)>3:
	if name[-3:]=='.ld':
		name=name[:-3]
#
if 'history' in info.keys():
	if type(info['history'])==list:
		info['history'].append(command)
		info['file_time'].append(time.strftime('%Y-%m-%dT%H:%M:%S',time.gmtime()))
	else:
		info['history']=[info['history'],command]
		info['file_time']=[info['file_time'],time.strftime('%Y-%m-%dT%H:%M:%S',time.gmtime())]
else:
	info['history']=command
	info['file_time']=time.strftime('%Y-%m-%dT%H:%M:%S',time.gmtime())
#
def shift(y,x):
	ffts=y*np.exp(x*1j)
	fftr=fft.irfft(ffts)
	return fftr
#
dm_zone=np.max([0.1,np.float64(info['dm'])/100])
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
		data=shift(fftdata,const*ddm)
		return data,ddm,ddmerr
	else:
		return ddm,ddmerr
#
rchan=np.array(list(set(range(chanstart,chanend))-set(list(zchan))))-chanstart
tpdata=np.zeros([nchan_new,nbin0])
data0=d0.period_scrunch()[chanstart:chanend,:,0]
if args.dm_corr:
	psr=pm.psr_timing(pr.psr(info0['psr_par']),te.times(te.time(np.float64(info0['stt_time'])+np.float64(info0['length'])/86400,0)),freq_start)
	freq_real=np.linspace(freq_start,freq_end,nchan0+1)[:-1]*psr.vchange.mean()
	if 'best_dm' in info0.keys():
		ddm0=np.float64(info0['best_dm'])-np.float64(info0['dm'])
		fftdata0=fft.rfft(data0,axis=1)
		tmp=np.shape(fftdata0)[-1]
		const=(1/freq_real[:-1]**2*4148.808/np.float64(info0['period'])*np.pi*2.0).repeat(tmp).reshape(-1,tmp)*np.arange(tmp)
		data0=shift(fftdata,const*ddm0)
	else:
		data0[rchan],ddm0,ddm0err=dmcor(data0,freq_real,rchan,np.float64(info0['period']))
		print(ddm0,ddm0err,freq_real,rchan,np.float64(info0['period']),psr.vchange.mean())
i_new=0
nchan0=chanend-chanstart
res=nchan0
for i in np.arange(chanstart,chanend):
	if res>nchan_new:
		res-=nchan_new
		if i in zchan: continue
		tpdata[i_new]+=data0[i-chanstart]
	else:
		chan_data=data0[i-chanstart]
		if i in zchan: chan_data*=0
		tpdata[i_new]+=chan_data*(res*1.0/nchan_new)
		i_new+=1
		if i==chanend-1: break
		tpdata[i_new]=chan_data*((nchan_new-res)*1.0/nchan_new)
		res=nchan0-(nchan_new-res)
#
def lin(x,k):
	return k*x
#
def poa(tpdata0,tpdata):
	nb=int(min(nbin0,nbin)//2+1)
	tpdata-=tpdata.mean(1).reshape(-1,1)
	tpdata/=tpdata.max(1).reshape(-1,1)
	tpdata0-=tpdata0.mean(1).reshape(-1,1)
	tpdata0/=tpdata0.max(1).reshape(-1,1)
	f0=fft.rfft(tpdata0,axis=1)[rchan,:nb]
	d0=fft.irfft(f0,axis=1)
	f=fft.rfft(tpdata,axis=1)[rchan,:nb]
	d=fft.irfft(f,axis=1)
	tmpnum=np.argmax(fft.irfft(f0*f.conj(),axis=1).mean(0))
	d0=np.concatenate((d0[:,tmpnum:],d0[:,:tmpnum]),axis=1)
	f0=fft.rfft(d0,axis=1)
	errbinnum=np.min([int(nb/6),20])
	sr0=np.std(f0.real[:,-errbinnum:],1)
	si0=np.std(f0.imag[:,-errbinnum:],1)
	sr=np.std(f.real[:,-errbinnum:],1)
	si=np.std(f.imag[:,-errbinnum:],1)
	df=f0/f
	err=np.sqrt((sr**2+si**2).reshape(-1,1)/np.abs(f)**2+(sr0**2+si0**2).reshape(-1,1)/np.abs(f0)**2)
	ang=np.angle(df)
	ang1=(ang/err**2).sum(0)/(1/err**2).sum(0)
	err1=1/((1/err**2).sum(0))**0.5
	fitnum=lnumber
	popt,pcov=so.curve_fit(lin,np.arange(1,fitnum),ang1[1:fitnum],p0=[0.0],sigma=err1[1:fitnum])
	dt=-popt/(2*np.pi)+tmpnum/((nb-1)*2)
	dterr=pcov[0,0]**0.5/(2*np.pi)
	return [dt,dterr]
#
freq=(freq0+freq1)/2.0
bandwidth=freq1-freq0
channel_width=bandwidth/nchan
if args.freqrange:
	freq_start,freq_end=np.float64(args.freqrange.split(','))
	chanstart,chanend=np.int16(np.round((np.array([freq_start,freq_end])-freq)/channel_width+0.5*nchan))
freq_new=np.linspace(freq_start,freq_end,nchan_new+1)
freq=np.linspace(freq_start,freq_end,nchan0+1)
d1=ld.ld(name+'.ld')
d1.write_shape([1,nsub_new,2,1])
if nsub_new>1:
	p0=np.zeros([nsub_new,2])
	phase0=int(info['phase0'])
	sub_nperiod=int(info['sub_nperiod'])*np.ones(nperiod,dtype=np.float64)
	sub_nperiod[-1]=int(info['sub_nperiod_last'])
	middle=sub_nperiod.cumsum()-sub_nperiod/2
	time0=np.linspace(0,np.float64(info['length']),nperiod)
	psr=pm.psr_timing(pr.psr(info['psr_name']),te.times(te.time(np.float64(info['stt_date'])*np.ones(nperiod,dtype=np.float64),np.float64(info['stt_sec'])+time0)),np.float64(info['freq_end']))
	phase1=psr.phase
	middle_time=np.interp(middle,phase1.date-phase0+phase1.second,time0)[sub_start:sub_end]
	psr1=pm.psr_timing(pr.psr(info['psr_name']),te.times(te.time(np.float64(info['stt_date'])*np.ones(nsub_new,dtype=np.float64),np.float64(info['stt_sec'])+middle_time)),np.float64(info['freq_end']))
	vchange=psr1.vchange
	data1=d.period_scrunch(sub_start,sub_end)[chanstart:chanend,:,0]
	if args.dm_corr:
		if 'best_dm' in info.keys():
			ddm=np.float64(info['best_dm'])-np.float64(info['dm'])
		else:
			ddm,ddmerr=dmcor(data1,freq[:-1]*vchange.mean(),rchan,np.float64(info['period']),output=0)
			print(ddm,ddmerr,freq[:-1]*vchange.mean(),rchan,np.float64(info['period']),vchange)
		print('The relative DM from the template file to data file is '+str(ddm-ddm0))
	for s in np.arange(nsub_new):
		tpdata0=np.zeros([nchan_new,nbin])
		data=d.read_period(s+sub_start)[chanstart:chanend,:,0]
		if args.dm_corr:
			fftdata=fft.rfft(data,axis=1)
			tmp=np.shape(fftdata)[-1]
			const=(1/freq[:-1]**2*4148.808/np.float64(info['period'])*np.pi*2.0).repeat(tmp).reshape(-1,tmp)*np.arange(tmp)
			data=shift(fftdata,const*ddm)
		for i in np.arange(chanend-chanstart):
			if (i+chanstart) in zchan: continue
			chan0,chan1=np.sum(freq[i]-freq_new>0),np.sum(freq[i+1]-freq_new>0)
			if chan0<chan1:
				res=(freq[i+1]-freq_new[chan0])/(freq[i+1]-freq[i])
				if chan0>0:
					tpdata0[chan0-1]+=data[i]*(1-res)
				tpdata0[chan1-1]+=data[i]*res
			elif chan0==chan1:
				tpdata0[chan0-1]+=data[i]
			else:
				parser.error('The unexpected error.')
		p0[s]=poa(tpdata0,tpdata)
		#d1.write_period(np.array([np.float64(info['stt_date']),np.float64(info['stt_sec'])+middle_time[s],p0[s]]),s)
		print(np.float64(info['stt_date']),np.float64(info['stt_sec'])+middle_time[s],*p0[s])
else:
	tpdata0=np.zeros([nchan_new,nbin])
	phase0=int(info['phase0'])
	sub_nperiod=int(info['sub_nperiod'])*np.ones(nperiod,dtype=np.float64)
	sub_nperiod[-1]=int(info['sub_nperiod_last'])
	sub_nperiod_cumsum=sub_nperiod.cumsum()
	if sub_start==0: middle=sub_nperiod_cumsum[sub_end]/2
	else: middle=(sub_nperiod_cumsum[sub_start-1]+sub_nperiod_cumsum[sub_end])/2
	time0=np.linspace(0,np.float64(info['length']),12)
	phase1=pm.psr_timing(pr.psr(info['psr_name']),te.times(te.time(np.float64(info['stt_date'])*np.ones(12,dtype=np.float64),np.float64(info['stt_sec'])+time0)),(np.float64(info['freq_start'])+np.float64(info['freq_end']))/2).phase
	middle_time=np.interp(middle,phase1.date-phase0+phase1.second,time0)
	data=d.period_scrunch()[chanstart:chanend,:,0]
	data,ddm,ddmerr=dmcor(data,freq[:-1]*vchange[s],rchan)
	for i in np.arange(chanend-chanstart):
		if (i+chanstart) in zchan: continue
		chan0,chan1=np.sum(freq[i]-freq_new>0),np.sum(freq[i+1]-freq_new>0)
		if chan0<chan1:
			res=(freq[i+1]-freq_new[chan0])/(freq[i+1]-freq[i])
			if chan0>0:
				tpdata0[chan0-1]+=data[i]*(1-res)
			tpdata0[chan1-1]+=data[i]*res
		elif chan0==chan1:
			tpdata0[chan0-1]+=data[i]
		else:
			parser.error('The unexpected error.')
	p0=poa(tpdata0,tpdata)
	d1.write_period(np.array([np.float64(info['stt_date']),np.float64(info['stt_sec'])+middle_time,p0]),s)
	print(np.float64(info['stt_date']),np.float64(info['stt_sec'])+middle_time,p0)
#
info['mode']='ToA'
d1.write_info(info)

