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
parser=ap.ArgumentParser(prog='ldtoa',description='Get the relative pulse rotating phase and DM of the ld file.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",nargs='+',help="input ld file or files")
parser.add_argument('-t',dest='template',required=True,help="template ld file")
parser.add_argument('-f',action='store_true',default=False,dest='freq_align',help="use same frequency band to obtain the ToA")
parser.add_argument('-T','--tscrunch',action='store_true',default=False,dest='tscrunch',help='time scrunch to one subint to obtain result')
parser.add_argument('-r','--frequency_range',default=0,dest='freqrange',help='calculate in the frequency range (FREQ0,FREQ1)')
parser.add_argument('-s','--subint_range',default=0,dest='subint_range',help='calculate in the subint range (SUBINT0,SUBINT1)')
parser.add_argument('-z',"--zap",dest="zap_file",default=0,help="file recording zap channels")
parser.add_argument('-Z',action='store_true',default=False,dest="zap_template",help="zap same channels for the template file")
parser.add_argument('-o',"--output",dest="output",default="",help="outputfile name")
parser.add_argument('-d',action='store_true',default=False,dest='dm_corr',help='the progress will not correcting the DM deviation before calculating rotating phase')
parser.add_argument('-l',"--linear_number",dest="lnumber",type=np.int8,default=0,help="the number of frequency-domain points for linear fitting")
parser.add_argument('-n',action='store_true',default=False,dest='norm',help='normalized the data at each channel before cal')
args=(parser.parse_args())
command=['ldtoa.py']
sys.stdout=open(os.devnull,'w')
#
if not os.path.isfile(args.template):
	parser.error('A valid ld file name is required.')
#
d0=ld.ld(args.template)
info0=d0.read_info()
psrname=pr.psr(info0['psr_par']).name
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
	if np.min(zchan)<0:
		parser.error('The zapped channel number is overrange.')
	if args.zap_template:
		command.append('-Z')
		if np.max(zchan)>=nchan0:
			parser.error('The zapped channel number is overrange.')
		if 'zchan' in info0.keys():
			info0['zchan']=str(list(set(map(int,info0['zchan'].split(','))).union(zchan)))[1:-1]
		else:
			info0['zchan']=str(list(zchan))[1:-1]
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
if args.freq_align:
	command.append('-f')
#
freq_start0,freq_end0=np.float64(info0['freq_start']),np.float64(info0['freq_end'])
freq0=(freq_start0+freq_end0)/2.0
channel_width0=(freq_end0-freq_start0)/nchan0
#
if args.freqrange:
	command.append('-r '+args.freqrange)
	freq_s,freq_e=np.float64(args.freqrange.split(','))
	if freq_s>freq_e:
		parser.error("Starting frequency larger than ending frequency.")
else:
	freq_s,freq_e=freq_start0,freq_end0
#
filelist=args.filename
filenum=len(filelist)
file_times=[]
file_len=[]
nsub_new=[]
def ld_check(fname,filetype='Ld file'):
	global freq_s,freq_e
	if not os.path.isfile(fname):
		parser.error(filetype+' name '+fname+' '+'is invalid.')
	try:
		f=ld.ld(filelist[i])
		finfo=f.read_info()
	except:
		parser.error(filetype+' '+fname+' is invalid.')
	tmpname=pr.psr(finfo['psr_par']).name
	if psrname!=tmpname:
		parser.error('The pulsar recorded in '+fname+' is different from the template.')
	#
	if 'compressed' in finfo.keys():
		nchan=int(finfo['nchan_new'])
		nperiod=int(finfo['nsub_new'])
	else:
		nchan=int(finfo['nchan'])
		nperiod=int(finfo['nsub'])
	#
	if args.zap_file:
		if args.zap_template:
			if nchan0!=nchan: parser.error('The channel numbers of data and template are different, and the zap file should not be same.')
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
	if freq_start0>=freq_end or freq_end0<=freq_start: parser.error('The template has different frequency band from the data.')
	freq=(freq_start+freq_end)/2.0
	channel_width=(freq_end-freq_start)/nchan
	if args.zap_file and args.zap_template:
		if max(np.abs(freq_end-freq_end0),np.abs(freq_start-freq_start0))>min(channel_width,channel_width0):
			parser.error('The frequency ranges of data and template are different, and the zap file should not be same.')
	if args.freq_align:
		if freq_start0>freq_end or freq_end<freq_start0: parser.error("Frequency bands of data and template have no overlap.")
		if args.freqrange: 
			if freq_s<freq_start or freq_e>freq_end: parser.error("Input frequency is overrange.")
		else:
			freq_s=max(freq_s,freq_start)
			freq_e=min(freq_e,freq_end)
	elif args.freqrange:
		if freq_s<freq_start or freq_e>freq_end: parser.error("Input frequency is overrange.")
#
for i in np.arange(filenum):
	ld_check(filelist[i])
#
nsub_new=np.array(nsub_new)
if args.freq_align:
	chanstart0,chanend0=np.int16(np.round((np.array([freq_s,freq_e])-freq0)/channel_width0+0.5*nchan0))
else:
	chanstart0,chanend0=0,nchan0
#
name=args.output
if not name:
	output='screen'
else:
	name_tmp='    '+name
	if name_tmp[-3:]=='.ld':
		output='ld'
		name=name[:-3]
	elif name_tmp[-4:]=='.txt': 
		if os.path.isfile(name):
			parser.error('The name of output file already existed. Please provide a new name.')
		output='txt'
	else: output='ld'
	if output=='ld':
		if os.path.isfile(name+'.ld'):
			parser.error('The name of output file already existed. Please provide a new name.')
#
def shift(y,x):
	ffts=y*np.exp(x*1j)
	fftr=fft.irfft(ffts)
	return fftr
#
dm_zone=np.max([0.1,np.float64(info0['dm'])/100])
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
if args.zap_file and args.zap_template:
	rchan=np.array(list(set(range(chanstart0,chanend0))-set(list(zchan))))-chanstart0
else:
	rchan=np.arange(chanend0-chanstart0)
data0=d0.period_scrunch()[chanstart0:chanend0,:,0]
if not args.dm_corr and chanend0-chanstart0>1:
	psr=pm.psr_timing(pr.psr(info0['psr_par']),te.times(te.time(np.float64(info0['stt_time'])+np.float64(info0['length'])/86400,0)),freq_start0)
	freq_real0=np.linspace(freq_start0,freq_end0,nchan0+1)[:-1]*psr.vchange.mean()
	if 'best_dm' in info0.keys():
		ddm0=np.float64(info0['best_dm'][0])-np.float64(info0['dm'])
		fftdata0=fft.rfft(data0,axis=1)
		tmp=np.shape(fftdata0)[-1]
		const=(1/freq_real0**2*4148.808/np.float64(info0['period'])*np.pi*2.0).repeat(tmp).reshape(-1,tmp)*np.arange(tmp)
		data0=shift(fftdata0,const*ddm0)
	else:
		data0[rchan],ddm0,ddm0err=dmcor(data0,freq_real0,rchan,np.float64(info0['period']))
data_tmp=data0[rchan]
if args.norm:
	command.append('-n')
	data_tmp-=data_tmp.mean(1).reshape(-1,1)
	data_tmp/=data_tmp.std(1).reshape(-1,1)
tpdata0=data0[rchan].sum(0)
#
if args.lnumber:
	command.append('-l '+str(args.lnumber))
	lnumber=args.lnumber
	if lnumber>100 or lnumber<=2:
		parser.error('The number of frequency-domain points for linear fitting is invalid.')
else:
	tmp=tpdata0-af.baseline(tpdata0)
	ew=tmp.sum()/tmp.max()
	lnumber=int(round(nbin0/ew/2))
	lnumber=max(3,lnumber)
#
def lin(x,k):
	return k*x
#
def poa(tpdata0,tpdata):
	nb=int(min(nbin0,nbin)//2+1)
	tpdata-=tpdata.mean()
	tpdata/=tpdata.max()
	tpdata0-=tpdata0.mean()
	tpdata0/=tpdata0.max()
	f0=fft.rfft(tpdata0)[:nb]
	d0=fft.irfft(f0)
	f=fft.rfft(tpdata)[:nb]
	d=fft.irfft(f)
	tmpnum=np.argmax(fft.irfft(f0*f.conj()))
	d0=np.append(d0[tmpnum:],d0[:tmpnum])
	f0=fft.rfft(d0)
	errbinnum=np.min([int(nb/6),20])
	sr0=np.std(f0.real[-errbinnum:])
	si0=np.std(f0.imag[-errbinnum:])
	sr=np.std(f.real[-errbinnum:])
	si=np.std(f.imag[-errbinnum:])
	df=f0/f
	err=np.sqrt((sr**2+si**2)/np.abs(f)**2+(sr0**2+si0**2)/np.abs(f0)**2)
	ang=np.angle(df)
	fitnum=lnumber
	popt,pcov=so.curve_fit(lin,np.arange(1,fitnum),ang[1:fitnum],p0=[0.0],sigma=err[1:fitnum])
	dt=-popt[0]/(2*np.pi)+tmpnum/((nb-1)*2)
	dterr=pcov[0,0]**0.5/(2*np.pi)
	return [dt,dterr]
#
command=' '.join(command)
#
result=np.zeros([nsub_new.sum(),6])
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
			zchan=np.array(list(set(map(int,info['zchan'].split(','))).union(zchan)))
	elif 'zchan' in info.keys():
		zchan=np.int32(info['zchan'].split(','))
	else:
		zchan=np.int32([])
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
	if args.freq_align:
		chanstart,chanend=np.int16(np.round((np.array([freq_s,freq_e])-freq)/channel_width+0.5*nchan))
	else:
		chanstart,chanend=0,nchan
	#
	nchan_new=chanend-chanstart
	rchan=np.array(list(set(range(chanstart,chanend))-set(list(zchan))))-chanstart
	if not args.tscrunch:
		phase0=int(info['phase0'])
		sub_nperiod=int(info['sub_nperiod'])*np.ones(nperiod,dtype=np.float64)
		sub_nperiod[-1]=int(info['sub_nperiod_last'])
		middle=sub_nperiod.cumsum()-sub_nperiod/2
		time0=np.linspace(0,np.float64(info['length']),nperiod)
		psr=pm.psr_timing(pr.psr(info['psr_par']),te.times(te.time(np.float64(info['stt_date'])*np.ones(nperiod,dtype=np.float64),np.float64(info['stt_sec'])+time0)),freq_s)
		phase1=psr.phase
		phase_start=phase1.date[0]
		chebc=nc.chebfit(time0,phase1.date-phase1.date[0]+phase1.second,int(nperiod/2)+1)
		chebd=nc.chebder(chebc)
		middle_time=np.interp(middle,phase1.date-phase0+phase1.second,time0)[sub_s:sub_e]
		psr1=pm.psr_timing(pr.psr(info['psr_par']),te.times(te.time(np.float64(info['stt_date'])*np.ones(nsub_new[k],dtype=np.float64),np.float64(info['stt_sec'])+middle_time)),freq_s)
		middle_phase=psr1.phase
		data1=d.period_scrunch(sub_s,sub_e)[chanstart:chanend,:,0]
		if not args.dm_corr:
			freq_real=np.linspace(freq_start,freq_end,nchan+1)[:-1]*psr.vchange.mean()
			if 'best_dm' in info.keys():
				ddm=np.float64(info['best_dm'][0])-np.float64(info['dm'])
			else:
				ddm,ddmerr=dmcor(data1,freq_real,rchan,np.float64(info['period']),output=0)
			#print('The relative DM from the template file to data file is '+str(ddm-ddm0))
			dm_new=ddm+np.float64(info['dm'])
		else:
			dm_new=np.float64(info['dm'])
		for s in np.arange(nsub_new[k]):
			data=d.read_period(s+sub_s)[chanstart:chanend,:,0]
			if not args.dm_corr:
				fftdata=fft.rfft(data,axis=1)
				tmp=np.shape(fftdata)[-1]
				const=(1/freq_real**2*4148.808/np.float64(info['period'])*np.pi*2.0).repeat(tmp).reshape(-1,tmp)*np.arange(tmp)
				data1=shift(fftdata,const*ddm)
			tpdata=data1[rchan].mean(0)
			middle_int=middle_phase.date[s]-phase_start
			middle_offs=middle_phase.second[s]
			dp,dpe=poa(tpdata0,tpdata)
			dp0=dp+np.round(middle_offs-dp)
			roots=nc.chebroots(chebc-([dp0+middle_int]+[0]*(int(nperiod/2)+1)))
			roots=np.real(roots[np.isreal(roots)])
			root=roots[np.argmin(np.abs(roots))]
			toa=te.time(np.float64(info['stt_date']),root+np.float64(info['stt_sec']))
			period0=1/nc.chebval(root,chebd)
			toae=dpe*period0
			if output=='screen':
				sys.stdout=sys.__stdout__
				print(toa.date[0],toa.second[0],dpe,freq,dm_new,period0)
				sys.stdout=open(os.devnull,'w')
			else:
				result[cumsub[k]+s]=[toa.date[0],toa.second[0],toae,freq,dm_new,period0]
	else:
		phase0=int(info['phase0'])
		sub_nperiod=int(info['sub_nperiod'])*np.ones(nperiod,dtype=np.float64)
		sub_nperiod[-1]=int(info['sub_nperiod_last'])
		sub_nperiod_cumsum=sub_nperiod.cumsum()
		if sub_s==0: middle=sub_nperiod_cumsum[sub_e-1]/2
		else: middle=(sub_nperiod_cumsum[sub_s-1]+sub_nperiod_cumsum[sub_e-1])/2
		time0=np.linspace(0,np.float64(info['length']),12)
		phase1=pm.psr_timing(pr.psr(info['psr_par']),te.times(te.time(np.float64(info['stt_date'])*np.ones(12,dtype=np.float64),np.float64(info['stt_sec'])+time0)),freq_s).phase
		chebc=nc.chebfit(phase1.date-phase0+phase1.second,time0,7)
		chebd=nc.chebder(chebc)
		middle_time=nc.chebval(middle,chebc)
		psr1=pm.psr_timing(pr.psr(info['psr_par']),te.times(te.time(np.float64(info['stt_date']),np.float64(info['stt_sec'])+middle_time)),freq_s)
		data=d.period_scrunch()[chanstart:chanend,:,0]
		if not args.dm_corr:
			freq_real=np.linspace(freq_start,freq_end,nchan+1)[:-1]*psr1.vchange.mean()
			if 'best_dm' in info.keys():
				ddm=np.float64(info['best_dm'][0])-np.float64(info['dm'])
				fftdata=fft.rfft(data,axis=1)
				tmp=np.shape(fftdata)[-1]
				const=(1/freq_real**2*4148.808/np.float64(info['period'])*np.pi*2.0).repeat(tmp).reshape(-1,tmp)*np.arange(tmp)
				data=shift(fftdata,const*ddm)[rchan]
			else:
				data,ddm,ddmerr=dmcor(data,freq_real,rchan,np.float64(info['period']))
			dm_new=ddm+np.float64(info['dm'])
		else:
			dm_new=np.float64(info['dm'])
		data-=data.mean(1).reshape(-1,1)
		data/=data.std(1).reshape(-1,1)
		tpdata=data.mean(0)
		dp,dpe=poa(tpdata0,tpdata)
		middle_int,middle_offs=np.divmod(middle,1)
		dp0=dp+np.round(middle_offs-dp)
		root=nc.chebval(dp0+middle_int,chebc)
		toa=te.time(np.float64(info['stt_date']),root+np.float64(info['stt_sec']))
		period0=nc.chebval(root,chebd)
		toae=dpe*period0
		if output=='screen':
			sys.stdout=sys.__stdout__
			print(toa.date[0],toa.second[0],toae,freq,dm_new,period0)
			sys.stdout=open(os.devnull,'w')
		else:
			result[cumsub[k]]=[toa.date[0],toa.second[0],toae,freq,dm_new,period0]
#
if output=='ld':
	d1=ld.ld(name+'.ld')
	d1.write_shape([1,nsub_new.sum(),6,1])
	d1.write_chan(result,0)
	toainfo={'psrname':psrname,'history':command,'file_time':time.strftime('%Y-%m-%dT%H:%M:%S',time.gmtime()),'mode':'ToA'}
	d1.write_info(toainfo)
elif output=='txt':
	fout=open('a.txt','w')
	fout.write(psrname+' ToA\n')
	for i in result:
		fout.write('{:28s} {:10s} {:13f}'.format(str(int(i[0]))+str(i[1]/86400)[1:],"%.3e"%(i[2]/86400),i[5])+'\n')
	fout.close()
