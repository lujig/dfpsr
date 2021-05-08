#!/usr/bin/env python
import numpy as np
import numpy.polynomial.chebyshev as nc
import argparse as ap
import numpy.fft as fft
import os,ld,time
import scipy.optimize as so
#
version='JigLu_20201202'
parser=ap.ArgumentParser(prog='ldtoa',description='Get the ToA and DM of the ld file.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",help="input ld file")
parser.add_argument('-t',dest='template',help="template ld file")
parser.add_argument('-f',dest='fscrunch',default=0,type=np.int16,help="obtain results after frequency scrunch by this factor")
parser.add_argument('-T',action='store_true',default=False,dest='tscrunch',help='time scrunch to one subint to obtain result')
parser.add_argument('-r','--frequency_range',default=0,dest='freqrange',help='calculate in the frequency range (FREQ0,FREQ1)')
parser.add_argument('-s','--subint_range',default=0,dest='subint_range',help='calculate in the subint range (SUBINT0,SUBINT1)')
parser.add_argument("-z","--zap",dest="zap_file",default=0,help="file recording zap channels")
parser.add_argument("-a","--toa",dest="toa",action='store_true',default=False,help="save all ToA at different frequency.")
parser.add_argument("-o","--output",dest="output",default="toa",help="outputfile name")
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
nchan0=chanend-chanstart
res=nchan0
tpdata=np.zeros([nchan_new,nbin0])
i_new=0
data0=d0.period_scrunch()[:,:,0]
for i in np.arange(chanstart,chanend):
	if res>nchan_new:
		res-=nchan_new
		tpdata[i_new]+=data0[i]
	else:
		chan_data=data0[i]
		tpdata[i_new]+=chan_data*(res*1.0/nchan_new)
		i_new+=1
		if i==chanend-1: break
		tpdata[i_new]=chan_data*((nchan_new-res)*1.0/nchan_new)
		res=nchan0-(nchan_new-res)
#
freq=(freq0+freq1)/2.0
bandwidth=freq1-freq0
channel_width=bandwidth/nchan
chanstart,chanend=np.int16(np.round((np.array([freq0,freq1])-freq)/channel_width+0.5*nchan))
#
def lin(x,k):
	return k*x
#
def dmdt(df,ddm,dt0):
	return ddm/df**2*4148.808/np.float64(info['period'])+dt0
#
def poa(tpdata0,tpdata,fulloutput=False):
	import matplotlib.pyplot as plt
	plt.figure(1)
	nb=int(min(nbin0,nbin)//2+1)
	dt,dterr=np.zeros(nchan_new),np.zeros(nchan_new)
	for i in np.arange(nchan_new):
		d0=tpdata0[i]
		d=tpdata[i]
		d-=d.mean()
		d/=d.max()
		d0-=d0.mean()
		d0/=d0.max()
		f0=np.fft.rfft(d0)[:nb]
		f=np.fft.rfft(d)[:nb]
		d0=np.fft.irfft(f0)
		d=np.fft.irfft(f)
		tmpnum=np.argmax(np.correlate(d0,d,mode='full'))-(nb-1)*2+1
		d0=np.append(d0[tmpnum:],d0[:tmpnum])
		f0=np.fft.rfft(d0)
		df=f0/f
		sr0=np.std(f0.real[-10:])
		si0=np.std(f0.imag[-10:])
		sr=np.std(f.real[-10:])
		si=np.std(f.imag[-10:])
		err=np.sqrt((sr**2+si**2)/np.abs(f)**2+(sr0**2+si0**2)/np.abs(f0)**2)
		ang=np.angle(df)
		num=20#np.where(err>0.2)[0][0]
		ang=ang[1:num]
		err=err[1:num]
		popt,pcov=so.curve_fit(lin,np.arange(1,num),ang,p0=[0.0],sigma=err)
		dt[i]=-popt/(2*np.pi)+tmpnum/(nb-1)/2
		dterr[i]=pcov[0,0]**0.5/(2*np.pi)
	popt,pcov=so.curve_fit(dmdt,(freq_new[:-1]+freq_new[1:])/2,dt,p0=[0.0,dt.mean()],sigma=dterr)
	dm,t1=popt
	dmerr,t1err=np.diagonal(pcov)**0.5
	t0=(dt/dterr**2).sum()/(1/dterr**2).sum()
	t0err=1/((1/dterr**2).sum())**0.5
	if fulloutput:
		return [[dm,dmerr],[t1,t1err],[t0,t0err]],dt,dterr
	else:
		return [dm,dmerr],[t1,t1err],[t0,t0err]
#
if args.freqrange:
	freq_start,freq_end=np.float64(args.freqrange.split(','))
	chanstart,chanend=np.int16(np.round((np.array([freq_start,freq_end])-freq)/channel_width+0.5*nchan))
freq_new=np.linspace(freq_start,freq_end,nchan_new+1)
freq=np.linspace(freq_start,freq_end,nchan0+1)
d1=ld.ld(name+'.ld')
if args.toa:
	d1.write_shape([nchan_new,nsub_new,2,1])
else:
	d1.write_shape([1,nsub_new,8,1])
if nsub_new>1:
	p0=np.zeros([nsub_new,3,2])
	for s in np.arange(nsub_new):
		tpdata0=np.zeros([nchan_new,nbin])
		data=d.read_period(s+sub_start)[chanstart:chanend,:,0]
		for i in np.arange(chanend-chanstart):
			if i in zchan: continue
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
		if args.toa:
			p0[s],dt,dterr=poa(tpdata0,tpdata,fulloutput=True)
			d1.write_period(np.array([dt,dterr]).T,s)
		else:
			p0[s]=poa(tpdata0,tpdata)
			d1.write_period(np.array([np.float64(info['stt_date'][1:-1]),np.float64(info['stt_sec'][1:-1])+(s+0.5)*np.float64(info['sub_nperiod'])*np.float64(info['period']),*p0[s].reshape(-1)]),s)
		print(np.float64(info['stt_date'][1:-1]),np.float64(info['stt_sec'][1:-1])+(s+0.5)*np.float64(info['sub_nperiod'])*np.float64(info['period']),*p0[s].reshape(-1))
else:
	tpdata0=np.zeros([nchan_new,nbin])
	data=d.period_scrunch()[chanstart:chanend,:,0]
	for i in np.arange(chanend-chanstart):
		if i in zchan: continue
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
	if args.toa:
		p0,dt,dterr=poa(tpdata0,tpdata,fulloutput=True)
		d1.write_period(np.array([dt,dterr]).T,s)
	else:
		p0=poa(tpdata0,tpdata)
		d1.write_period(np.array([np.float64(info['stt_date'][1:-1]),np.float64(info['stt_sec'][1:-1])+(s+0.5)*np.float64(info['sub_nperiod'])*np.float64(info['period']),*p0[s].reshape(-1)]),s)
	print(np.float64(info['stt_date'][1:-1]),np.float64(info['stt_sec'][1:-1])+0.5*nperiod*np.float64(info['period']),*p0.reshape(-1))
#
#print(p0)
if args.toa:
	info['ToA']=p0.reshape(-1)
info['mode']='ToA'
d1.write_info(info)
#!/usr/bin/env python
