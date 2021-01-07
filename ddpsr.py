#!/usr/bin/env python
import numpy as np
import numpy.fft as fft
import numpy.polynomial.chebyshev as nc
import argparse as ap
import os,sys,time,ld
import mpmath as mm
mm.mp.dps=30
import astropy.io.fits as ps
from multiprocessing import Pool
import gc
#
version='JigLu_20201125'
#
parser=ap.ArgumentParser(prog='dfpsr',description='Dedisperse and Fold the psrfits data.',epilog='Ver '+version)
parser.add_argument('-v','--version', action='version', version=version)
parser.add_argument('--verbose', action="store_true",default=False,help="print detailed information")
parser.add_argument("filename",nargs='+',help="name of file or filelist")
parser.add_argument("-a","--cal",dest="cal",nargs='+',help="name of calibration file or calibration filelist")
parser.add_argument("--cal_period",dest="cal_period",default=0,type=np.float64,help="period of the calibration fits file (s)")
parser.add_argument("-o","--output",dest="output",default="psr",help="output file name")
parser.add_argument("-f","--frequency",dest='freqrange',default=0,help="output frequency range (MHz) in form start_freq,end_freq")
parser.add_argument('-d','--dm',dest='dm',default=0,type=np.float64,help="dispersion measure")
parser.add_argument('-n','--pulsar_name',default=0,dest='psr_name',help='input pulsar name')
parser.add_argument('-e','--pulsar_ephemeris',default=0,dest='par_file',help='input pulsar parameter file')
parser.add_argument("-z","--zap",dest="zap_file",default=0,help="file recording zap channels")
parser.add_argument("-s","--fscrunch",type=int,default=0,help="frequency scrunch by this factor, if not set, the file will be scrunched to one channel")
parser.add_argument("-r","--reverse",action="store_true",default=False,help="reverse the band")
parser.add_argument("-l","--large_mem",action="store_true",default=False,help="large RAM")
parser.add_argument("-m","--multi",dest="multi",default=0,type=int,help="number of processes")
args=(parser.parse_args())
command=['dfpsr.py']
#
if args.verbose:
	sys.stdout.write('Analyzing the arguments...\n')
filelist=args.filename
filenum=len(filelist)
file_t0=[]
file_time=[]
file_len=[]
#
def file_error(para,filetype):
	parser.error("Fits "+filetype+" have different parameters: "+para+".")
#
telename,pol_type,npol,nchan,freq,bandwidth,tsamp,nsblk,bw_sign,stt_imjd,stt_smjd,stt_offs,nsub,offs_sub='','',0,0,0,0.0,0.0,0,True,0,0,0.0,0,0.0
def file_check(fname,notfirst=True,filetype='data'):
	if not os.path.isfile(fname):
		parser.error('Fits '+filetype+' name is invalid.')
	try:
		f=ps.open(filelist[i],mmap=True)
	except:
		parser.error('Fits '+filetype+' is invalid.')
	head=f['PRIMARY'].header
	subint=f['SUBINT']
	subint_header=subint.header
	subint_data=subint.data[0]
	global telename,pol_type,npol,nchan,freq,bandwidth,tsamp,nsblk,bw_sign,stt_imjd,stt_smjd,stt_offs,nsub
	if not notfirst:
		telename=head['TELESCOP']
		npol=subint_header['NPOL']
		nchan=head['OBSNCHAN']
		freq=head['OBSFREQ']
		bandwidth=subint_header['NCHAN']*subint_header['CHAN_BW']
		bw_sign=(bandwidth>0)
		bandwidth=np.abs(bandwidth)
		tsamp=subint_header['TBIN']
		nsblk=subint_header['NSBLK']
		pol_type=subint_header['POL_TYPE']
	else:
		if telename!=head['TELESCOP']:
			file_error('telescope name',filetype)
		if pol_type!=subint_header['POL_TYPE']:
			file_error('polarisation type',filetype)
		if npol!=subint_header['NPOL']:
			file_error('number of polorisations',filetype)
		if nchan!=head['OBSNCHAN']:
			file_error('number of channels',filetype)
		if freq!=head['OBSFREQ']:
			file_error('central frequency',filetype)
		if bandwidth!=np.abs(subint_header['NCHAN']*subint_header['CHAN_BW']):
			file_error('bandwidth',filetype)
		if tsamp!=subint_header['TBIN']:
			file_error('sampling time',filetype)
		#
	stt_imjd=head['STT_IMJD']
	stt_smjd=head['STT_SMJD']
	stt_offs=head['STT_OFFS']
	nsub=subint_header['NAXIS2']
	offs_sub=subint_data['OFFS_SUB']
	del subint_data
	f.close()
#
for i in np.arange(filenum):
	file_check(filelist[i],notfirst=i)
	#
	subint_t0=(offs_sub-tsamp*nsblk/2.0+stt_smjd+stt_offs)/86400.0+stt_imjd
	file_time.append([offs_sub-tsamp*nsblk/2.0,stt_smjd,stt_offs,stt_imjd])
	file_len.append(nsub)
	file_t0.append(subint_t0)
#
file_len,file_t0,filelist,file_time=np.array(file_len),np.array(file_t0),np.array(filelist),np.array(file_time)
sorts=np.argsort(file_t0)
file_len,file_t0,filelist,file_time=file_len[sorts],np.sort(file_t0),filelist[sorts],file_time[sorts]
if len(file_len)>1:
	if np.max(np.abs((file_len*nsblk*tsamp/86400.0+file_t0)[:-1]-file_t0[1:]))>(tsamp/86400.0):
		parser.error("Data files are not continuous.")
#
if args.cal:
	command.append('-a ')
	if len(args.cal)==1:
		if args.cal[0][-3:]=='.ld':
			noise_mark='ld'
			noise=ld.ld(args.cal[0])
			noise_info=noise.read_info()
			if noise_info['mode']!='cal':
				parser.error("LD file is not caliration file.")
			elif telename!=noise_info['telename']:
				parser.error("LD calibration file has different telescope name.")
			elif nchan!=int(noise_info['nchan']):
				parser.error("LD calibration file has different channel number.")
		else:
			noise_mark='fits'
	else:
		noise_mark='fits'
	if noise_mark=='fits':
		if not args.cal_period:
			parser.error("Noise period is not given.")
		noiselist=args.cal
		noisenum=len(noiselist)
		noise_t0,noise_len=[],[]
		for i in np.arange(noisenum):
			file_check(noiselist[i],filetype='noise')
			subint_t0=(offs_sub-tsamp*nsblk/2.0+stt_smjd+stt_offs)/86400.0+stt_imjd
			noise_len.append(nsub)
			noise_t0.append(subint_t0)
		#
		noise_t0,noise_len,noiselist=np.array(noise_t0),np.array(noise_len),np.array(noiselist)
		sorts=np.argsort(noise_t0)
		noise_t0,noise_len,noiselist=noise_t0[sorts],noise_len[sorts],noiselist[sorts]
		if len(noise_len)>1:
			if np.max(np.abs((noise_len*nsblk*tsamp/86400.0+noise_t0)[:-1]-noise_t0[1:]))>(tsamp/86400.0):
				parser.error("Noise files are not continuous.")
else:
	noise_mark=''
#
channel_width=bandwidth*1.0/nchan
if args.freqrange:
	command.append('-f '+args.freqrange)
	freq_start,freq_end=np.float64(args.freqrange.split(','))
	chanstart,chanend=np.int16(np.round((np.array([freq_start,freq_end])-freq)/channel_width+0.5*nchan))
	if freq_start>freq_end:
		parser.error("Starting frequency larger than ending frequency.")
	elif freq_start<(freq-bandwidth/2.0) or freq_end>(freq+bandwidth/2.0):
		parser.error("Input frequency is overrange.")
else:
	chanstart,chanend=0,nchan
nchan_new=chanend-chanstart
#
nbin=file_len.sum()*nsblk
stt_time=file_t0[0]
freq_start,freq_end=(np.array([chanstart,chanend])-0.5*nchan)*channel_width+freq
info={'nbin_origin':nbin,'telename':telename,'freq_start':freq_start,'freq_end':freq_end,'nchan':chanend-chanstart,'tsamp_origin':tsamp,'stt_time_origin':stt_time,'stt_time':stt_time,'npol':npol}
#
if args.psr_name and args.par_file:
	parser.error('At most one of flags -n and -p is required.')
elif args.psr_name or args.par_file:
	if args.psr_name:
		command.append('-n '+args.psr_name)
		psr_name=args.psr_name
		os.system('psrcat -e '+psr_name+' > psr.par')
		par_file=open('psr.par','r')
		psr_par=par_file.readlines()
		par_file.close()
		if len(psr_par)<3:
			parser.error('A valid pulsar name is required.')
		par_file='psr.par'
	else:
		command.append('-e')
		par_file=open(args.par_file,'r')
		psr_par=par_file.readlines()
		par_file.close()
		par_file=args.par_file
	info['psr_par']=psr_par
	pepoch=False
	for line in psr_par:
		elements=line.split()
		if elements[0]=='PSRJ':
			psr_name=elements[1].strip('\n')
			info['psr_name']=psr_name
		elif elements[0]=='F0':
			period=1./np.float64(elements[1])
		elif elements[0]=='DM':
			if not args.dm:
				dm=np.float64(elements[1])
			else:
				dm=args.dm
		elif elements[0]=='PEPOCH':
			pepoch=True
else:
	if not args.dm:
		parser.error('Pulsar DM should be provided.')
	dm=args.dm
	command.append(' -d '+str(args.dm))
#
info['mode']='dedisperse'
#
info['dm']=dm
#
if args.fscrunch:
	fscrunch=args.fscrunch
	if fscrunch<0: 
		parser.error('The scrunch factor should be positive.')
	elif fscrunch>nchan_new:
		parser.error('The scrunch factor should be smaller than the channel number.')
	elif nchan_new%fscrunch!=0:
		parser.error('The channel number should be divisible by the scrunch factor.')
	nchan1=nchan_new/fscrunch
	command.append('-s '+str(fscrunch))
#
if args.zap_file:
	command.append('-z')
	if not os.path.isfile(args.zap_file):
		parser.error('The zap channel file is invalid.')
	zchan=np.loadtxt(args.zap_file,dtype=np.int32)
	if np.max(zchan)>=nchan or np.min(zchan)<0:
		parser.error('The zapped channel number is overrange.')
	info['zchan']=str(list(zchan))[1:-1]
else:
	zchan=[]
name=args.output
if os.path.isfile(name):
	parser.error('The name of output file already existed. Please provide a new name.')
if len(name)>3:
	if name[-3:]=='.ld':
		name=name[:-3]
#
if args.reverse:
	command.append('-r')
#
if args.multi:
	if args.multi>20:
		parser.error('The processes number is too large!')
	command.append('-m '+str(args.multi))
#
command=' '.join(command)
info['history']=command
#
if noise_mark=='fits':
	if args.verbose:
		sys.stdout.write('Processing the noise file...\n')
	noisen=np.int64(args.cal_period//tsamp)
	noise_data=np.zeros([noisen,npol,nchan_new],dtype=np.float64)
	noise_cum=np.zeros(noisen)
	for n in np.arange(noisenum):
		f=ps.open(noiselist[n],mmap=True)
		fsub=f['SUBINT'].header['naxis2']
		for i in np.arange(fsub):
			dtmp=f['SUBINT'].data[i]
			data=np.int16(dtmp['DATA'].reshape(nsblk,npol,nchan)*dtmp['dat_scl'].reshape(1,npol,nchan)+dtmp['dat_offs'].reshape(1,npol,nchan))
			del f['SUBINT'].data
			if args.reverse or (not bw_sign):
				if nchan==chanend:
					data=data[:,:,(nchan-chanstart-1)::-1]
				else:
					data=data[:,:,(nchan-chanstart-1):(nchan-chanend-1):-1]
			else:
				data=data[:,:,chanstart:chanend]
			cumsub=noise_len[:n].sum()+i
			noise_t=np.int64(np.round((np.arange(nsblk)+cumsub*nsblk)*tsamp%args.cal_period/tsamp))
			for k in np.arange(nsblk):
				tmp_noise_t=noise_t[k]
				if tmp_noise_t>=noisen:
					continue
				noise_data[tmp_noise_t]+=data[k]
				noise_cum[tmp_noise_t]+=1
		f.close()
	tmp_noise=noise_data[:,0].sum(1)/noise_cum
	sorts=np.argsort(tmp_noise)
	noise_data,noise_cum=noise_data[sorts],noise_cum[sorts]
	noisen_center=np.int64(noisen//2)
	noise_off=noise_data[3:(noisen_center-1)].sum(0)/noise_cum[3:(noisen_center-1)].sum().reshape(-1,1)
	noise_on=noise_data[(noisen_center+2):-3].sum(0)/noise_cum[(noisen_center+2):-3].sum().reshape(-1,1)-noise_off
	noise_a12,noise_a22=noise_on[:2]
	noise_dphi=np.arctan2(noise_on[3],noise_on[2])
	noise_cos,noise_sin=np.cos(noise_dphi),np.sin(noise_dphi)
elif noise_mark=='ld':
	noise_a12,noise_a22,noise_cos,noise_sin=noise.read_period(0)[chanstart:chanend,0].T
if args.cal:
	info['cal']=list(map(str,np.array([noise_a12,noise_a22,noise_cos,noise_sin]).T.reshape(-1)))
	noise_a12=np.where(noise_a12>0,1./noise_a12,0)
	noise_a22=np.where(noise_a22>0,1./noise_a22,0)
	noise_a1a2=noise_a12*noise_a22
	noise_a1a2=np.sqrt(noise_a1a2)
	noise_cos=noise_cos*noise_a1a2
	noise_sin=noise_sin*noise_a1a2
	noise_a1p2=(noise_a12+noise_a22)/2.0
	noise_a1m2=(noise_a12-noise_a22)/2.0
#
if args.verbose:
	sys.stdout.write('Constructing the output file...\n')
#
freq0=freq_start
freq1=freq_end
#
roach_delay=8.192e-6*3
gps_delay=1e-5
transline_delay=2e-5
light_delay=(300.0+141.96)/3.0e8
delay=transline_delay+light_delay-gps_delay-roach_delay
#
time0=file_time[0]
stt_sec=time0[:-1].sum()-delay
stt_date=time0[-1]+stt_sec//86400
stt_sec=stt_sec%86400
#
info['stt_sec']=stt_sec
info['stt_date']=stt_date
info['stt_time']=stt_date+stt_sec/86400.0
#
info['length']=nbin*tsamp
#
d=ld.ld(name+'.ld')
info['file_time']=time.strftime('%Y-%m-%dT%H:%M:%S',time.gmtime())
d.write_shape([nchan_new,1,nbin,npol])
#
def write_data(ldfile,data,startbin,channum):
	if args.cal:
		if pol_type=='AABBCRCI':
			a12,a22,ncos,nsin=noise_a12[channum],noise_a22[channum],noise_cos[channum],noise_sin[channum]
			aa,bb,cr,ci=data
			aa0,bb0,cr0,ci0=a12*aa,a22*bb,ncos*cr+nsin*ci,-nsin*cr+ncos*ci
			i,q,u,v=aa0+bb0,aa0-bb0,cr0*2,ci0*2
			data=np.array([i,q,u,v])
		elif pol_type=='IQUV':
			a1p2,a1m2,ncos,nsin=noise_a1p2[channum],noise_a1m2[channum],noise_cos[channum],noise_sin[channum]
			i,q,u,v=data
			i,q,u,v=a1p2*i-a1m2*q,a1p2*q-a1m2*i,ncos*u+nsin*v,-nsin*u+ncos*v
			data=np.array([i,q,u,v])
		info['pol_type']='IQUV'
	else:
		info['pol_type']=pol_type
	d.__write_chanbins_add__(data.T,startbin,channum)
#
freq=np.arange(chanstart,chanstart+nchan_new)*channel_width+freq_start
if args.dm:
	dispb=4148.808*dm/freq**2/tsamp
	dispb-=4148.808*dm/freq_end**2/tsamp
else:
	stt_time=info['stt_time_origin']
	end_time=stt_time+tsamp*nbin/86400.0+60./86400
	stt_time-=60./86400
	stt_time=str(int(stt_time))+str(stt_time%1)[1:]
	end_time=str(int(end_time))+str(end_time%1)[1:]
	os.popen('tempo2 -f '+par_file+' -pred \"'+telename+' '+stt_time+' '+end_time+' '+str(freq_start)+' '+str(freq_end)+' 12 2 '+str(int(tsamp*nbin)+150)+'\"').close()
	predictor_file='t2pred.dat'
	#
	polyco=open(predictor_file,'r')
	lines=polyco.readlines()
	polyco.close()
	os.remove(predictor_file)
	os.remove('pred.tim')
	coeff=[]
	predictor=[]
	for line in lines:
		predictor.append(line)
		elements=line.split()
		if elements[0]=='TIME_RANGE':
			t0=mm.mpf(elements[1])
			t1=mm.mpf(elements[2])
		elif elements[0]=='DISPERSION_CONSTANT':
			dispc=np.float64(elements[1])
		elif elements[0]=='COEFFS':
			coeff.append(list(map(mm.mpf,elements[1:])))
		elif line=="ChebyModel END\n" or line=="ChebyModel END":
			break
	info['predictor']=predictor[1:]
	coeff=np.array(coeff)
	coeff[0,:]/=2.0
	coeff[:,0]/=2.0
	tmp=int(coeff[0,0])
	coeff[0,0]-=tmp
	coeff=np.float64(coeff)
	time0=file_time[0]
	t0=np.float64(t0-time0[-1])*86400.0
	t1=np.float64(t1-time0[-1])*86400.0
	dt=(np.array([0,nbin])*tsamp+time0[:-1].sum()-delay-t0)/(t1-t0)*2-1
	coeff1=coeff.sum(1)
	phase=nc.chebval(dt,coeff1)
	period=tsamp*nbin/(phase[1]-phase[0])
	dispb=-dispc*period/freq**2/tsamp
	dispb-=(-dispc*period)/freq_end**2/tsamp
dispb_i=np.int64(np.floor(dispb))
dispb_f=dispb-dispb_i
dispb_r=1-dispb_f
#
def gendata(cums,nsub,data,last=False):
	data=data.transpose(2,1,0)
	if args.reverse or (not bw_sign):
		if nchan==chanend:
			data=data[(nchan-chanstart-1)::-1]
		else:
			data=data[(nchan-chanstart-1):(nchan-chanend-1):-1]
	else:
		data=data[chanstart:chanend]
	for f in np.arange(nchan_new):
		if f+chanstart in zchan: continue
		data0=np.zeros([npol,nsub+1])
		data0[:,:nsub]=dispb_r[f]*data[f]
		data0[:,1:]+=dispb_f[f]*data[f]
		newphase=cums*nsblk+np.arange(nsub+1)-dispb_i[f]
		if newphase[-1]<0 or newphase[0]>=nbin:
			continue
		tpdata=data0[:,(newphase>=0) & (newphase<nbin)]
		newphase=newphase[(newphase>=0) & (newphase<nbin)]
		startphase=newphase[0]
		write_data(d,tpdata,startphase,f)
#
if args.verbose:
	sys.stdout.write('Dedispersing and folding the data...\n')
#
def dealdata(filelist,n):
	if args.verbose:
		sys.stdout.write('Processing the '+str(n+1)+'th fits file...\n')
		timemark=time.time()
	f=ps.open(filelist[n],mmap=True)
	fsub=f['SUBINT'].header['naxis2']
	if args.large_mem:
		cumsub=np.int64(file_len[:n].sum())
		dtmp=f['SUBINT'].data
		data=np.float64((dtmp['DATA']).reshape(fsub,nsblk,npol,nchan)*dtmp['dat_scl'].reshape(fsub,1,npol,nchan)+dtmp['dat_offs'].reshape(fsub,1,npol,nchan)).reshape(fsub*nsblk,npol,nchan)
		del f['SUBINT'].data
		f.close()
		gendata(cumsub,file_len[n]*nsblk,data,tpsub,tpsubn,last=True)
	else:
		for i in np.arange(fsub):
			cumsub=np.int64(file_len[:n].sum()+i)
			dtmp=f['SUBINT'].data[i]
			data=np.int16(dtmp['DATA'].reshape(nsblk,npol,nchan)*dtmp['dat_scl'].reshape(1,npol,nchan)+dtmp['dat_offs'].reshape(1,npol,nchan))
			del f['SUBINT'].data
			gendata(cumsub,nsblk,data,last=(i==(fsub-1)))
		f.close()
	gc.collect()
	if args.verbose:
		sys.stdout.write('Processing the '+str(n+1)+'th fits file takes '+str(time.time()-timemark)+' second.\n')
#
if args.multi:
	pool=Pool(processes=args.multi)
#
for n in np.arange(filenum):
	if args.multi:
		pool.apply_async(dealdata,(filelist,n))
	else:
		dealdata(filelist,n)
if args.multi:
	pool.close()
	pool.join()
#
if fscrunch:
	for i in np.arange(nchan1):
		fdata=d.chan_scrunch(np.int32(np.arange(fscrunch)+i*fscrunch))
		d.write_chan(fdata,np.int32(i))
	d.__size__[1:]=np.int32([nchan1,1,nbin,npol])
	d.__size__[0]=nchan1*nbin*npol*8
	d.__write_size__(d.__size__)
	d.file=open(d.name,'rb+')
	d.file.seek(24+d.__size__[0],0)
	d.file.truncate()
	d.file.flush()
	d.file.close()
	info['nchan']=np.int32(nchan1)
else:
	fdata=d.chan_scrunch()
	d.write_shape([1,1,nbin,npol])
	d.write_chan(fdata,0)
	info['nchan']=1

#
d.write_info(info)
#

