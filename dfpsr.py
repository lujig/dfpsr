#!/usr/bin/env python
import numpy as np
import numpy.fft as fft
import numpy.polynomial.chebyshev as nc
import argparse as ap
import os,time,ld
import threading
try:
	import astropy.io.fits as ps
except:
	import pyfits as ps
#
version='JigLu_20180507'
#
parser=ap.ArgumentParser(prog='dfpsr',description='Dedisperse and Fold the psrfits data.',epilog='Ver '+version)
parser.add_argument('-v','--version', action='version', version=version)
parser.add_argument("filename",nargs='+',help="name of file or filelist")
parser.add_argument("-o","--output",dest="output",default="psr",help="output file name")
parser.add_argument("-f","--frequency",dest='freqrange',default=0,help="output frequency range (MHz) in form start_freq,end_freq")
parser.add_argument('-d','--dm',dest='dm',default=0,type=np.float64,help="dispersion measure")
parser.add_argument('-p','--period',dest='period',default=0,type=np.float64,help="dispersion measure")
parser.add_argument('-n','--pulsar_name',default=0,dest='psr_name',help='input pulsar name')
parser.add_argument('-e','--pulsar_ephemeris',default=0,dest='par_file',help='input pulsar parameter file')
parser.add_argument("-c","--coefficients_num",dest="ncoeff",default=12,type=int,help="numbers of Chebyshev polynomial coefficients on time axis")
parser.add_argument("-b","--nbin",dest="nbin",default=0,type=int,help="number of phase bins in each period")
parser.add_argument("-z","--zap",dest="zap_file",default=0,help="file recording zap channels")
parser.add_argument("-r","--reverse",action="store_true",default=False,help="reverse the band")
parser.add_argument("-t","--test",action="store_true",default=False,help="generate a subint scrunched test datafile test.npy to help in judging noisy channel")
parser.add_argument("-m","--multithreading",dest='threads',default=1,type=int,help="deal with multi threading")
args=(parser.parse_args())
command=['dfpsr.py']
#
filelist=args.filename
filenum=len(filelist)
file_t0=[]
file_len=[]
#
def file_error(para):
	parser.error("Files have different parameters: "+para+".")
#
for i in np.arange(filenum):
	if not os.path.isfile(filelist[i]):
		parser.error('Fits file name is invalid.')
	f=ps.open(filelist[i],mmap=True)
	head=f['PRIMARY'].header
	subint=f['SUBINT']
	subint_header=subint.header
	subint_data=subint.data[0]
	if i==0:
		telename=head['TELESCOP']
		npol=subint_header['NPOL']
		nchan=head['OBSNCHAN']
		freq=head['OBSFREQ']
		bandwidth=subint_header['NCHAN']*subint_header['CHAN_BW']
		tsamp=subint_header['TBIN']
		nsblk=subint_header['NSBLK']
	else:
		# if telename!=head['TELESCOP']:
			# file_error('telescope name')
		# elif npol!=subint_header['NPOL']:
			# file_error('number of polorisations')
		# elif nchan!=head['OBSNCHAN']:
			# file_error('number of channels')
		if freq!=head['OBSFREQ']:
			file_error('central frequency')
		elif bandwidth!=subint_header['NCHAN']*subint_header['CHAN_BW']:
			file_error('bandwidth')
		elif tsamp!=subint_header['TBIN']:
			file_error('sampling time')
		#
	bw_sign=(bandwidth>0)
	bandwidth=np.abs(bandwidth)
	stt_imjd=head['STT_IMJD']
	stt_smjd=head['STT_SMJD']
	stt_offs=head['STT_OFFS']
	nsub=subint_header['NAXIS2']
	#
	subint_t0=(subint_data['OFFS_SUB']-tsamp*nsblk/2.0+stt_smjd+stt_offs)/86400.0+stt_imjd
	file_len.append(nsub*nsblk)
	file_t0.append(subint_t0)
	del subint_data
	f.close()
#
file_len,file_t0,filelist=np.array(file_len),np.array(file_t0),np.array(filelist)
file_len,file_t0,filelist=file_len[np.argsort(file_t0)],np.sort(file_t0),filelist[np.argsort(file_t0)]
if len(file_len)>1:
	if np.max(np.abs((file_len*tsamp/86400.0+file_t0)[:-1]-file_t0[1:]))>(tsamp/86400.0):
		parser.error("Files are not continuous.")
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
#
nbin=file_len.sum()
stt_time=file_t0[0]
freq_start,freq_end=(np.array([chanstart,chanend])-0.5*nchan)*channel_width+freq
info={'nbin_origin':nbin,'telename':telename,'freq_start':freq_start,'freq_end':freq_end,'nchan':chanend-chanstart,'tsamp_origin':tsamp,'stt_time_origin':stt_time,'stt_time':stt_time,'npol':npol}
#
if args.psr_name and args.par_file:
	parser.error('At most one of flags -n and -p is required.')
elif args.psr_name or args.par_file:
	if args.period:
		parser.error('With pulsar name or ephemeris, period is needless.')
	elif args.psr_name:
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
	if not (args.dm or args.period):
		parser.error('Pulsar Parameter should be provided.')
	if not (args.dm and args.period):
		parser.error('Both DM and period should be provided.')
	period=args.period
	dm=args.dm
	command.append('-p '+str(args.period)+' -d '+str(args.dm))
#
command.append('-c '+str(args.ncoeff))
if args.nbin:
	command.append('-b '+str(args.nbin))
	if args.nbin>(period/tsamp):
		parser.error('Provided phase bin number in each period is too large.')
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
d=ld.ld(name+'.ld')
info['file_time']=time.strftime('%Y-%m-%dT%H:%M:%S',time.gmtime())
info['dm']=dm
nchan_new=chanend-chanstart
d.write_shape([nchan_new,nbin/nsblk,nsblk,1])
#
if args.test:
	command.append('-t')
	testdata=np.zeros([nchan_new,nsblk])
#
cumsub=0
for n in np.arange(filenum):
	f=ps.open(filelist[n],mmap=True)
	fsub=f['SUBINT'].header['naxis2']
	for i in np.arange(fsub):
		data=(np.int16(f['SUBINT'].data[i]['DATA'].reshape(nsblk,npol,nchan).sum(1)).transpose())
		if args.reverse or (not bw_sign):
			if nchan==chanend:
				data=data[(nchan-chanstart-1)::-1,:]
			else:
				data=data[(nchan-chanstart-1):(nchan-chanend-1):-1,:]
			command.append('-r')
		else:
			data=data[chanstart:chanend,:]
		d.write_period(data,cumsub)
		if args.test:
			testdata+=data
		cumsub+=1
		del f['SUBINT'].data
	f.close()
if args.test:
	test=ld.ld('test.ld')
	test.write_shape([nchan_new,1,nsblk,1])
	test.write_period(testdata,0)
	testinfo={'nbin':nbin,'test':True,'freq_start':freq_start,'freq_end':freq_end,'nchan':chanend-chanstart}
	test.write_info(testinfo)
#
command=' '.join(command)
info['history']=command
#
nbin01=np.int32(2**np.ceil(np.log2(nbin)))
nbin_old=nbin
freq0=freq_start+channel_width/2.0
freq1=freq_end-channel_width/2.0
nbin0=nbin-np.int32(np.ceil((1/freq0**2-1/freq1**2)*dm*4148.808/tsamp)+1)
info['stt_time']-=1/freq1**2*dm*4148.808/tsamp/86400.0
disp_time=1/np.linspace(freq0,freq1,nchan_new)**2*dm*4148.808
disp=disp_time/tsamp*np.pi*2.0/nbin01
disp=disp-np.min(disp)
def shift(y,x):
	fftp=fft.rfft(y)
	ffts=fftp*np.exp(x*1j*np.arange(len(fftp)))
	fftr=fft.irfft(ffts)
	return fftr
#
stt_time=info['stt_time']
end_time=stt_time+tsamp*nbin0/86400.0+60./86400
stt_time-=60./86400
stt_time=str(int(stt_time))+str(stt_time%1)[1:]
end_time=str(int(end_time))+str(end_time%1)[1:]
if args.period or (not pepoch):
	if args.period:
		period=args.period
	phase=np.arange(nbin0)*tsamp/period
	info['phase0']=0
	nperiod=int(np.ceil(np.max(phase)))
else:
	os.popen('tempo2 -f '+par_file+' -pred \"'+telename+' '+stt_time+' '+end_time+' '+str(freq_start)+' '+str(freq_end)+' '+str(args.ncoeff)+' 2 '+str(int(tsamp*nbin0)+150)+'\"')
	predictor_file='t2pred.dat'
	#
	polyco=open(predictor_file,'r')
	lines=polyco.readlines()
	polyco.close()
	os.remove(predictor_file)
	os.remove('pred.tim')
	if not args.par_file:
		os.remove(par_file)
	coeff=[]
	predictor=[]
	for line in lines:
		predictor.append(line)
		elements=line.split()
		if elements[0]=='TIME_RANGE':
			t0=np.float64(elements[1])
			t1=np.float64(elements[2])
		elif elements[0]=='FREQ_RANGE':
			f0=np.float64(elements[1])
			f1=np.float64(elements[2])
		elif elements[0]=='DISPERSION_CONSTANT':
			dispc=np.float64(elements[1])
		elif elements[0]=='COEFFS':
			coeff.append(np.float64(elements[1:]))
		elif line=="ChebyModel END\n" or line=="ChebyModel END":
			break
	info['predictor']=predictor[1:]
	coeff=np.array(coeff)
	coeff[0,:]/=2.0
	coeff[:,0]/=2.0
	dt=(((np.arange(nbin0)+0.5)*tsamp-np.min(disp_time))/86400+info['stt_time']-t0)/(t1-t0)*2-1
	phase=nc.chebval(dt,coeff[:,0])
	phase0=np.ceil(phase[0])
	nperiod=int(np.floor(phase[-1]-phase0))
	period=(nbin0-1)*tsamp/(phase[-1]-phase[0])
	coeff[0,0]-=phase0
	roots=nc.chebroots(coeff[:,0])
	roots=np.real(roots[np.isreal(roots)])
	info['stt_time']=(roots[np.argmin(np.abs(roots-dt[0]))]+1)/2.0*(t1-t0)+t0
	info['phase0']=phase0
	phase-=phase0
#
info['nperiod']=nperiod
info['period']=period
#
nbin_max=(nbin0-1)/(np.max(phase)-np.min(phase))
if args.nbin:
	nbin=args.nbin
	temp_multi=2**(np.int16(np.log2(nbin_max/nbin))+1)
else:
	nbin=2**np.int16(np.log2(nbin_max))
	temp_multi=2
info['nbin']=nbin
info['length']=period*nperiod
#
phase*=(nbin*temp_multi)
phase-=0.5
phasenum=np.int32(phase)
totalbin=nperiod*nbin*temp_multi
judge0=len(phasenum)-np.sum(phasenum>=0)
judge1=np.sum(phasenum<totalbin)
phasenum=phasenum[judge0:judge1]
phase=phase[judge0:judge1]
phaseres=(phase-phasenum)
binint=np.zeros(totalbin+1)
binint[phasenum]=1-phaseres
binint[phasenum+1]+=phaseres
binint=binint[:-1].reshape(-1,temp_multi).sum(1)
binint[binint==0]=1
#
d.write_shape([nchan_new,nperiod,nbin,1])
t=time.time()
if args.threads<=0:
	parser.error('Threads number should be positive.')
elif args.threads==1:
	for k in np.arange(nchan_new):
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
else:
	np.save(name+'_df_assist_temp',{'nchan_new':nchan_new,'zchan':zchan,'totalbin':totalbin,'temp_multi':temp_multi,'nbin0':nbin0,'nbin01':nbin01,'nbin_old':nbin_old,'disp':disp,'judge0':judge0,'judge1':judge1,'phasenum':phasenum,'phaseres':phaseres,'binint':binint})
	if nchan_new<args.threads:
		threads=nchan_new
		start_list=np.arange(0,nchan_new)
		len_list=np.ones(threads)
	elif nchan_new%args.threads==0:
		threads=args.threads
		space=nchan_new/threads
		start_list=np.arange(0,nchan_new,space)
		len_list=np.int32(np.ones(threads))*space
	else:
		space=nchan_new/args.threads+1
		threads=nchan_new/space+1
		start_list=np.arange(0,nchan_new,space)
		len_list=np.int32(np.ones(threads))*space
		len_list[-1]-=space*threads-nchan_new
	def df(i):
		if os.path.isfile('df_assist.py'):
			os.system('python df_assist.py -s '+str(start_list[i])+' -n '+str(len_list[i])+' '+name)
		else:
			os.system('df_assist.py -s '+str(start_list[i])+' -n '+str(len_list[i])+' '+name)
	thlist=[]
	testlist=[]
	for i in np.arange(threads):
		th=threading.Thread(target=df, args=(i,))
		thlist.append(th)
		th.start()
	for i in thlist:
		i.join()
	os.remove(name+'_df_assist_temp.npy')
d.write_info(info)
#

