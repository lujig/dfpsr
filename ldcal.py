#!/usr/bin/env python
import numpy as np
import numpy.fft as fft
import numpy.polynomial.chebyshev as nc
import argparse as ap
import os,sys,time,ld
try:
	import astropy.io.fits as ps
except:
	import pyfits as ps
#
version='JigLu_20200923'
#
parser=ap.ArgumentParser(prog='ldcal',description='Dedisperse and Fold the psrfits data.',epilog='Ver '+version)
parser.add_argument('-v','--version', action='version', version=version)
parser.add_argument('--verbose', action="store_true",default=False,help="print detailed information")
parser.add_argument("filename",nargs='+',help="name of file or filelist")
parser.add_argument("--cal_period",dest="cal_period",default=0,type=np.float64,help="period of the calibration fits file (s)")
parser.add_argument("-o","--output",dest="output",default="psr",help="output file name")
parser.add_argument("-r","--reverse",action="store_true",default=False,help="reverse the band")
args=(parser.parse_args())
command=['ldcal.py']
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
	file_len.append(nsub*nsblk)
	file_t0.append(subint_t0)
#
file_len,file_t0,filelist=np.array(file_len),np.array(file_t0),np.array(filelist)
sorts=np.argsort(file_t0)
file_len,file_t0,filelist=file_len[sorts],np.sort(file_t0),filelist[sorts]
if len(file_len)>1:
	if np.max(np.abs((file_len*tsamp/86400.0+file_t0)[:-1]-file_t0[1:]))>(tsamp/86400.0):
		parser.error("Data files are not continuous.")
#
channel_width=bandwidth*1.0/nchan
#
nbin=file_len.sum()
stt_time=file_t0[0]+tsamp/2.0
freq_start,freq_end=np.array([-0.5,0.5])*nchan*channel_width+freq
info={'nbin_origin':nbin,'telename':telename,'freq_start':freq_start,'freq_end':freq_end,'nchan':nchan,'tsamp_origin':tsamp,'stt_time':stt_time,'npol':npol,'mode':'cal','length':nbin*tsamp}
#
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
command=' '.join(command)
info['history']=command
#
sys.stdout.write('Processing the noise file...\n')
noisen=np.int64(args.cal_period//tsamp)
noise_data=np.zeros([noisen,npol,nchan])
noise_cum=np.zeros(noisen)
cumsub=0
for n in np.arange(filenum):
	f=ps.open(filelist[n],mmap=True)
	fsub=f['SUBINT'].header['naxis2']
	for i in np.arange(fsub):
		dtmp=f['SUBINT'].data[i]
		data=np.int16(dtmp['DATA'].reshape(nsblk,npol,nchan)*dtmp['dat_scl'].reshape(1,npol,nchan)+dtmp['dat_offs'].reshape(1,npol,nchan))
		del f['SUBINT'].data
		if args.reverse or (not bw_sign):
			data=data[:,:,::-1]
		noise_t=np.int64((np.arange(nsblk)+cumsub*nsblk)*tsamp%args.cal_period//tsamp)
		for k in np.arange(nsblk):
			tmp_noise_t=noise_t[k]
			if tmp_noise_t==noisen:
				continue
			noise_data[tmp_noise_t]+=data[k]
			noise_cum[tmp_noise_t]+=1
		cumsub+=1
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
#
if args.verbose:
	sys.stdout.write('Constructing the output file...\n')
#
d=ld.ld(name+'.ld')
d.write_shape([nchan,1,1,4])
d.write_period(np.array([noise_a12,noise_a22,noise_cos,noise_sin]).T,0)
d.write_info(info)
