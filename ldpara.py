#!/usr/bin/env python
import numpy as np
import argparse as ap
import os,time,ld,sys
#
version='JigLu_20200930'
parser=ap.ArgumentParser(prog='ldpara',description='Show the parameters of ld file.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",help="input ld file")
parser.add_argument('-c',dest='paras',default=0,help="parameter name list, include nsub, nchan, nbin, npol, stt_time, file_time, psr_name, period, nperiod, dm, freq, bw and length")
#
args=(parser.parse_args())
if not os.path.isfile(args.filename):
	parser.error('A valid ld file name is required.')
d=ld.ld(args.filename)
info=d.read_info()
#
if not args.paras:
	quit()
plist=args.paras.split(',')
for pname in plist:
	if pname=='nsub':
		if 'compressed' in info.keys():
			sys.stdout.write(pname+' '+info['nsub_new']+'\n')
		else:
			sys.stdout.write(pname+' '+info['nsub']+'\n')
	elif pname=='nchan':
		if 'compressed' in info.keys():
			sys.stdout.write(pname+' '+info['nchan_new']+'\n')
		else:
			sys.stdout.write(pname+' '+info['nchan']+'\n')
	elif pname=='nbin':
		if 'compressed' in info.keys():
			sys.stdout.write(pname+' '+info['nbin_new']+'\n')
		else:
			sys.stdout.write(pname+' '+info['nbin']+'\n')
	elif pname=='nbin':
		if 'compressed' in info.keys():
			sys.stdout.write(pname+' '+info['npol_new']+'\n')
		else:
			sys.stdout.write(pname+' '+info['npol']+'\n')
	elif pname=='shape':
		sys.stdout.write(pname+' '+str(tuple(d.read_shape()))+'\n')
	elif pname in ['stt_time', 'file_time', 'psr_name', 'nperiod', 'period', 'dm', 'length', 'mode']:
		sys.stdout.write(pname+' '+info[pname]+'\n')
	elif pname=='freq':
		sys.stdout.write(pname+' '+str((np.float64(info['freq_end'])+np.float64(info['freq_start']))/2)+'\n')
	elif pname=='bw':
		sys.stdout.write(pname+' '+str(np.float64(info['freq_end'])-np.float64(info['freq_start']))+'\n')
	else:
		sys.stdout.write('Parameter '+pname+' can not be found.'+'\n')
