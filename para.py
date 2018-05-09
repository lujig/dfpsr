#!/usr/bin/env python
import numpy as np
import argparse as ap
import os,time,ld
#
version='JigLu_20180507'
parser=ap.ArgumentParser(prog='para',description='Show the parameters of ld file.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",help="input ld file")
parser.add_argument('-c',dest='paras',default=0,help="parameter name list, include nsub, nchan, nbin, stt_time, file_time, psr_name, period, dm, freq, bw and length")
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
			print pname,info['nsub_new']
		else:
			print pname,info['nperiod']
	elif pname=='nchan':
		if 'compressed' in info.keys():
			print pname,info['nchan_new']
		else:
			print pname,info['nchan']
	elif pname=='nbin':
		if 'compressed' in info.keys():
			print pname,info['nbin_new']
		else:
			print pname,info['nbin']
	elif pname in ['stt_time', 'file_time', 'psr_name', 'period', 'dm', 'length']:
		print pname,info[pname]
	elif pname=='freq':
		print pname,str((np.float64(info['freq_end'])+np.float64(info['freq_start']))/2)
	elif pname=='bw':
		print pname,str(np.float64(info['freq_end'])-np.float64(info['freq_start']))
	else:
		print 'Parameter',pname,'can not be found.'
