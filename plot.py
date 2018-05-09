#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import argparse as ap
import os,time,ld
import warnings as wn
#
version='JigLu_20180508'
parser=ap.ArgumentParser(prog='plot',description='Plot the ld file.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",help="input ld file")
parser.add_argument('-f',action='store_true',default=False,dest='fdomain',help='show the frequency domain image')
parser.add_argument('-t',action='store_true',default=False,dest='tdomain',help='show the time domain image')
parser.add_argument('-p',action='store_true',default=False,dest='profile',help='show the pulse profile')
parser.add_argument('-b','--phase_range',default=0,dest='phase',help='limit the phase range, PHASE0,PHASE1')
parser.add_argument('-r','--frequency_range',default=0,dest='frequency',help='limit the frequency rangeFREQ0,FREQ1')
parser.add_argument('-s','--subint_range',default=0,dest='subint',help='limit the subint range SUBINT0,SUBINT1')
args=(parser.parse_args())
wn.filterwarnings('ignore')
#
if not os.path.isfile(args.filename):
	parser.error('A valid ld file name is required.')
d=ld.ld(args.filename)
info=d.read_info()
if 'compressed' in info.keys():
	nchan=int(info['nchan_new'])
	nbin=int(info['nbin_new'])
	nsub=int(info['nsub_new'])
else:
	nchan=int(info['nchan'])
	nbin=int(info['nbin'])
	nsub=int(info['nperiod'])
freq_start=np.float64(info['freq_start'])
freq_end=np.float64(info['freq_end'])
freq=(freq_start+freq_end)/2.0
channel_width=(freq_end-freq_start)/nchan
#
plotflag=np.sum(map(np.bool,[args.fdomain,args.tdomain,args.profile]))
if plotflag>1:
	parser.error('At most one of flags -f, -t and -p is required.')
elif plotflag==0:
	parser.error('At least one of flags -f, -t and -p is required.')
#
if args.phase:
	phase=np.float64(args.phase.split(','))
else:
	phase=np.arange(2)
#
if args.frequency:
	frequency=np.float64(args.frequency.split(','))
	if len(frequency)!=2:
		parser.error('A valid frequency range should be given.')
	if frequency[0]>frequency[1]:
		parser.error("Starting frequency larger than ending frequency.")
	freq_start=max(frequency[0],freq_start)
	freq_end=min(frequency[1],freq_end)
	chanstart,chanend=np.int16(np.round((np.array([freq_start,freq_end])-freq)/channel_width+0.5*nchan))
	chan=np.arange(chanstart,chanend)
	if len(chan)==0:
		parser.error('Input bandwidth is too narrow.')
else:
	frequency=np.array([freq_start,freq_end])
	chan=[]
#
if args.subint:
	subint=np.float64(args.subint.split(','))
	if len(subint)!=2:
		parser.error('A valid subint range should be given.')
	if subint[0]>subint[1]:
		parser.error("Starting subint larger than ending subint.")
	subint_start=max(int(subint[0]),0)
	subint_end=min(int(subint[1]+1),nsub)
else:
	subint_start=0
	subint_end=nsub
	subint=np.array([subint_start,subint_end])
#
if args.phase:
	phase=np.float64(args.phase.split(','))
	if len(phase)!=2:
		parser.error('A valid phase range should be given.')
	if phase[0]>phase[1]:
		parser.error("Starting phase larger than ending phase.")
else:
	phase=np.array([0,1])
#
fig=Figure(figsize=(40,30))
fig.set_facecolor('white')
ax=fig.add_subplot(111)
if args.fdomain:
	data=d.period_scrunch(subint_start,subint_end,chan).sum(2)
	ax.imshow(data,aspect='auto',interpolation='nearest',extent=(0,1,freq_start,freq_end),cmap='jet')
	ax.set_xlabel('Pulse Phase',fontsize=15)
	ax.set_ylabel('Frequency (MHz)',fontsize=15)
	ax.set_xlim(phase[0],phase[1])
	ax.set_ylim(frequency[0],frequency[1])
if args.tdomain:
	data=d.chan_scrunch(chan,subint_start,subint_end).sum(2)
	ax.imshow(data,aspect='auto',interpolation='nearest',extent=(0,1,subint_start,subint_end),cmap='jet')
	ax.set_xlabel('Pulse Phase',fontsize=15)
	ax.set_ylabel('Subint Number',fontsize=15)
	ax.set_xlim(phase[0],phase[1])
	ax.set_ylim(subint[0],subint[1])
if args.profile:
	data=d.chan_scrunch(chan,subint_start,subint_end).sum(2).sum(0)
	low=min(data)*1.1-max(data)*0.1
	high=max(data)*1.1-min(data)*0.1
	x=np.linspace(0,1,nbin)
	ax.plot(x,data,'k-')
	ax.set_xlabel('Pulse Phase',fontsize=15)
	ax.set_ylabel('Flux (Arbitrary Unit)',fontsize=15)
	ax.set_yticks([])
	ax.set_xlim(phase[0],phase[1])
	ax.set_ylim(low,high)
try:
	import gtk
	from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg
	window=gtk.Window()
	window.set_title(args.filename)
	window.set_size_request(800,600)
	box=gtk.VBox()
	canvas=FigureCanvasGTKAgg(fig)
	box.pack_start(canvas)
	window.add(box)
	window.modify_bg('normal',gtk.gdk.Color('#fff'))
	window.show_all()
	window.connect('destroy',gtk.main_quit)
	gtk.main()
except:
	import matplotlib as mpl
	from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
	import Tkinter as tk
	mpl.use('TkAgg')
	root=tk.Tk()
	root.title(args.filename)
	root.geometry('800x600+100+100')
	canvas=FigureCanvasTkAgg(fig,master=root)
	canvas.get_tk_widget().grid()
	canvas.get_tk_widget().pack(fill='both')
	canvas.show()
	root.mainloop()
