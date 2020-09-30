#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import numpy.fft as fft
from matplotlib.figure import Figure
import argparse as ap
import os,time,ld
import warnings as wn
#
version='JigLu_20180508'
parser=ap.ArgumentParser(prog='ldplt',description='Plot the ld file. Press \'s\' in figure window to save figure.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",help="input ld file")
parser.add_argument('-f',action='store_true',default=False,dest='fdomain',help='show the frequency domain image')
parser.add_argument('-t',action='store_true',default=False,dest='tdomain',help='show the time domain image')
parser.add_argument('-p',action='store_true',default=False,dest='profile',help='show the pulse profile')
parser.add_argument('-b','--phase_range',default=0,dest='phase',help='limit the phase range, PHASE0,PHASE1')
parser.add_argument('-r','--frequency_range',default=0,dest='frequency',help='limit the frequency rangeFREQ0,FREQ1')
parser.add_argument('-s','--subint_range',default=0,dest='subint',help='limit the subint range SUBINT0,SUBINT1')
parser.add_argument('-o','--polynomial_order',default=0,dest='n',type=int,help='fit the back ground with Nth order polynomial')
parser.add_argument('--polar',default=0,dest='polar',type=int,help='plot the specified polarization')
parser.add_argument('-c','--rotation',default=0,dest='rotation',type=np.float64,help='rotate the plot phase')
parser.add_argument('-n',action='store_true',default=False,dest='norm',help='normalized the data at each channel or subint')
parser.add_argument('-i',action='store_false',default=True,dest='title',help='hide file information above the figure')
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
	npol=int(info['npol_new'])
else:
	nchan=int(info['nchan'])
	nbin=int(info['nbin'])
	nsub=int(info['nperiod'])
	npol=int(info['npol'])
freq_start=np.float64(info['freq_start'])
freq_end=np.float64(info['freq_end'])
freq=(freq_start+freq_end)/2.0
bw=freq_end-freq_start
channel_width=(freq_end-freq_start)/nchan
#
plotflag=np.sum(map(np.bool,[args.fdomain,args.tdomain,args.profile]))
if plotflag>1:
	parser.error('At most one of flags -f, -t and -p is required.')
elif plotflag==0:
	parser.error('At least one of flags -f, -t and -p is required.')
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
def shift(y,x):
	fftp=fft.rfft(y,axis=0)
	ffts=fftp*np.exp(-2*np.pi*x*1j*np.arange(np.shape(fftp)[-1])).reshape(-1,1)
	fftr=fft.irfft(ffts,axis=0)
	return fftr
#
fig=Figure(figsize=(40,30),dpi=80)
fig.set_facecolor('white')
ax=fig.add_axes([0.12,0.1,0.82,0.83])
ax.patch.set_facecolor('w')
if args.fdomain:
	data=d.period_scrunch(subint_start,subint_end,chan)[:,:,0]
	if 'zchan' in info.keys():
		if len(chan):
			zchan=np.array(list(set(np.int32(info['zchan'].split(','))).intersection(chan)))-chanstart
		else:
			zchan=np.int32(info['zchan'].split(','))
		zaparray=np.zeros_like(data)
		zaparray[zchan]=True
		data=ma.masked_array(data,mask=zaparray)
	if args.n:
		data-=np.polyval(np.polyfit(np.arange(nbin),data.T,args.n),np.array([range(nbin)]*len(data)).T).T
	else:
		data-=data.mean(1).reshape(-1,1)
	if args.norm:
		data/=data.max(1).reshape(-1,1)
	if args.rotation:
		data=shift(data,args.rotation)
	ax.imshow(data[::-1],aspect='auto',interpolation='nearest',extent=(0,1,freq_start,freq_end),cmap='jet')
	ax.set_ylabel('Frequency (MHz)',fontsize=30)
	ax.set_ylim(frequency[0],frequency[1])
	texty=frequency[1]*1.01-frequency[0]*0.01
if args.tdomain:
	data=d.chan_scrunch(chan,subint_start,subint_end)[:,:,0]
	if args.n:
		data-=np.polyval(np.polyfit(np.arange(nbin),data.T,args.n),np.array([range(nbin)]*len(data)).T).T
	else:
		data-=data.mean(1).reshape(-1,1)
	if args.norm:
		data/=data.max(1).reshape(-1,1)
	if args.rotation:
		data=shift(data,args.rotation)
	ax.imshow(data[::-1],aspect='auto',interpolation='nearest',extent=(0,1,subint_start,subint_end),cmap='jet')
	ax.set_ylabel('Subint Number',fontsize=30)
	ax.set_ylim(subint[0],subint[1])
	texty=subint[1]*1.01-subint[0]*0.01
if args.profile:
	data=d.chan_scrunch(chan,subint_start,subint_end).sum(0)
	if args.n:
		data-=np.polyval(np.polyfit(np.arange(nbin),data,args.n),np.arange(nbin))
	bins,mn=10,nbin/10
	stat,val,tmp=ax.hist(data[:,0],bins)
	while stat.max()>mn and 2*bins<mn:
		bins*=2
		stat,val,tmp=ax.hist(data[:,0],bins)
	val=(val[1:]+val[:-1])/2.0
	argmaxstat=np.argmax(stat)
	if argmaxstat==0:
		base=data[(data[:,0]>(val[1]*-0.5+val[0]*1.5))&(data[:,0]<(val[1]*0.5+val[0]*0.5))].mean(0)
	elif argmaxstat==1:
		poly=np.polyfit(val[:3],stat[:3],2)
		#base=-poly[1]/poly[0]/2.0
		base=data[(data[:,0]>(-poly[1]/poly[0]/2.0-(val[1]-val[0])*0.5))&(data[:,0]<(-poly[1]/poly[0]/2.0+(val[1]-val[0])*0.5))].mean(0)
	else:
		poly=np.polyfit(val[(argmaxstat-2):(argmaxstat+3)],stat[(argmaxstat-2):(argmaxstat+3)],2)
		#base=-poly[1]/poly[0]/2.0
		base=data[(data[:,0]>(-poly[1]/poly[0]/2.0-(val[1]-val[0])*0.5))&(data[:,0]<(-poly[1]/poly[0]/2.0+(val[1]-val[0])*0.5))].mean(0)
	data-=base
	data/=np.max(data[:,0])
	if args.rotation:
		data=shift(data,args.rotation)
	low=min(data[:,0])*1.1-max(data[:,0])*0.1
	high=max(data[:,0])*1.1-min(data[:,0])*0.1
	x=np.linspace(0,1,len(data))
	ax.plot(x,data,'k-')
	ax.set_ylabel('Flux (Arbitrary Unit)',fontsize=30)
	ax.set_ylim(low,high)
	texty=high*1.01-low*0.01
#
ax.set_xlabel('Pulse Phase',fontsize=30)
ax.set_xlim(phase[0],phase[1])
if args.title:
	ax.text(0.0,texty,'Freq: '+str(np.round(freq,1))+' BW: '+str(np.round(bw,1)),horizontalalignment='left',verticalalignment='bottom',fontsize=15)
	ax.text(1.0,texty,'Length: '+str(np.round(np.float64(info['length']),1)),horizontalalignment='right',verticalalignment='bottom',fontsize=15)
	if 'psr_name' in info.keys():
		ax.text(0.5,texty,info['psr_name'],horizontalalignment='center',verticalalignment='bottom',fontsize=20)
	else:
		ax.text(0.5,texty,args.filename,horizontalalignment='center',verticalalignment='bottom',fontsize=20)		
#
def save_fig():
	figname=raw_input("Please input figure name:")
	if figname.split('.')[-1] not in ['ps','eps','png','pdf','pgf']:
		figname+='.pdf'
	fig.savefig(figname)
	print 'Figure file',figname,'has been saved.'
#
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
	def save_gtk(window,event):
		if gtk.gdk.keyval_name(event.keyval)=='s':
			save_fig()
	window.connect('key-press-event',save_gtk)
	gtk.main()
except:
	import matplotlib as mpl
	from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
	import Tkinter as tk
	mpl.use('TkAgg')
	ax.tick_params(axis='x',labelsize=15)
	ax.tick_params(axis='y',labelsize=15)
	root=tk.Tk()
	root.title(args.filename)
	root.geometry('800x600+100+100')
	canvas=FigureCanvasTkAgg(fig,master=root)
	canvas.get_tk_widget().grid()
	canvas.get_tk_widget().pack(fill='both')
	canvas.show()
	def save_tk(event):
		if event.keysym=='s':
			save_fig()
	root.bind('<KeyPress>',save_tk)
	root.mainloop()
