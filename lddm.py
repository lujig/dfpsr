#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import numpy.fft as fft
from matplotlib.figure import Figure
import argparse as ap
import os,time,ld
import warnings as wn
import matplotlib.pyplot as plt
#
version='JigLu_20180515'
parser=ap.ArgumentParser(prog='lddm',description='Calculate the best DM value. Press \'s\' in figure window to save figure.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("filename",help="input ld file")
parser.add_argument('-r','--frequency_range',default=0,dest='frequency',help='limit the frequency rangeFREQ0,FREQ1')
parser.add_argument('-s','--subint_range',default=0,dest='subint',help='limit the subint range SUBINT0,SUBINT1')
parser.add_argument('-n',action='store_true',default=False,dest='norm',help='normalized the data at each channel before optimization')
parser.add_argument('-d','--dm_center',dest='dm',default=0,type=np.float64,help="dispersion measure")
parser.add_argument('-z','--dm_zone',dest='zone',default=0,type=np.float64,help="dispersion measure")
parser.add_argument('-o','--polynomial_order',default=0,dest='n',type=int,help='fit the dm-maxima curve with Nth order polynomial')
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
dm0=np.float64(info['dm'])
period=np.float64(info['period'])
if args.dm:
	ddm=args.dm-dm0
else:
	ddm=0
#	
length=100
if args.zone:
	zone=args.zone/2
else:
	zone=np.max([0.1,dm0/100])
	zone=np.min([0.5,zone])
dm=np.linspace(ddm-zone,ddm+zone,length)
freq_start0=np.float64(info['freq_start'])
freq_end0=np.float64(info['freq_end'])
channel_width=(freq_end0-freq_start0)/nchan
freq0=np.arange(freq_start0,freq_end0,channel_width)+channel_width/2.0
if args.frequency:
	frequency=np.float64(args.frequency.split(','))
	if len(frequency)!=2:
		parser.error('A valid frequency range should be given.')
	if frequency[0]>frequency[1]:
		parser.error("Starting frequency larger than ending frequency.")
	freq_start=max(frequency[0],freq_start0)
	freq_end=min(frequency[1],freq_end0)
	chanstart,chanend=np.int16(np.round((np.array([freq_start,freq_end])-(freq_start+freq_end)/2.0)/channel_width+0.5*nchan))
	chan=np.arange(chanstart,chanend)
	if len(chan)==0:
		parser.error('Input bandwidth is too narrow.')
	freq=freq0[chan]
else:
	freq=freq0
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
data0=d.period_scrunch(subint_start,subint_end).sum(2)
if len(chan):
	data=data0[chan]
else:
	data=data0.copy()
#
fig=Figure(figsize=(40,30),dpi=80)
fig.set_facecolor('white')
x0,x1,x2=0.1,0.6,0.9
y0,y1=0.11,0.96
ax=fig.add_axes([x0,y0,x1-x0,y1-y0])
ax.patch.set_facecolor('w')
ax1=fig.add_axes([x1,y0,x2-x1,(y1-y0)/2])
ax1.set_xlabel('Pulse Phase',fontsize=25)
ax.set_yticks([])
ax1.set_yticks([])
ax1=ax1.twinx()
ax2=fig.add_axes([x1,(y1+y0)/2,x2-x1,(y1-y0)/2])
ax2.set_xticks([])
ax2.set_yticks([])
ax2=ax2.twinx()
ax1.patch.set_facecolor('w')
ax2.patch.set_facecolor('w')
ax.set_xlabel('DM',fontsize=30)
ax.set_ylabel('Relative Maxima',fontsize=30)
ax1.set_ylabel('Frequency (MHz)',fontsize=25)
ax2.set_ylabel('Frequency (MHz)',fontsize=25)
#
if 'zchan' in info.keys():
	if len(chan):
		zchan=np.array(list(set(np.int32(info['zchan'].split(','))).intersection(chan)))-chanstart
	else:
		zchan=np.int32(info['zchan'].split(','))
	zaparray=np.zeros_like(data0)
	zaparray[zchan]=True
	data0=ma.masked_array(data0,mask=zaparray)
#
data-=data.mean(1).reshape(-1,1)
data0-=data0.mean(1).reshape(-1,1)
if args.norm:
	maxima=data.max(1)
	maxima[maxima==0]=1
	data/=maxima.reshape(-1,1)
	maxima=data0.max(1)
	maxima[maxima==0]=1
	data0/=maxima.reshape(-1,1)
#
#
ax2.imshow(data0[::-1],aspect='auto',interpolation='nearest',extent=(0,1,freq_start0,freq_end0),cmap='jet')
if args.frequency:
	ax2.plot([0,1],[freq_start,freq_start],'k--')
	ax2.plot([0,1],[freq_end,freq_end],'k--')
#
def shift(y,x):
	ffts=y*np.exp(x*1j)
	fftr=fft.irfft(ffts)
	return fftr
#
fftdata=fft.rfft(data,axis=1)
tmp=np.shape(fftdata)[-1]
const=(1/freq**2*4148.808/period*np.pi*2.0).repeat(tmp).reshape(-1,tmp)*np.arange(tmp)
#
value=np.zeros(length)
for i in np.arange(length):
	disp=dm[i]*const
	value[i]=np.max(shift(fftdata,disp).sum(0))
ax.plot(dm+dm0,value,'b-')
x0=dm.min()+dm0
x1=dm.max()+dm0
y0=value.min()*1.1-value.max()*0.1
y1=value.max()*1.1-value.min()*0.1
ax.set_xlim(x0,x1)
ax.set_ylim(y0,y1)
#
if args.n:
	order=args.n
else:
	order=7
polyfunc=np.polyfit(dm,value,order)
fitvalue=np.polyval(polyfunc,dm)
roots=np.roots(np.polyder(polyfunc))
roots=np.real(roots[np.isreal(roots)])
ax.plot(dm+dm0,fitvalue,'k--')
if len(roots):
	dmmax=roots[np.argmin(np.abs(roots-dm[np.argmax(value)]))]
	if dmmax<dm[-1] and dmmax>dm[0] and np.polyval(np.polyder(polyfunc,2),dmmax)<0:
		error=np.std(value-fitvalue)
		errfunc=np.append(polyfunc[:-1],polyfunc[-1]-np.polyval(polyfunc,dmmax)+error)
		roots=np.roots(errfunc)
		roots=np.real(roots[np.isreal(roots)])
		dmerr=np.mean(np.abs([roots[np.argmin(np.abs(roots-dmmax+zone/10))],roots[np.argmin(np.abs(roots-dmmax-zone/10))]]-dmmax))
		ax.plot([dm0+dmmax,dm0+dmmax],[y0,y1],'k:')
		ax.text(dm0+ddm,y0*0.95+y1*0.05,'DM$_0$='+str(dm0)+'\nBest DM='+str(np.round(dmmax+dm0,3))+'$\pm$'+str(np.round(dmerr,3)),horizontalalignment='center',verticalalignment='bottom',fontsize=25)
		fftdata=fft.rfft(data0,axis=1)
		tmp=np.shape(fftdata)[-1]
		frac=1/freq0**2*4148.808/period*dmmax
		const=(frac*np.pi*2.0).repeat(tmp).reshape(-1,tmp)*np.arange(tmp)
		data0=shift(fftdata,const)
		if 'zchan' in info.keys():
			data0=ma.masked_array(data0,mask=zaparray)
		ax1.imshow(data0[::-1],aspect='auto',interpolation='nearest',extent=(0,1,freq_start0,freq_end0),cmap='jet')
		ax1.plot(np.ones_like(freq0)*0.5,freq0,'r--')
		ax2.plot((frac+0.5)%1,freq0,'r--')
		if args.frequency:
			ax1.plot([0,1],[freq_start,freq_start],'k--')
			ax1.plot([0,1],[freq_end,freq_end],'k--')
	else:
		ax.text(dm0+ddm,y0*0.95+y1*0.05,'The best DM cannot be found',horizontalalignment='center',verticalalignment='bottom',fontsize=25)
else:
	ax.text(dm0+ddm,y0*0.95+y1*0.05,'The best DM cannot be found',horizontalalignment='center',verticalalignment='bottom',fontsize=25)
#
ax2.set_ylim(freq_start0,freq_end0)
ax2.set_xlim(0,1)
ax1.set_ylim(freq_start0,freq_end0)
ax1.set_xlim(0,1)
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
	window.set_size_request(1000,600)
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
	root.geometry('1000x600+100+100')
	canvas=FigureCanvasTkAgg(fig,master=root)
	canvas.get_tk_widget().grid()
	canvas.get_tk_widget().pack(fill='both')
	canvas.show()
	def save_tk(event):
		if event.keysym=='s':
			save_fig()
	root.bind('<KeyPress>',save_tk)
	root.mainloop()
