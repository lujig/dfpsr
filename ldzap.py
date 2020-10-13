#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import numpy.fft as fft
import argparse as ap
from matplotlib.figure import Figure
import matplotlib.lines as ln
import ld,os,copy
#
version='JigLu_20200923'
parser=ap.ArgumentParser(prog='ldzap',description='Zap the frequency domain interference in ld file.',epilog='Ver '+version)
parser.add_argument('-v','--version',action='version',version=version)
parser.add_argument("-z","--zap",dest="zap_file",default=0,help="file recording zap channels")
parser.add_argument("filename",help="input ld file")
args=(parser.parse_args())
#
if not os.path.isfile(args.filename):
	parser.error('A valid ld file name is required.')
d=ld.ld(args.filename)
info=d.read_info()
#
if 'compressed' in info.keys():
	nchan=int(info['nchan_new'])
	nbin=int(info['nbin_new'])
	nperiod=int(info['nsub_new'])
else:
	nchan=int(info['nchan'])
	if info['mode']=='test':
		nbin=1
		nperiod=int(d.read_shape()[1])
	else:
		nbin=int(info['nbin'])
		nperiod=int(info['nsub'])
npol=int(info['npol'])
if nbin!=1:
	data=d.period_scrunch()[:,:,0]
else:
	data=d.__read_bin_segment__(0,nperiod)[:,:,0]
if nbin>128 or ((nbin==1)&(nperiod>128)):
	data=fft.irfft(fft.rfft(data,axis=1)[:,:65],axis=1)
testdata=copy.deepcopy(data)
testdata=ma.masked_where(testdata<0,testdata)
if 'zchan' in info.keys():
	zaplist=[map(int,info['zchan'].split(','))]
	zapnum=zaplist[0]
	zaparray=np.zeros_like(testdata)
	zaparray[zapnum,:]=True
	testdata.mask=zaparray
	zap0=1
else:
	zaplist=[]
	zap0=0
#
if args.zap_file:
	if not os.path.isfile(args.zap_file):
		parser.error('The zap channel file is invalid.')
	zchan=np.loadtxt(args.zap_file,dtype=np.int32)
	if np.max(zchan)>=nchan or np.min(zchan)<0:
		parser.error('The zapped channel number is overrange.')
	zap0+=1
	zaplist.append(zchan)
	zapnum=set()
	for i in zaplist:
		zapnum.update(i)
	zapnum=np.array(list(zapnum))
	zaparray=np.zeros_like(testdata)
	zaparray[zapnum,:]=True
	testdata.mask=zaparray
#
spec=testdata.sum(1)
spec=spec-np.min(spec)
spec0=np.append(0,np.append(spec.repeat(2),0))
spec1=copy.deepcopy(spec0)
freq_start,freq_end=np.float64(info['freq_start']),np.float64(info['freq_end'])
ylim0=[freq_start,freq_end]
channelwidth=(freq_end-freq_start)/nchan
halfwidth=channelwidth/2
freq=np.linspace(ylim0[0]-halfwidth,ylim0[1]-halfwidth,len(spec)+1).repeat(2)
#
def plotimage(ylim):
	ax.imshow(testdata[::-1,:],aspect='auto',interpolation='nearest',extent=(0,1,ylim0[0]-halfwidth,ylim0[1]-halfwidth),cmap=colormap)
	ax1.plot(spec1,freq,'k-')
	ax.set_ylim(ylim[0]-halfwidth,ylim[1]-halfwidth)
	ax1.set_ylim(ylim[0]-halfwidth,ylim[1]-halfwidth)
	ax1.set_xlim(0,np.max(spec1)*1.1)
	ax1.set_xticks([])
	ax1.set_yticks([])
	if nbin==1:
		ax.set_xlabel('Integration Length',fontsize=30)
	else:
		ax.set_xlabel('Pulse Phase',fontsize=30)
	ax.set_ylabel('Frequency (MHz)',fontsize=30)
	canvas.draw()
#
def ycal(y,ylim):
	return (fig.bbox.extents[3]-y-ax.bbox.extents[1])/ax.bbox.bounds[3]*(ylim[1]-ylim[0])+ylim[0]
#
def chancal(y):
	return np.int32((y-freq_start)/channelwidth)
#
fig=Figure(figsize=(40,30),dpi=80)
fig.clf()
colormap='jet'
x0,x1,x2=0.13,0.8,0.95
y0,y1=0.11,0.96
ax=fig.add_axes([x0,y0,x1-x0,y1-y0])
ax1=fig.add_axes([x1,y0,x2-x1,y1-y0])
ylim=0
ylimlist=[]
l1 = ln.Line2D([0,1],[0.5,0.5],color='k',transform=fig.transFigure,figure=fig)
fig.lines.append(l1)
#
def leftclick(event):
	global ylim
	if event.x<ax.bbox.extents[0] or event.x>ax1.bbox.extents[2] or event.y<fig.bbox.extents[1] or event.y>fig.bbox.extents[3]: return
	if ylimlist:
		ylim1=ylimlist[-1]
	else:
		ylim1=ylim0
	if ylim==0:
		y=(fig.bbox.extents[3]-event.y)/fig.bbox.extents[3]
		l2=ln.Line2D([0,1],[y,y],color='k',transform=fig.transFigure,figure=fig)
		fig.lines.append(l2)
		canvas.draw()
		ylim=ycal(event.y,ylim1)
	else:
		ylim=[ylim,ycal(event.y,ylim1)]
		ylimlist.append(ylim)
		ax.set_ylim(ylim[0]-halfwidth,ylim[1]-halfwidth)
		ax1.set_ylim(ylim[0]-halfwidth,ylim[1]-halfwidth)
		fig.lines=[l1]
		canvas.draw()
		ylim=0
#
def rightclick(event):
	global ylim
	if event.x<ax.bbox.extents[0] or event.x>ax1.bbox.extents[2] or event.y<fig.bbox.extents[1] or event.y>fig.bbox.extents[3]: return
	if ylimlist:
		ylim1=ylimlist[-1]
	else:
		ylim1=ylim0
	chan=chancal(ycal(event.y,ylim1))
	if chan>=0 and chan<nchan:
		if ylim==0:
			zaplist.append([chan])
		else:
			if chancal(ylim)<=chan:
				zaplist.append(range(chancal(ylim),chan+1))
			else:
				zaplist.append(range(chan,chancal(ylim)+1))
			ylim=0
		update_image()
#
def update_image():
	global ylim,spec1
	zapnum=set()
	for i in zaplist:
		zapnum.update(i)
	zapnum=np.array(list(zapnum))
	zaparray=np.zeros_like(testdata)
	zaparray[zapnum,:]=True
	testdata.mask=zaparray
	if ylimlist:
		ylim=ylimlist[-1]
	else:
		ylim=ylim0
	ax.cla()
	ax1.cla()
	spec1=copy.deepcopy(spec0)
	spec1[2*zapnum+1]=0
	spec1[2*zapnum+2]=0
	fig.lines=[l1]
	plotimage(ylim)
	ylim=0
#
def move_tk(event):
    if event.x>ax.bbox.extents[0] and event.x<ax1.bbox.extents[2]: 
		y=(fig.bbox.extents[3]-event.y)/fig.bbox.extents[3]
		l1.set_ydata([y,y])
		canvas.draw()
#
def move_gtk(window,event):
    if event.x>ax.bbox.extents[0] and event.x<ax1.bbox.extents[2]: 
		y=(fig.bbox.extents[3]-event.y)/fig.bbox.extents[3]
		l1.set_ydata([y,y])
		canvas.draw()
#
def press_gtk(window,event):
	keymotion(gtk.gdk.keyval_name(event.keyval))
#
def press_tk(event):
	keymotion(event.keysym)
#
def keymotion(a):
	global ylim
	if a=='q':
		root.destroy()
		zapnum=set()
		for i in zaplist:
			zapnum.update(i)
		np.savetxt('zap.txt',np.sort(list(zapnum)),fmt='%i')
	elif a=='s':
		root.destroy()
		zapnum=set()
		for i in zaplist:
			zapnum.update(i)
		zapnum=list(zapnum)
		info['zchan']=str(zapnum)[1:-1]
		save=ld.ld('.'.join(args.filename.split('.')[:-1])+'_zap.ld')
		save.write_shape([nchan,nperiod,nbin,npol])
		for i in np.arange(nchan):
			if i in zapnum:
				save.write_chan(np.zeros(nperiod*nbin),i)
				continue
			save.write_chan(d.read_chan(i),i)
		save.write_info(info)
	elif a=='r':
		if ylimlist:
			ylimlist.pop()
		else: return
		if ylimlist:
			ylim=ylimlist[-1]
		else:
			ylim=ylim0
		ax.set_ylim(ylim[0]-halfwidth,ylim[1]-halfwidth)
		ax1.set_ylim(ylim[0]-halfwidth,ylim[1]-halfwidth)
		canvas.draw()
		ylim=0
	elif a=='u':
		if len(zaplist)>zap0:
			zaplist.pop()
		else: return
		update_image()
	elif a=='h':
		print "\nldzap interactive commands\n"
		print "Mouse:"
		print "  Left-click selects the start of a range"
		print "    then left-click again to zoom, or right-click to zap."
		print "  Right-click zaps current cursor location.\n"
		print "Keyboard:"
		print "  h  Show this help"
		print "  u  Undo last zap command"
		print "  r  Reset zoom and update dynamic spectrum"
		print "  s  Save zapped version as (filename)_zap.ld and quit"
		print "  q  Exit program\n"
#
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import Tkinter as tk
import matplotlib as mpl
mpl.use('TkAgg')
root=tk.Tk()
root.title(args.filename)
root.geometry('800x600+100+100')
canvas=FigureCanvasTkAgg(fig,master=root)
canvas.get_tk_widget().grid()  
canvas.get_tk_widget().pack(fill='both')
root.bind('<KeyPress>',press_tk)
root.bind('<ButtonPress-1>',leftclick)
root.bind('<ButtonPress-3>',rightclick)
root.bind('<Motion>',move_tk)
canvas.draw()
plotimage(ylim0)
root.mainloop()
#except:
	#raise(Exception,"Only Tkinter available for this programme.")
	#import gtk
	#from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg
	#root=gtk.Window()
	#root.set_title(args.filename)
	#root.set_size_request(800,600)
	#box=gtk.VBox()
	#canvas=FigureCanvasGTKAgg(fig)
	#box.pack_start(canvas)
	#root.add(box)
	#root.modify_bg('normal',gtk.gdk.Color('#fff'))
	#root.show_all()
	#root.connect('destroy',gtk.main_quit)
	#root.connect('key-press-event',press_gtk)
	#root.connect('move-cursor',move_gtk)
	#gtk.main()

