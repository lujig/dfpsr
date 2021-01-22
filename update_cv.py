import urllib.request as ur
import read as rd
import numpy as np
#
# utc2tai, leap second
a=ur.urlopen('https://hpiers.obspm.fr/eop-pc/earthor/utc/UTC-offsets_tab.html').readlines()
b=list(map(lambda x:x.strip().decode(),a))
lenb=len(b)
for i in range(lenb):
	num=lenb-1-i
	t=b[num]
	if t:
		if not t[0:4].isdigit(): 
			b.pop(num)
		elif int(t[0:4])<1972:
			b.pop(num)
		else:
			b[num]=t.replace('Jan.','1').replace('Jul.','7').replace('1s','').strip()+'1\n'
			if b[num].split()[:2]==['1972','1']:
				b.pop(num)
	else:
		b.pop(num)
f=open('convensions/'+'leap.txt','w')
f.writelines(b)
f.close()
#
# utc2ut1
a=ur.urlopen('https://hpiers.obspm.fr/eoppc/eop/eopc04/eopc04.62-now').readlines()
b=list(map(lambda x:x.strip().decode(),a))
lenb=len(b)
for i in range(lenb):
	num=lenb-1-i
	t=b[num]
	if t:
		if not t[0:4].isdigit(): 
			b.pop(num)
		else:
			b[num]=t+'\n'
	else:
		b.pop(num)
f=open('convensions/'+'eopc.txt','w')
f.writelines(b)
f.close()
#
# tai2ut1
mjd,deltat=np.float64(list(map(lambda x:x.split(),b)))[:,[3,6]].T
taimjd=rd.utc2tai(mjd)
f=np.loadtxt('convensions/'+'leap.txt')
leap_time0=np.array(list(map(lambda x:rd.datetime2mjd([x[0],x[1],x[2],0,0,0]),f))).reshape(-1)
leap_time=np.array(list(map(lambda x:f[leap_time0<=x,3].sum(),mjd)))
deltat+=leap_time
f=open('convensions/'+'tai2ut1.txt','w')
f.writelines(list(map(lambda x: str(x[0])+' '+str(x[1])+' '+str(x[2])+'\n',zip(mjd,taimjd,deltat))))
f.close()
#
# tai2tt
a=ur.urlopen('ftp://ftp2.bipm.org/pub/tai/ttbipm/TTBIPM.2019').readlines()
b=list(map(lambda x:x.decode(),a))
lenb=len(b)
for i in range(lenb):
	num=lenb-1-i
	t=b[num]
	if t:
		if not t[0:5].isdigit(): 
			b.pop(num)
	else:
		b.pop(num)
f=open('convensions/'+'tai2tt.txt','w')
f.writelines(b)
f.close()
#
# gps2utc
a=ur.urlopen('ftp://ftp2.bipm.org/pub/tai/other-products/utcgnss/utc-gnss').readlines()[20:]
b=list(map(lambda x:x.decode().strip(),a))
lenb=len(b)
for i in range(lenb):
	num=lenb-1-i
	t=b[num]
	if t:
		if not t[0:5].isdigit(): 
			b.pop(num)
		else:
			b[num]=t[:23]+'\n'
	else:
		b.pop(num)
f=open('convensions/'+'gps2utc.txt','w')
f.writelines(b)
f.close()













