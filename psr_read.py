import numpy as np
import time_eph as te
import subprocess as sp
#
class psr:
	def __init__(self,name):
		self.name=name
		self.readpara()
		if self.units=='TDB':
			self.change_units()
		self.cal_posvel()
	#
	def cal_posvel(self):
		alpha=self.raj
		delta=self.decj
		ca=np.cos(alpha)
		sa=np.sin(alpha)
		cd=np.cos(delta)
		sd=np.sin(delta)
		epsilon = 23.439292*np.pi/180 #???
		sinb=sd*np.cos(epsilon)-cd*np.sin(epsilon)*sa
		self.elat=np.arcsin(sinb)
		y=sa*np.cos(epsilon)+sd/cd*np.sin(epsilon)
		self.elong=np.arctan2(y,ca)%(np.pi*2)
		self.pos_equ=te.vector(ca*cd,sa*cd,sd,center='bary',scale='si',coord='equ',unit=te.sl,type0='pos')
		self.pos=self.pos_equ.copy()
		self.pos.equ2ecl()
		convert = 1.0/1000.0/60.0/60.0*np.pi/180.0*100.0
		if hasattr(self,'pmra'):
			# note that the mas/yr in PMRA is not the true mas, but to be the mas of RA. similar in PMRA2
			self.vel_equ=te.vector(convert*(-self.pmra*sa-self.pmdec*ca*sd),convert*(self.pmra*ca-self.pmdec*sa*sd),convert*self.pmdec*cd,center='bary',scale='si',coord='equ',unit=te.sl,type0='vel')
			self.vel=self.vel_equ.copy()
			self.vel.equ2ecl()
		if hasattr(self,'pmra2'):
			convert = 1.0/1000.0/60.0/60.0*np.pi/180.0*100.0*100.0
			self.acc_equ=te.vector(convert*(-self.pmra2*sa-self.pmdec2*ca*sd),convert*(self.pmra2*ca-self.pmdec2*sa*sd),convert*self.pmdec2*cd,center='bary',scale='si',coord='equ',unit=te.sl,type0='vel') # something wrong???
			self.acc=self.acc_equ.copy()
			self.acc.equ2ecl()
	#
	def change_units(self):
		if self.units=='TCB':
			return
		for i in self.paras:
			val=self.__getattribute__(i)
			if i in paras_p1:
				factor=np.float64(te.iftek)
			elif i in paras_m1:
				factor=1/np.float64(te.iftek)
			elif i in paras_m2:
				factor=1/np.float64(te.iftek)**2
			elif i in paras_m3:
				factor=1/np.float64(te.iftek)**3
			elif i in paras_m4:
				factor=1/np.float64(te.iftek)**4
			elif i in paras_time:
				factor=1
				self.__setattr__(i,val.add(val.tdb2tcb(),scale='tcb'))
			elif i in paras_time_array:
				factor=1
				val_new=[]
				for k in val:
					val_new.append(k.add(k.tdb2tcb(),scale='tcb'))
				self.__setattr__(i,val_new)
			elif i in paras_eph:
				factor=1
				self.__setattr__(i,val+te.time(val,np.zeros_like(val),scale='tdb').tdb2tcb().mjd)
			elif i in paras_mix:
				factor=1/np.float64(te.iftek)**(paras_mix[i]-np.arange(0,len(val)))
			else:
				continue
			self.__setattr__(i,self.__getattribute__(i)*factor)
	#
	def __str__(self):
		string=''
		for i in self.paras:
			if i in para_with_err:
				val=self.__getattribute__(i)
				err=self.__getattribute__(i+'_err')
				if type(val) is np.ndarray:
					for k in np.arange(val.size):
						if k==0: string+='{:12s} '.format(i.upper())
						else: string+='{:12s} '.format('')
						val_str=str(val[k])
						err_str=str(err[k])
						string+='{:25s} {:25s}'.format(val_str,err_str)+'\n'
				else:
					val_str=str(val)
					err_str=str(err)
					string+='{:12s} {:25s} {:25s}'.format(i.upper(),val_str,err_str)+'\n'
			else:
				val=self.__getattribute__(i)
				if type(val)==te.time:
					val=val.mjd
				if type(val) is np.ndarray:
					for k in np.arange(val.size):
						if k==0: string+='{:12s} '.format(i.upper())
						else: string+='{:12s} '.format('')
						val_str=str(val[k])
						string+='{:25s} '.format(val_str)+'\n'
				else:
					val_str=str(val)
					string+='{:12s} {:25s}'.format(i.upper(),val_str)+'\n'
		return string
	#
	def __repr__(self):
		return self.__str__()
	#
	def readpara(self):
		paras=sp.getoutput('psrcat -e '+self.name).split('\n')
		paras=dict(list(map(lambda x: [x[0],np.array(x[1].split())],map(lambda x: x.split(None,1),paras))))
		paras_key=paras.keys()
		self.paras=[]
		if 'PSRJ' in paras_key:
			self.name=paras['PSRJ'][0]
			self.paras.append('name')
		else:
			raise Exception('No pulsar name in par file.')
		#
		if 'UNITS' in paras_key:
			i=paras['UNITS']
			self.units=i[0].lower()
			self.paras.append('units')
		else:
			raise Exception('No para units in par file.')
		#
		if 'RAJ' in paras_key:
			i=paras['RAJ']
			self.raj=(np.float64(i[0].split(':'))*np.array([3600,60,1])).sum()/3600/12*np.pi
			self.raj_err=np.float64(i[1])/3600/12*np.pi
			self.paras.append('raj')
		else:
			raise Exception('No pulsar RA in par file.')
		#
		if 'DECJ' in paras_key:
			i=paras['DECJ']
			if i[0][0]=='-':
				self.decj=(np.float64(i[0].split(':'))*np.array([3600,-60,-1])).sum()/3600/12*np.pi
			else:
				self.decj=(np.float64(i[0].split(':'))*np.array([3600,60,1])).sum()/3600/12*np.pi
			self.decj_err=np.float64(i[1])/3600/12*np.pi
			self.paras.append('decj')
		else:
			raise Exception('No pulsar DEC in par file.')
		#
		if 'F0' in paras_key:
			i=paras['F0']
			self.f0=np.float64(i[0])
			self.f0_err=np.float64(i[1])
			self.p0=1/self.f0
			self.p0_err=self.f0_err/self.f0**2
			self.paras.extend(['f0','p0'])
		elif 'P0' in paras_key:
			i=paras['P0']
			self.p0=np.float64(i[0])
			self.p0_err=np.float64(i[1])
			self.f0=1/self.p0
			self.f0_err=self.p0_err/self.p0**2
			self.paras.extend(['f0','p0'])
		else:
			raise Exception('No pulsar period in par file.')
		#
		if 'F1' in paras_key:
			i=paras['F1']
			self.f1=np.float64(i[0])
			self.f1_err=np.float64(i[1])
			self.p1=-self.f1*self.p0**2
			self.p1_err=np.sqrt((self.f1_err*self.p0**2)**2+(2*self.f1*self.p0*self.p0_err)**2)
			self.paras.extend(['f1','p1'])
		elif 'P1' in paras_key:
			i=paras['P1']
			self.p1=np.float64(i[0])
			self.p1_err=np.float64(i[1])
			self.f1=-self.p1*self.f0**2
			self.f1_err=np.sqrt((self.p1_err*self.f0**2)**2+(2*self.p1*self.f0*self.f0_err)**2)
			self.paras.extend(['f1','p1'])
		else:
			self.f1=0
			self.p1=0
		#
		if 'F2' in paras_key:
			if 'f1' not in self.paras:
				raise Exception('The parameter F2 is in parfile without F1.')
			i=paras['F2']
			self.f2=np.float64(i[0])
			self.f2_err=np.float64(i[1])
			self.p2=-self.f2*self.p0**2-2*self.f1*self.p0*self.p1
			self.p2_err=np.sqrt((self.f2_err*self.p0**2)**2+((self.f2*self.p0+self.f1*self.p1)*2*self.p0_err)**2+(2*self.f1*self.p0*self.p1_err)**2+(2*self.p0*self.p1*self.f1_err)**2)
			self.paras.extend(['f2','p2'])
		elif 'P2' in paras_key:
			if 'f1' not in self.paras:
				raise Exception('The parameter P2 is in parfile without P1.')
			i=paras['P2']
			self.p2=np.float64(i[0])
			self.p2_err=np.float64(i[1])
			self.f2=-self.p2*self.f0**2-2*self.p1*self.f0*self.f1
			self.f2_err=np.sqrt((self.p2_err*self.f0**2)**2+((self.p2*self.f0+self.p1*self.f1)*2*self.p0_err)**2+(2*self.p1*self.f0*self.f1_err)**2+(2*self.f0*self.f1*self.p1_err)**2)
			self.paras.extend(['f2','p2'])
		else:
			self.f2=0
			self.p2=0
		#
		if 'F3' in paras_key:
			if 'f2' not in self.paras:
				raise Exception('The parameter F3 is in parfile without F2.')
			i=paras['F3']
			self.f3=np.float64(i[0])
			self.f3_err=np.float64(i[1])
			self.p3=-self.f3*self.p0**2-4*self.f2*self.p0*self.p1-2*self.f1*self.p1**2-2*self.f1*self.p0*self.p2
			self.p3_err=np.sqrt((self.f3_err*self.p0**2)**2+((self.f3*self.p0+2*self.f2*self.p1+self.self.f1*self.p2)*2*self.p0_err)**2+((self.f2*self.p0+self.f1*self.p1)*4*self.p1_err)**2+(2*self.p0*self.p1*self.f2_err)**2+((self.p1**2+self.p0*self.p2)*2*self.f1_err)**2+(2*self.f1*self.p0*self.p2_err)**2)
			self.paras.extend(['f3','p3'])
		elif 'P3' in paras_key:
			if 'f2' not in self.paras:
				raise Exception('The parameter P3 is in parfile without P2.')
			i=paras['P3']
			self.p3=np.float64(i[0])
			self.p3_err=np.float64(i[1])
			self.f3=-self.p3*self.f0**2-4*self.p2*self.f0*self.f1-2*self.p1*self.f1**2-2*self.p1*self.f0*self.f2
			self.f3_err=np.sqrt((self.p3_err*self.f0**2)**2+((self.p3*self.f0+2*self.p2*self.f1+self.self.p1*self.f2)*2*self.f0_err)**2+((self.p2*self.f0+self.p1*self.f1)*4*self.f1_err)**2+(2*self.f0*self.f1*self.p2_err)**2+((self.f1**2+self.f0*self.f2)*2*self.p1_err)**2+(2*self.p1*self.f0*self.f2_err)**2)
			self.paras.extend(['f3','p3'])
		else:
			self.f3=0
			self.p3=0
		#
		self.deal_para('pepoch',paras,paras_key,exce='No PEPOCH in par file.')
		#
		if 'DM' in paras_key:
			i=paras['DM'].reshape(-1,2).T
			self.dm=np.float64(i[0]).reshape(-1)
			self.dm_err=np.float64(i[1]).reshape(-1)
			self.paras.append('dm')
			self.dmmodel=0
			if 'DMX' in paras_key:
				i=paras['DM'].reshape(-1,4).T
				self.dmxr1=np.array(list(map(lambda x:te.time(np.float64(x),0,scale=self.units),np.reshape(i[0],-1))))
				self.dmxr2=np.array(list(map(lambda x:te.time(np.float64(x),0,scale=self.units),np.reshape(i[1],-1))))
				self.dmx=np.float64(i[2]).reshape(-1)
				self.dmx_err=np.float64(i[3]).reshape(-1)
				self.paras.extend(['dmx','dmxr1','dmxr2'])
			else:
				self.dmx=np.array([])
				self.dmxr1=np.array([])
				self.dmxr2=np.array([])
		elif 'DMMODEL' in paras_key:
			if 'DMOFF' not in paras_key:
				raise Exception('The DMMODEL para is in parfile without DMOFF.')
			i=paras['DMMODEL']
			self.dmmodel=np.float64(i[0])
			self.dmmodel_err=np.float64(i[1])
			self.paras.append('dmmodel')
			i=paras['DMOFF'].reshape(-1,3).T
			self.dmoffs_mjd=np.array(list(map(lambda x:te.time(np.float64(x),0,scale=self.units),np.reshape(i[0],-1))))
			self.dmoffs=np.float64(i[1]).reshape(-1)
			self.dmoffs_err=np.float64(i[2]).reshape(-1)
			self.paras.extend(['dmmodel','dmoffs','dmoffs_mjd'])
		else:
			raise Exception('No pulsar DM in par file.')
		#
		self.deal_para('dm_s1yr',paras,paras_key,err_case=['DM_C1YR' not in self.paras],err_exc=['The parameter DM_S1YR is in parfile without DM_C1YR.'])
		self.deal_para('dm_c1yr',paras,paras_key,err_case=['DM_S1YR' not in self.paras],err_exc=['The parameter DM_C1YR is in parfile without DM_S1YR.'])
		self.deal_para('fddc',paras,paras_key,err_case=['FDDI' not in self.paras],err_exc=['The parameter FDDC is in parfile without FDDI.'])
		self.deal_para('fddi',paras,paras_key,err_case=['FDDC' not in self.paras],err_exc=['The parameter FDDI is in parfile without FDDC.'])
		self.deal_para('fd',paras,paras_key,err_case=['FDDI' not in self.paras],err_exc=['The parameter FD is in parfile without FDDI.'])
		self.deal_para('cm',paras,paras_key,err_case=['CMIDX' not in self.paras],err_exc=['The parameter CM is in parfile without CMIDX.'])
		self.deal_para('cmidx',paras,paras_key,err_case=['CM' not in self.paras],err_exc=['The parameter CMIDX is in parfile without CM.'])
		self.deal_para('dmepoch',paras,paras_key,exce='Warning: No DMEPOCH in the parfile, using PEPOCH instead.',value=self.posepoch)
		self.deal_para('pmra',paras,paras_key,err_case=['PMDEC' not in self.paras],err_exc=['The parameter PMRA is in parfile without PMDEC.'])
		self.deal_para('pmdec',paras,paras_key,err_case=['PMRA' not in self.paras],err_exc=['The parameter PMDEC is in parfile without PMRA.'])
		self.deal_para('pmra2',paras,paras_key,err_case=['pmra' not in self.paras,'PMDEC2' not in self.paras],err_exc=['The parameter PM2 is in parfile without PM.','The parameter PMRA2 is in parfile without PMDEC2.'])
		self.deal_para('pmdec2',paras,paras_key,err_case=['pmdec' not in self.paras,'PMRA2' not in self.paras],err_exc=['The parameter PM2 is in parfile without PM.','The parameter PMDEC2 is in parfile without PMRA2.'])
		self.deal_para('pmrv',paras,paras_key)
		self.deal_para('px',paras,paras_key)
		self.deal_para('posepoch',paras,paras_key,exce='Warning: No POSEPOCH in the parfile, using PEPOCH instead.',value=self.posepoch)
		self.deal_para('binary',paras,paras_key)
		if self.binary:
			if self.psr.binary=='BT2P' or self.psr.binary=='BT1P': self.psr.binary='T2'
			elif self.psr.binary='T2-PTA': self.psr.binary='T2_PTA'
			for i in eval('paras_'+self.binary)['necessary']:
				self.deal_para(i,paras,paras_key, exce='No '+i.upper()+' parameter for '+self.binary+'model')
			for i in eval('paras_'+self.binary)['optional']:
				self.deal_para(i,paras,paras_key)
			if not (self.pb or self.fb):
				raise Exception('The binary orbit period is not given in par file.')
		self.deal_para('rm',paras,paras_key)
		self.deal_para('dshk',paras,paras_key)
	#
	def deal_para(self,paraname,paras,paras_key,exce=False,value=0,err_case=[],err_exc=[]):
		paraname0=paraname.upper()
		if paraname in alias_keys:
			for i in alias[paraname]:
				if i in paras_key: paraname0=i
		if paraname0 in paras_key:
			for i,k in zip(err_case,err_exc): 
				if i: raise Exception(k)
			i=paras[paraname0]
			if paraname in paras_float:
				i=np.array(np.float64(i))
			elif paraname in paras_float_array:
				i=np.float64(i).reshape(-1,2).T
			elif paraname in paras_time:
				i=np.array(te.time(np.float64(i[0]),0,scale=self.units))
			elif paraname in paras_time_array:
				i=np.array([list(map(lambda x:te.time(np.float64(x),0,scale=self.units),np.reshape(i[0],-1)))])
			self.__setattr__(paraname,np.float64(i[0]))
			if paraname in para_with_err:
				self.__setattr__(paraname+'_err',np.float64(i[0]))
			self.paras.append(paraname)
		elif exce[:7]=='Warning': 
			self.__setattr__(paraname,value)
			self.paras.append(paraname)
			print(exce)
		elif exce: raise Exception(exce)
		else: 
			if paraname in paras_float:
				self.__setattr__(paraname,0)
			elif paraname in paras_float_array:
				self.__setattr__(paraname,np.array([]))
		#
#
all_paras={'p0', 'pmra', 'pmdec', 'pmrv', 'f0', 'p2', 'p1', 'pmra2', 'pmdec2', 'f1', 'p3', 'f2', 'f3', 'dm', 'cm', 'dmmodel', 'dmoffs', 'dmx', 'dm_s1ry', 'dm_c1yr', 'fddc', 'fd', 'raj', 'decj', 'px' ,'cmidx', 'fddi', 'rm', 'pepoch', 'posepoch', 'dmepoch', 'dmoffs_mjd', 'dmxr1', 'dmxr2', 'name', 'dshk', 'binary','t0', 'pb', 'ecc', 'pbdot', 'a1dot', 'a1', 'omdot', 'om', 'gamma','bpjep','bpjph','bpja1','bpjec','bpjom','bpjpb', 'fb', 'tasc', 'eps1', 'eps2', 'sini', 'm2', 'eps1dot', 'eps2dot', 'orbifunc', 'xpbdot', 'edot', 'kom', 'kin', 'mtot', 'a2dot', 'e2dot', 'orbpx', 'dr', 'dtheta', 'a0', 'b0', 'om2dot', 'xomdot','afac','daop', 'pb2dot'}
# with or without error
para_with_err={'p0', 'pmra', 'pmdec', 'pmrv', 'f0', 'p2', 'p1', 'pmra2', 'pmdec2', 'f1', 'p3', 'f2', 'f3', 'dm', 'cm', 'dmmodel', 'dmoffs', 'dmx', 'dm_s1ry', 'dm_c1yr', 'fddc', 'fd', 'raj', 'decj', 'px', 'cmidx', 'fddi', 'rm', 'dshk', 't0', 'pb', 'ecc', 'pbdot', 'a1dot', 'a1', 'omdot', 'om', 'gamma', 'fb', 'bpjph','bpja1','bpjec','bpjom','bpjpb', 'tasc', 'eps1', 'eps2', 'sini', 'm2', 'eps1dot', 'eps2dot', 'orbifuncV','h3', 'h4','nharm','stig', 'xpbdot', 'edot', 'kom', 'kin', 'mtot', 'a2dot', 'e2dot', 'orbpx', 'dr', 'dtheta', 'a0', 'b0', 'om2dot'}
para_without_err={'pepoch','posepoch','dmepoch','dmoffs_mjd','dmxr1','dmxr2','name', 'binary', 'orbifunc', 'bpjep', 'orbifunc', 'orbifuncT', 'xomdot','afac','daop', 'pb2dot'}
# tdb to tcb
paras_p1={'p0','dshk', 'pb', 'bpjpb', 'm2', 'h3', 'mtot', 'a0', 'b0'}
paras_m1={'pmra','pmdec','pmrv','f0','p2', 'omdot', 'a1dot', 'eps1dot', 'eps2dot', 'edot', 'xomdot', 'pb2dot'}
paras_m2={'pmra2','pmdec2','f1','p3', 'a2dot', 'e2dot', 'om2dot'}
paras_m3={'f2'}
paras_m4={'f3'}
paras_eph={'pepoch','posepoch','dmepoch','dmoffs_mjd','dmxr1','dmxr2', 't0', 'tasc', 'orbifuncT'}
paras_mix={'fb':-1,'dm':-1}
uncertain_pm={'dm','cm','dmmodel','dmoffs','dmx','dm_s1ry','dm_c1yr','fddc','fd','rm'}
# para type
paras_float={'p0', 'pmra', 'pmdec', 'pmrv', 'f0', 'p2', 'p1', 'pmra2', 'pmdec2', 'f1', 'p3', 'f2', 'f3', 'dmmodel', 'dm_s1ry', 'dm_c1yr', 'fddc', 'raj', 'decj', 'px' ,'cmidx', 'fddi', 'rm', 'dshk', 't0', 'pb', 'ecc', 'pbdot', 'a1dot', 'a1', 'omdot', 'om', 'gamma', 'tasc', 'eps1', 'eps2', 'sini', 'm2', 'eps1dot', 'eps2dot', 'orbifunc','h3', 'h4','nharm','stig', 'xpbdot', 'edot', 'kom', 'kin', 'mtot', 'a2dot', 'e2dot', 'orbpx', 'dr', 'dtheta', 'a0', 'b0', 'om2dot', 'xomdot','afac','daop', 'pb2dot'}
paras_float_array={'dm', 'cm', 'dmx', 'dmoffs', 'fd','bpjph','bpja1','bpjec','bpjom','bpjpb', 'fb', 'orbifuncV'}
paras_time={'pepoch', 'posepoch', 'dmepoch'}
paras_time_array={'dmoffs_mjd', 'dmxr1', 'dmxr2', 'bpjep', 'orbifuncT'}
paras_text={'name', 'binary'}
# binary
paras_BT={'necessary':['t0', 'pb', 'ecc', 'a1', 'om'], 'optional':['pbdot', 'a1dot', 'omdot', 'gamma']}
paras_BTJ={'necessary':['t0', 'pb', 'ecc', 'a1', 'om','bpjep','bpjph','bpja1','bpjec','bpjom','bpjpb'], 'optional':['pbdot', 'a1dot', 'omdot', 'gamma']}
paras_BTX={'necessary':['t0', 'fb', 'ecc', 'a1', 'om'], 'optional':['pbdot', 'a1dot', 'omdot', 'gamma']}
paras_ECC1={'necessary':['tasc', 'eps1', 'eps2', 'a1'], 'optional':['pb', 'fb', 'pbdot', 'sini', 'a1dot', 'm2', 'eps1dot', 'eps2dot', 'orbifunc', 'orbifuncT', 'orbifuncV']}
paras_ECC1H={'necessary':['tasc', 'eps1', 'eps2', 'a1', 'pb', 'h3'], 'optional':['pbdot', 'sini', 'a1dot', 'm2', 'eps1dot', 'eps2dot', 'h4','nharm','stig']}
paras_ECC1k={'necessary':['tasc', 'eps1', 'eps2', 'pb', 'a1'], 'optional':['pbdot', 'sini', 'a1dot', 'm2', 'omdot']}
paras_DD={'necessary':['t0', 'pb', 'ecc', 'a1', 'omdot'], 'optional':['pbdot', 'sini', 'om', 'a1dot', 'm2', 'xpbdot', 'edot', 'gamma']}
paras_DDH={'necessary':['t0', 'pb', 'ecc', 'a1', 'omdot', 'h3', 'stig'], 'optional':['pbdot', 'om', 'a1dot', 'xpbdot', 'edot', 'gamma']}
paras_DDK={'necessary':['t0', 'pb', 'ecc', 'a1', 'omdot', 'kom', 'kin'], 'optional':['pbdot', 'a1dot', 'm2', 'xpbdot', 'edot', 'gamma']}
paras_DDS={'necessary':['t0', 'pb', 'ecc', 'a1', 'omdot'], 'optional':['pbdot', 'shapmax', 'a1dot', 'om', 'm2', 'xpbdot', 'edot', 'gamma']}
paras_DDGR={'necessary':['t0', 'pb', 'ecc', 'a1'], 'optional':['pbdot', 'sini', 'om', 'a1dot', 'mtot', 'm2', 'xpbdot', 'edot']}
paras_MSS={'necessary':['t0', 'pb', 'ecc', 'a1','om'], 'optional':['pbdot', 'shapmax', 'a1dot', 'omdot', 'm2', 'sini', 'edot', 'gamma', 'a2dot', 'e2dot', 'orbpx', 'dr', 'dtheta', 'a0', 'b0', 'om2dot']}
paras_TT={'necessary':[],'optional':['t0', 'pb', 'ecc', 'a1', 'om', 'tasc', 'eps1', 'eps2', 'shapmax', 'kom', 'kin', 'h3', 'stig', 'h4','nharm', 'm2', 'mtot', 'xpbdot', 'sini', 'pbdot', 'a1dot', 'omdot', 'gamma', 'bpjep', 'bpjph', 'bpja1', 'bpjec', 'bpjom', 'bpjpb', 'edot', 'dr', 'dth', 'a0', 'b0', 'eps1dot', 'eps2dot', 'xomdot','afac','daop', 'pb2dot']}
#
aliase={'ecc':['e'],'a1dot':['xdot'],'a2dot':['x2dot'],'daop':['d_aop']}
aliase_keys=aliase.keys()

