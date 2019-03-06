#coding:utf-8

#
# A class of IIR Band Pass Filter Bank
#

import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal

# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.14.0 
#  matplotlib  2.1.1
#  scipy 1.0.0


class Class_BPF_bank(object):
	def __init__(self, fbase=100.0, fstep=10.0, fband=5, gain=1.0, Q=10.0, sampling_rate=48000):
		# initalize
		# number of BPF: fband
		# center frequency of 1st BPF: fbase [Hz] 
		# frequency step: fstep [Hz]
		self.fband= fband # number of filter bank
		self.fc_list= np.linspace(fbase,(fband-1) * fstep + fbase, fband)  # center frequency of Band Pass Filter by unit is [Hz]
		self.gain_list= np.ones(fband) * gain # magnification
		# check Q
		if Q <= 0.0:
			print ('error: Q must be > 0')
			sys.exit()
		else:
			self.Q_list= np.ones(fband) * Q   # Q factor
		
		self.sr= sampling_rate
		
		self.a=np.zeros((fband,3))
		self.b=np.zeros((fband,3))
		for i in range( self.fband) :
			a, b = self.bpf1(self.fc_list[i], self.gain_list[i], self.Q_list[i])
			self.a[i]=a
			self.b[i]=b
		
	def bpf1(self,fc,gain,Q):
		# primary digital filter
		a= np.zeros(3)
		b= np.zeros(3)
		
		wc= 2.0 * np.pi * fc / self.sr
		g0= 2.0 * np.tan( wc/2.0)
		
		a[0]=   4.0 +  2.0 * g0 / Q +  g0 * g0
		a[1]=  -8.0 + 2.0 * g0 * g0
		a[2]=   4.0 -  2.0 * g0 / Q +  g0 * g0
		
		b[0]=   2.0 * gain *  g0 / Q
		b[2]=  -2.0 * gain *  g0 / Q
		
		b /= a[0]
		a /= a[0]
		
		return  a,b
		
	def filtering(self, xin):
		# filtering process, using scipy
		yout=np.zeros( (self.fband, len(xin)) )
		for i in range (self.fband):
			zi= signal.lfilter(self.b[i], self.a[i], np.zeros(len(self.a[i])-1) ) 
			xout1, zf = signal.lfilter(self.b[i], self.a[i], xin, zi=zi)
			yout[i]=np.copy(xout1)
			
		return yout # output yout.shape( fband, len(xin) )
		
	def f_show(self, worN=1024*8):
		# show frequency response, using scipy
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		plt.title('BPF bank frequency response')
		plt.ylabel('Amplitude [dB]', color='b')
		plt.xlabel('Frequency [Hz]')
		
		for i in range (self.fband):
			wlist, fres = signal.freqz(self.b[i], self.a[i], worN=worN)
			flist = wlist / ((2.0 * np.pi) / self.sr)
			plt.semilogx(flist, 20 * np.log10(abs(fres)), 'b')  # plt.plot(flist, 20 * np.log10(abs(fres)), 'b')
			
			"""
			ax2 = ax1.twinx()
			angles = np.unwrap(np.angle(fres))
			angles = angles / ((2.0 * np.pi) / 360.0)
			plt.semilogx(flist, angles, 'y')  # plt.plot(flist, angles, 'y')
			plt.ylabel('Angle(deg)', color='y')
			"""
		plt.grid()
		plt.axis('tight')
		plt.show()

if __name__ == '__main__':
	
	# instance
	bpf_bank=Class_BPF_bank(fbase=77.0, fstep=10.0, fband=7, Q=10.0)
	
	# draw frequecny response, using scipy
	bpf_bank.f_show()
	
#This file uses TAB

