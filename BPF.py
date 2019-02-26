#coding:utf-8

#
# A class of iir Band Pass Filter
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


class Class_BPF(object):
	def __init__(self, fc=1000, gain=1.0, Q=1.0, sampling_rate=48000):
		# initalize
		self.fc= fc # center frequency of Band Pass Filter by unit is [Hz]
		self.gain= gain # magnification
		self.Q= Q   # Q factor
		# check Q
		if self.Q <= 0.0:
			print ('error: Q must be > 0')
			sys.exit()
		
		self.sr= sampling_rate
		self.a, self.b = self.bpf1()
		
	def bpf1(self,):
		# primary digital filter
		a= np.zeros(3)
		b= np.zeros(3)
		
		wc= 2.0 * np.pi * self.fc / self.sr
		g0= 2.0 * np.tan( wc/2.0)
		
		a[0]=   4.0 +  2.0 * g0 / self.Q +  g0 * g0
		a[1]=  -8.0 + 2.0 * g0 * g0
		a[2]=   4.0 -  2.0 * g0 / self.Q +  g0 * g0
		
		b[0]=   2.0 * self.gain *  g0 / self.Q
		b[2]=  -2.0 * self.gain *  g0 / self.Q
		
		b /= a[0]
		a /= a[0]
		
		return  a,b
		
	def iir2(self,x):
		# filtering process
		# calculate iir filter: x is input, y is output
		# y[0]= b[0] * x[0]  + b[1] * x[-1] + b[2] * x[-1]
		# y[0]= y[0] - a[1] * y[-1] - a[2] * y[-1]
		y= np.zeros(len(x))
		for n in range(len(x)):
			for i in range(len(self.b)):
				if n - i >= 0:
					y[n] += self.b[i] * x[n - i]
			for j in range(1, len(self.a)):
				if n - j >= 0:
					y[n] -= self.a[j] * y[n - j]
		return y
		
	def fone(self, xw):
		# calculate one point of frequecny response
		f= xw / self.sr
		yi= self.b[0] + self.b[1] * np.exp(-2j * np.pi * f) + self.b[2] * np.exp(-2j * np.pi * 2 * f)
		yb= self.a[0] + self.a[1] * np.exp(-2j * np.pi * f) + self.a[2] * np.exp(-2j * np.pi * 2 * f)
		val= yi/yb
		return np.sqrt(val.real ** 2 + val.imag ** 2)

	def H0(self, freq_low=100, freq_high=7500, Band_num=256):
		# get Log scale frequecny response, from freq_low to freq_high, Band_num points
		amp=[]
		freq=[]
		bands= np.zeros(Band_num+1)
		fcl=freq_low * 1.0    # convert to float
		fch=freq_high * 1.0   # convert to float
		delta1=np.power(fch/fcl, 1.0 / (Band_num)) # Log Scale
		bands[0]=fcl
		#print ("i,band = 0", bands[0])
		for i in range(1, Band_num+1):
			bands[i]= bands[i-1] * delta1
			#print ("i,band =", i, bands[i]) 
		for f in bands:
			amp.append(self.fone(f))
		return   np.log10(amp) * 20, bands # = amp value, freq list
		
	def H0_show(self,freq_low=100, freq_high=7500, Band_num=256):
		# draw frequecny response
		plt.xlabel('Hz')
		plt.ylabel('dB')
		plt.title('Band Pass Filter')
		amp, freq=self.H0(freq_low=freq_low, freq_high=freq_high, Band_num=Band_num)
		plt.plot(freq, amp)
		plt.grid()
		plt.show()
		
	def filtering(self, xin):
		# filtering process, using scipy
		return signal.lfilter(self.b, self.a, xin)
		
	def f_show(self, worN=1024):
		# draw frequency response, using scipy
		wlist, fres = signal.freqz(self.b, self.a, worN=worN)
		
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		flist = wlist / ((2.0 * np.pi) / self.sr)
		plt.title('frequency response')
		ax1 = fig.add_subplot(111)
		
		plt.semilogx(flist, 20 * np.log10(abs(fres)), 'b')  # plt.plot(flist, 20 * np.log10(abs(fres)), 'b')
		plt.ylabel('Amplitude [dB]', color='b')
		plt.xlabel('Frequency [Hz]')
		
		ax2 = ax1.twinx()
		angles = np.unwrap(np.angle(fres))
		angles = angles / ((2.0 * np.pi) / 360.0)
		plt.semilogx(flist, angles, 'g')  # plt.plot(flist, angles, 'g')
		plt.ylabel('Angle(deg)', color='g')
		plt.grid()
		plt.axis('tight')
		plt.show()

if __name__ == '__main__':
	
	# instance
	bpf=Class_BPF(fc=1000,  Q=10.0, sampling_rate=48000)
	
	# draw frequecny response, using scipy
	bpf.f_show()
	
	# draw frequecny response
	bpf.H0_show(freq_high=2000)
	
	
#This file uses TAB

