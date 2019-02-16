#coding:utf-8

# A class of iir peaking filter

# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.14.0 
#  scipy 1.0.0
#  matplotlib  2.1.1


import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


class Class_IIR_Peaking1(object):
    def __init__(self, fpeak=1000, gain=2.0, Q=1.0, sampling_rate=48000):
        # design iir peaking filter
        # initalize
        self.fpeak= fpeak # peak frequecny by unit is [Hz]
        self.sr= sampling_rate # sampling frequecny by unit is [Hz]
        self.gain= gain # amplification factor (magnification).   This must be > 0.0
        self.Q= Q # Q factor
        self.b, self.a= self.set_peaking( self.fpeak, self.gain, self.Q)
        print ('self.b,self.a', self.b, self.a)

    def filtering(self, xin):
    	# process of filtering
    	# input xin
    	# output filtered xin
        return signal.lfilter(self.b, self.a, xin)

    def f_show(self, worN=1024):
        # show frequency response
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

    def set_peaking(self, fpeak, gain0, Q0):
        
        omega= (fpeak / self.sr) * np.pi * 2.0
        sn= np.sin(omega)
        cs= np.cos(omega)
        alpha = sn / (2.0 * Q0)
        
        A=np.sqrt(gain0) 
        b=np.zeros(3) # umerator(bunsi)
        a=np.zeros(3) # denominator(bunbo)
        
        # if flat (gain is 1.0)
        if gain0 == 1.0:
        	a[0]=1.0
        	b[0]=1.0
        	return b,a
        
        a[0]= 1.0 + alpha / A
        a[1]= -2.0 * cs
        a[2]= 1.0 - alpha / A
        
        b[0]= 1.0 + alpha * A
        b[1]= -2.0 * cs
        b[2]= 1.0 - alpha * A
        
        b /= a[0]
        a /= a[0]
        
        return b, a


if __name__ == '__main__':
    
    # Boost sample 
    iir_pk1=Class_IIR_Peaking1(fpeak=500, gain=2.0, Q=1.5)
    # draw frequency response
    iir_pk1.f_show()
    
    
    # Drop sample 
    iir_pk1=Class_IIR_Peaking1(fpeak=2000, gain=0.5, Q=0.7)
    # draw frequency response
    iir_pk1.f_show()