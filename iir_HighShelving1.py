#coding:utf-8

# A class of iir High Shelving filter

# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.14.0 
#  scipy 1.0.0
#  matplotlib  2.1.1


import matplotlib.pyplot as plt
import numpy as np
from scipy import signal


class Class_IIR_highShelving1(object):
    def __init__(self, fc=2500, gain=2.0, slope=1.0, sampling_rate=48000):
        # design iir high Shelving filter
        # initalize
        self.fc= fc # midpoint frequecny by unit is [Hz]
        self.sr= sampling_rate # sampling frequecny by unit is [Hz]
        self.gain= gain # amplification factor (magnification).   This must be > 0.0
        self.gain= np.sqrt(self.gain)
        self.slope= slope # shelf slope (S=1 for steepest slope)
        self.b, self.a= self.set_highshelving()
        #print ('self.b,self.a', self.b, self.a)

    def filtering(self, xin):
    	# process filtering, using scipy
    	# input xin
    	# output filtered xin
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

    def set_highshelving(self,):
        
        omega= (self.fc / self.sr) * np.pi * 2.0
        sn= np.sin(omega)
        cs= np.cos(omega)
        alpha = sn / 2.0 * np.sqrt((self.gain + 1.0/self.gain) * (1.0/self.slope - 1.0) + 2.0)
        
        b=np.zeros(3) # umerator(bunsi)
        a=np.zeros(3) # denominator(bunbo)
        
        a[0]= ((self.gain + 1.0 ) - (self.gain - 1.0) * cs + 2.0 * np.sqrt(self.gain) * alpha)
        a[1]= 2.0 * (( self.gain - 1.0) - ( self.gain + 1.0) * cs)
        a[2]= ((self.gain + 1.0) - (self.gain - 1.0 ) * cs - 2.0 * np.sqrt( self.gain) * alpha )
        
        b[0]= self.gain * ((self.gain +1.0 ) + (self.gain - 1.0) * cs + 2.0 * np.sqrt(self.gain) * alpha )
        b[1]= -2.0 * self.gain * ((self.gain - 1.0) + (self.gain + 1.0 ) * cs)
        b[2]= self.gain * ((self.gain + 1.0) + (self.gain - 1.0) * cs - 2.0 * np.sqrt(self.gain) * alpha)
        
        b /= a[0]
        a /= a[0]
        
        return b, a


if __name__ == '__main__':
    
    # high shelf filter sample 
    iir_LS1=Class_IIR_highShelving1( fc=1250, gain=4.0, slope=0.75, sampling_rate=48000)
    # draw frequency response
    iir_LS1.f_show()
    
    
