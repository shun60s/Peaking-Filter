#coding:utf-8

#  
#  plot some filter frequenency characteristic using scipy
#  
#  some filter class
#  1: LPF4 100Hz  order 4
#  2: HPF4 5KHz  order 4
#  3: BPF4_butter 1KHz-2KHz order 4
#  4: BPF2_Q  1.5KHz Q=707 order 2
#  5: LPF1 moving average points number 1024  order 1


import numpy as np
from scipy import signal   # version > 1.2.0
from matplotlib import pyplot as plt

# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.18.4 
#  matplotlib  2.1.1
#  scipy 1.4.1

class HPF4(object):
    # fc cut off frequency
    # N filter order
    # sr sampling rate
    # num  f-points
    def __init__(self,fc=5000, N=4, sr=48000, num=1000 ):
        self.fc= fc
        self.N= N
        self.sr= sr
        self.num= num
        self.b, self.a = signal.iirfilter(self.N,  self.fc ,  btype='highpass', ftype='butter', fs=self.sr )
        self.w, self.h = signal.freqz(self.b, self.a, self.num, fs=self.sr)
        self.hlist= 20 * np.log10(abs(self.h))
        self.title= 'highpass butter ' + str(self.N) + ' fc ' + str(self.fc)
        print ( 'HPF4 b len, a len', len(self.b), len(self.a))
        
    def __call__(self, x_in):
        return signal.lfilter(self.b, self.a, x_in)

class LPF4(object):
    # fc cut off frequency
    # N filter order
    # sr sampling rate
    # num  f-points
    def __init__(self,fc=100, N=4, sr=48000, num=1000):
        self.fc= fc
        self.N= N
        self.sr= sr
        self.num= num
        self.b, self.a = signal.iirfilter(self.N,  self.fc ,  btype='lowpass', ftype='butter', fs=self.sr )
        self.w, self.h = signal.freqz(self.b, self.a, self.num, fs=self.sr )
        self.hlist= 20 * np.log10(abs(self.h))
        self.title= 'lowpass butter ' + str(self.N) + ' fc ' + str(self.fc)
        print ( 'LPF4 b len, a len', len(self.b), len(self.a))
        
    def __call__(self, x_in):
        return signal.lfilter(self.b, self.a, x_in)

class BPF4_butter(object):
    # fc1, fc2  cut off frequency
    # N filter order, N should be even number.
    # sr sampling rate
    # num  f-points
    def __init__(self,fc1=1000, fc2=2000, N=4, sr=48000, num=1000):
        self.fc1= fc1
        self.fc2= fc2
        self.N= N
        self.N2= int(N/2)
        self.sr= sr
        self.num= num
        self.b, self.a = signal.iirfilter(self.N2,  [self.fc1, self.fc2 ],  btype='bandpass', ftype='butter', fs=self.sr )
        self.w, self.h = signal.freqz(self.b, self.a, self.num, fs=self.sr )
        self.hlist= 20 * np.log10((np.maximum(abs(self.h), 1e-10)))
        self.title= 'bandpass butter ' + str(self.N) + ' fc1 ' + str(self.fc1) + ' fc2 ' + str(self.fc2)
        print ( 'BPF4_butter b len, a len', len(self.b), len(self.a))
        
    def __call__(self, x_in):
        return signal.lfilter(self.b, self.a, x_in)

class BPF2_Q(object):
    # fc center frequency
    # Q
    # N filter order, N should be even number.
    # sr sampling rate
    # num  f-points
    def __init__(self,fc=1500, Q=0.7071, N=2, sr=48000, num=1000):
        self.fc= fc
        self.Q= Q
        self.BW= self.fc / self.Q
        self.N= N
        self.N2 = int(N/2)
        self.sr= sr
        self.num= num
        
        z,p,k= signal.buttap(self.N2)  # analog prototype of Nth-order Butterworth filter
        b,a= signal.zpk2tf(z,p,k)     # polynomial transfer function representation from zeros and poles
        b,a= signal.lp2bp(b,a, wo=self.fc, bw=self.BW) # Transform a lowpass filter prototype to a bandpass filter
        self.ws, self.hs= signal.freqs(b, a, worN=np.logspace(1, 4.3, self.num)) # analog
        self.hslist= 20 * np.log10(abs(self.hs)) # analog
        self.titles= 'bandpass (analog) ' + str(self.N) + ' fc ' + str(self.fc) + ' Q ' + str(self.Q)
        
        self.b, self.a = signal.bilinear(b, a, fs=self.sr/( 2.0 * np.pi )) # Transform from the analog s-plane to the digital z-plane 
                                                                           # The reason (fs=fs.sr/!(2.0 * np.pi)!) is unknown.
        self.wz, self.hz = signal.freqz(self.b, self.a, worN=np.logspace(1, 4.3, self.num), fs=self.sr)
        self.hzlist= 20 * np.log10(np.maximum(abs(self.hz), 1e-10))
        self.title= 'bandpass (digital) ' + str(self.N) + ' fc ' + str(self.fc) + ' Q ' + str(self.Q)
        print ( 'BPF2_Q b len, a len', len(self.b), len(self.a))
        
    def __call__(self, x_in):
        return signal.lfilter(self.b, self.a, x_in)


class LPF1(object):
    # MAPN moving average points number
    # sr sampling rate
    # num  f-points
    def __init__(self, MAPN=1024, sr=48000, num=1000):
        self.MAPN= MAPN # moving average points number
        self.N= 1  # lpf order
        self.sr= sr
        self.num= num
        self.fc= 0.443/ np.sqrt(self.MAPN * self.MAPN - 1.0) *  self.sr
        print ('fc ', self.fc)
        self.b, self.a = signal.iirfilter(self.N,  self.fc ,  btype='lowpass', ftype='butter', fs=self.sr )
        self.w, self.h = signal.freqz(self.b, self.a, self.num, fs=self.sr )
        self.hlist= 20 * np.log10(abs(self.h))
        self.title= 'lowpass butter ' + str(self.N) + ' fc ' + str(self.fc)
        print ( 'LPF1 b', self.b, '   alfa ', self.b[0])
        print ( 'LPF1 a', self.a)
        
    def __call__(self, x_in):
        return signal.lfilter(self.b, self.a, x_in)


if __name__ == '__main__':
    #
    hpf= HPF4()
    lpf= LPF4()
    bpf4b= BPF4_butter()
    bpf2q= BPF2_Q()
    lpf1=LPF1()
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    ax.semilogx( hpf.w,   hpf.hlist, label=hpf.title)
    ax.semilogx( lpf.w,   lpf.hlist, label=lpf.title)
    ax.semilogx( bpf4b.w, bpf4b.hlist, label=bpf4b.title)
    ax.semilogx( bpf2q.ws, bpf2q.hslist, label=bpf2q.titles)
    ax.semilogx( bpf2q.wz, bpf2q.hzlist, label=bpf2q.title)
    ax.semilogx( lpf1.w, lpf1.hlist, label=lpf1.title)
    
    ax.set_title('Filter Frequenency Characteristic')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Amplitude [dB]')
    ax.set( xlim=(10,20000), ylim=(-120,10))
    ax.grid(which='both', axis='both')
    plt.legend(bbox_to_anchor=(0.02, 0.05), loc='lower left', borderaxespad=0, fontsize=10)
    plt.show()
    
    
    