#coding:utf-8

#  test filter to wav  using scipy


import sys
import argparse
import numpy as np
from scipy.io.wavfile import read as wavread
from scipy.io.wavfile import write as wavwrite

from filter_class1 import *

# Check version
#  Python 3.6.4 on win32 (Windows 10)
#  numpy 1.18.4 
#  scipy 1.4.1


def read_wav( file_path ):
    try:
        sr, w = wavread( file_path)
    except:
        print ('error: wavread ', file_path)
        sys.exit()
    else:
        w= w / (2 ** 15)
        if w.ndim == 2:  # if stereo, convert to mono
            w= np.average(w, axis=1)
        print ('sampling rate ', sr)
        print ('size', w.shape[0])
    return w, sr

def save_wav( file_path, data, sr=48000):
    amplitude = np.iinfo(np.int16).max
    try:
        wavwrite( file_path , sr, np.array( amplitude * data , dtype=np.int16))
    except:
        print ('error: wavwrite ', file_path)
        sys.exit()
    print ('wrote ', file_path)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='test filter to wav using scipy')
    parser.add_argument('--input_wav', '-i', default='wav/MIX100-1500-5000-10dB_ST_10sec.wav', help='wav file name(16bit)')
    args = parser.parse_args()
    
    file_path_in= args.input_wav
    
    w,sr= read_wav( file_path_in )
    
    hpf= HPF4(sr=sr)  # sr sampling rate
    lpf= LPF4(sr=sr)
    bpf4b= BPF4_butter(sr=sr)
    bpf2q= BPF2_Q(sr=sr)
    
    y= hpf(w)
    save_wav('wav/hpf4_out.wav', y, sr=sr)
    
    y= lpf(w)
    save_wav('wav/lpf4_out.wav', y, sr=sr)
    
    y= bpf4b(w)
    save_wav('wav/bpf4b_out.wav', y, sr=sr)
    
    y= bpf2q(w)
    save_wav('wav/bpf2q_out.wav', y, sr=sr)
    
    
    # moving average
    N=1024  # MAPN moving average points number
    lpf1=LPF1(MAPN=N,sr=sr)
    y= lpf1(np.abs(w))
    save_wav('wav/moving_average_out.wav', y, sr=sr)
    yma=y[:N * int(len(w)/N)][::N] # get every N point data
    print (yma.shape)
    
    # N points mean
    yn=np.abs(w[:N * int(len(w)/N)]).reshape(-1,N).mean(axis=1)
    print ( yn.shape)
    save_wav('wav/N_point_average_out.wav', yn, sr=sr)
    
    # comparison
    yma_vs_yn= np.stack([yma, yn], 1)
    save_wav('wav/moving_average_vs_N_point_average_out.wav', yma_vs_yn, sr=sr)
   
