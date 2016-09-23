#!/usr/bin/env python

''' This script (VERSION 3.3) is the improved version based on Dr. Kulig's first
version of ir_total_QC.py, which is a pure Python script, without using Numpy, 
Scipy and no visulization of the results.
    
    The author renamed it into IR_DACF_KW.py, where "K" refers to Dr. Kulig, 
and "W" refers to Dr. Wang.

The main improvement are:

(1) Implementation of the powerful Numpy module, which facilitating the scientific
calculation in data array.
    By verticalizing all data lists to data array, the Numpy module dramatically 
accelerated the calculations.

(2) Built a "zero_padding" function. This function dynamically add a series of zeros
to the end of the Dipole moment array before FFT. The length of the zero-series is 
the power-of-two (2^n).
    *[Note] FFT (Fast Fourier Transform) refers to a way the discrete Fourier 
    Transform (DFT) can be calculated efficiently, by using symmetries in the 
    calculated terms.The symmetry is highest when n is a power of 2, and the 
    transform is therefore most efficient for these sizes.

(3) Using built-in fftconvolve function in scipy.signal module for fast calculating 
the auto-correlation function.

(4) Window Function was taken into consideration for suppressing noise. The window 
function is imported from Scipy module. 

(5) Built a Visualization Function for plotting the results.

Contribution:
Dr. Huan Wang         (The 3rd and 2nd version)
Dr. Waldemar Kulig    (The 1st version)

E-mail address for contacting the authors:

huan.wang@mail.huji.ac.il  or  wanghuan@iccas.ac.cn (China)

Copyright:
The Hebrew University of Jerusalem, Givat Ram, Jerusalem, 91904, Israel.
'''

from __future__ import division
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import time


##### PLEASE READ THE FOLLOWING INSTRUCTIONs BEFORE RUNNING SCRIPT #####
####                                                                ####
####  The Format for Running This Script:                           ####
####  python IR_total_KW.py INPUT_FILE DELTA_T WINDOW OUTPUT_FILE   ####
####                                                                ####
####  The values need to input manually when runing this script     ####
####                                                                #### 
####  (1) INPUT_FILE_NAME: The Total_Dipole_Moment_*.Diople file    ####
####           (NOTE: do NOT need to re-split the Dipole file)      ####
####                                                                #### 
####  (2) DELTA_T: The Time_step set in simulation, in unit of fs   ####
####                                                                ####
####  (3) WINDOW: The Name of the Window Function                   ####
####                                                                ####  
####  (4) OUTPUT_FILE_NAME: The Name of the Output File.            ####
####           (NOTE: do NOT need to type > sign!)                  ####
####                                                                ####
#############################  Let's Try It! ###########################






#### The values need to input manually when running this script ####
fname = sys.argv[1]                  # The name of the input file
delta_t = float(sys.argv[2])*1.0e-15 # The time step in unit of femtoseconds
window = sys.argv[3]                 # The name of the window function
fout = sys.argv[4]                   # The name of the output file



#### The constants will be used in this script ####
c = 2.9979245899e10     # speed of light in vacuum in [cm/s], from Wikipedia.
kB = 0.6950347          # Boltzman constant in [cm^-1/K], from Wikipedia.
h_bar = 6.283185        # Reduced Planck constant in atomic unit, where h = 2*pi
beta = 1.0/(kB * T)                         
start = time.clock()	# The system clock, for checking the running speed. 



#### Functions will used in this script ####

def read_data(fname):
    with open(fname, "r") as fo:
        dipole = np.genfromtxt(fo, dtype=np.float64,
                               delimiter=None, usecols=(1,2,3))
    return dipole




def calc_derivative(data, delta_t):
    dy = np.zeros(np.shape(data))
    for i in xrange(3):
        dy[:,i] = np.gradient(data[:,i])
    print "dy = ", dy
    dy = dy[~(np.absolute(dy) > 0.1).any(1),:]
    return np.divide(dy, delta_t)




def zero_padding(sample_data):
    '''
      A series of Zeros will be padded to the end of the dipole moment array 
    (before FFT performed), in order to obtain a array with the length which
    is the "next power of two" of numbers.
    #### Next power of two is calculated as: 2**np.ceil(log2(x))
    #### or Nfft = 2**int(math.log(len(data_array)*2-1, 2))
    '''
    N = 2**int(math.log(len(sample_data)*2-1, 2))
    return N




def calc_ACF(array):
'''
    This function deals with the auto-correlation function (ACF) of the total
    dipole moment derivatives.

    With the Wiener-Khintchine theorem, the autocorrelation function is
    http://en.wikipedia.org/wiki/Wiener%E2%80%93Khinchin_theorem

####
####  http://stackoverflow.com/questions/4503325/autocorrelation-of-a-multidimensional-array-in-numpy
####
####  for fast convolution 
####  http://sebug.net/paper/books/scipydoc/frequency_process.html#id5
'''
    # normalization
    yunbiased = array - np.mean(array, axis=0)
    ynorm = np.sum(np.power(yunbiased,2), axis=0)
#    print "the average value of input data array", ynorm
    autocor = np.zeros(np.shape(array))

    for i in xrange(3):
        autocor[:,i] = signal.fftconvolve(array[:,i],
                                          array[:,i][::-1],
                                          mode='full')[len(array)-1:]/ynorm[i]
    print "shape of the result3 from signal.FFTcorrelate()", np.shape(autocor)
    return autocor




def choose_window(data, kind='string'):
    if kind == 'Gaussian':
        sigma = 2 * math.sqrt(2 * math.log(2))
        window = signal.gaussian(len(data), std=4000.0/sigma, sym=False)
    elif kind == 'BH':
        window = signal.blackmanharris(len(data), sym=False)
    elif kind == 'Hamming':
        window = signal.hamming(len(data), sym=False)
    elif kind == 'Hann':
        window = signal.hann(len(data), sym=False)
    return window




def calc_FFT(data, window):
    '''
    This function is for calculating the "intensity" of the ACF at each 
    frequency by using the discrete fast Fourier transform.
    
####
#### http://stackoverflow.com/questions/20165193/fft-normalization
####
    '''
    window = choose_window(data, kind=window)
    WE = sum(window) / len(data)
    wf = window / WE
    # convolve the window function. 
    sig = data * wf[None,:].T

    # A series of number of zeros will be padded to the end of the DACF \
    # array before FFT.
    N = zero_padding(sig)
	
    yfft = np.fft.fft(sig, N, axis=0) / len(sig)
# without window function
#    yfft = np.fft.fft(data, n=int(N_fft), axis=0) / len(data)
    return np.square(np.absolute(yfft))




######## Save The Results to A TEXT File ########
def save_results(fout, wavenumber, intensity):
    title = ("Wavenumber", "IR Intensity", "cm^-1", "a.u.")
    with open(fout, "w") as fw:
        np.savetxt(fout, np.c_[wavenumber[0:5000], intensity[0:5000]],
                   fmt="%10.5f %15.5e",
                   header="{0:>10}{1:>16}\n{2:^11}{3:^20}".format(*title),
                   comments='')




######## Plot The Spectrum by Using Matplotlib module ########
def visualization(D_p, DACF, wavenumber, intensity):
    plt.subplot(3,1,1)
    L1 = np.arange(len(D_p))
    plt.plot(L1, D_p[:,0], color='red', linewidth=1.5)
    plt.plot(L1, D_p[:,1], color='green', linewidth=1.5)
    plt.plot(L1, D_p[:,2], color='blue', linewidth=1.5)
    plt.axis([0, len(D_p), 1.1*np.min(D_p), 1.1*np.max(D_p)], fontsize=15)
    plt.xlabel("Data Points", fontsize=15)
    plt.ylabel("Derivative of Dipole (a.u.)", fontsize=15)

    plt.subplot(3,1,2)
    L2 = np.arange(len(DACF))
    plt.plot(L2, DACF[:,0], color='red', linewidth=1.5)
    plt.plot(L2, DACF[:,1], color='green', linewidth=1.5)
    plt.plot(L2, DACF[:,2], color='blue', linewidth=1.5)
    plt.axis([0, len(DACF), 1.1*np.min(DACF), 1.1*np.max(DACF)], fontsize=15)
    plt.xlabel("Data Points", fontsize=15)
    plt.ylabel("DACF (a.u.)", fontsize=15)

    plt.subplot(3,1,3)
    plt.plot(wavenumber, intensity, color='black', linewidth=1.5)
    plt.axis([0, 4000,
             -1.1*np.min(intensity), 1.1*np.max(intensity)],
             fontsize=15)
    plt.xlabel("Wavenumber (cm$^{-1}$)", fontsize=15)
    plt.ylabel("Intensity (a.u.)", fontsize=15)
    plt.subplots_adjust(hspace = 0.5)
    plt.show()




######## The main program ########
if __name__ == '__main__':
    dipole = read_data(fname)
    print "dipole \n", dipole, np.shape(dipole)
    
    D_p = calc_derivative(dipole, delta_t)
    DACF = calc_ACF(D_p)
    yfft = calc_FFT(DACF)
    print "SHAPE OF YFFT = ", np.shape(yfft)

    wavenumber = np.fft.fftfreq(len(yfft), delta_t*c)[0:int(len(yfft)/2)]
    intensity = np.sum(yfft, axis=1)[0:int(len(yfft)/2)]

#### Normalized the intensity
#    intensity = intensity/max(intensity)
    save_results(fout, wavenumber, intensity)
    print "Work Completed! Used Time = ", time.clock() - start
    visualization(D_p, DACF, wavenumber, intensity)
