# Dipole-ACF
This little project was one of my works in *ab initio* molecular dynamics (AIMD) simulations, in which a Python (version 2.7) script was composed for calculating the IR spectrum based on the fast Fourier transform (FFT) of the autocorrelation function. 

This script first read the **total dipole moment (Dipole in short)** data, which is generated from the [CP2K/QuickStep] (https://www.cp2k.org/quickstep "CP2K") simulations, and then calculates the time derivative of the **Dipole**, yielding the **dipole prime (D_p in short)**. After computing the autocorrelation of the **D_p**, the **DACF** data array was obtained. By performing the FFT on the **DACF**, the final IR spectrum is produced. And then plotted on the graph panel by using the Matplotlib module.

## Modules required:
- Numpy (version 1.9.1 or above)
- Scipy (version 0.17.0 or above) 
- Matplotlib (version 1.4 or above)

## Short history
  This script (VERSION 3.3) is the improved version based on Dr. Kulig's first version of `ir_total_QC.py`, which was a pure Python script without using Numpy, Scipy, and no visualization of the results.
    
  The author renamed it into IR_DACF_KW.py, where "K" refers to Dr. Kulig, and "W" refers to Dr. Wang.

## The main improvements are:

1. Implementation of the powerful Numpy module, which facilitates the fast calculation of the data array. The Numpy module accelerates the calculations dramatically by converting all data lists into data arrays. 
    Usually, the calculations would complete within 3 seconds.

2. Built a `zero_padding()` function. This function dynamically adds a series of zeros to the end of the Dipole moment array before FFT. The length of the whole data series is the power-of-two (2<sup>n</sup>).
    + *[Note] FFT (Fast Fourier Transform) refers to a way the discrete Fourier Transform (DFT) can be calculated efficiently, by using symmetries in the calculated terms. The symmetry is highest when n is a power of 2, and the transform is, therefore, most efficient for these sizes.*

3. Using built-in `signal.fftconvolve()` function in scipy.signal module for accelerating the auto-correlation function calculation.

4. A `choose_window()` function was taken into consideration for suppressing noise. The window function is imported from scipy.signal module. 

5. Built a `visualization()` function for plotting the results.


## Contribution:
Dr. Huan Wang         (The 3rd and 2nd versions)  
Dr. Waldemar Kulig    (The 1st version)

## E-mail address for contacting the authors:
`huan.wang@mail.huji.ac.il`  or  `wanghuan@iccas.ac.cn (China)`

## Copyright:
The Hebrew University of Jerusalem, Givat Ram, Jerusalem, 91904, Israel.

## PLEASE READ THE FOLLOWING INSTRUCTIONS BEFORE RUNNING SCRIPT
###  The Format for Running This Script:
`python` `IR_total_KW.py` *`INPUT_FILE`* *`DELTA_T`* *`WINDOW`* *`OUTPUT_FILE`*

###  The values need to be input manually when running this script
  1. *`INPUT_FILE`*: The Total_Dipole_Moment.Diople file
      + <sub>*NOTE: do NOT need to re-split the Dipole file)*</sub>
  2. *`DELTA_T`*: The Time_step set in simulation, in unit of fs
  3. *`WINDOW`*: The Name of the Window Function
  4. *`OUTPUT_FILE`*: The Name of the Output File
      + <sub>*NOTE: do NOT need to type '>' sign!*</sub>
#### ########################  Let's Try It! ###################### ####
