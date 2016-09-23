# Dipole-ACF
This little project was one of my works in *ab initio* molecular dynamics (AIMD) simulations, in which a Python (version 2.7) script was composed for calculating the IR spectrum based on the fast Fourier transform (FFT) of the autocorrelation function. 

This script first read the **total dipole moment (Dipole in short)** data, whcih generated from the [CP2K/QuickStep] (https://www.cp2k.org/quickstep "CP2K") simulations, and then calculate the time derivative of the **Dipole**, yielding the **dipole prime (D_p in short)**. After computing the autocorrelation of the **D_p**, the **DACF** data array obtained. By performing the FFT on the **DACF**, the final IR spectrum produced. And then plotted on the graph panel by using the Matplotlib module.

Modules required:
- Numpy (version 1.8.2 or above)
- Scipy (version 0.17.0 or above) 
- Matplotlib (version 1.4 or above)