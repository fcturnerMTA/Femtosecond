# -*- coding: utf-8 -*-
"""
FFT Practice

Try playing around with parameters to understand how it works, then load the 
tide data and see if you can determine the relevant frequencies
"""

from scipy.fft import fft, ifft, fftfreq
import numpy as np
import matplotlib.pyplot as plt

# Measurement Characteristics 
wtime = 5.0     # Measurement Window in s
sfreq = 100.0   # Sampling Frequency in measurements/s
Nfloat = np.ceil(wtime*sfreq)   # Sampling Points
N = Nfloat.astype(int)          # Change to int
sT = 1/sfreq    # Temporal Resolution

# Laser Characteristics
pulse = 0       # Set to 1 for Pulsed laser; 0 for continuous-wave (cw)
cfreq = 4.0     # Carrier Frequency in cycles/s
pw = .5         # Pulse Width in s
pc = 2.5        # Pulse Centre in s relative to start of window

# Arrays 
t = np.linspace(0.0, N*sT, N, endpoint=False)   # Time
A = np.sin(cfreq * 2.0*np.pi*t)*np.exp(-pulse*(t-pc)**2/(2*(pw)**2))     # Amplitude
Af = fft(A)                  # Spectral Amplitude
f = fftfreq(N, sT)[:N//2]    # The Positive Frequencies: Note the half

# Plots
plt.subplot(211)
plt.plot(t, A)
plt.subplot(212)
plt.plot(f, 2.0/N * np.abs(Af[0:N//2]))
#plt.semilogy(f, 2.0/N * np.abs(Af[0:N//2]))

#plt.subplot(313)
#plt.plot(f, 2.0/N * np.angle(Af[0:N//2]))
plt.show()