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
pulse = 1       # Set to 1 for Pulsed laser; 0 for continuous-wave (cw)
cfreq = 10.0    # Carrier Frequency in cycles/s
chirp = 4.0     # Linear frequency shift in (cycles/s)/s
pw = 0.2        # Pulse Width in s
pc = 0.0        # Pulse Centre in s relative to t = 0

# Arrays 
t = np.linspace(-N/2*sT, N/2*sT, N, endpoint=False)   # Time
f = fftfreq(N, sT)[:N//2]   # The Positive Frequencies: Note the half
g = fftfreq(N, sT)          # Optional: All frequencies; use this to visualize what's happening with the phase

# Pulse characteristics and transforms
E = np.sin((cfreq - chirp*(t-pc)) * 2.0*np.pi*(t-pc))*np.exp(-pulse*(t-pc)**2/(2*(pw)**2))  # Amplitude of the electric field in time
Ef = fft(E)                             # Transform to frequency space
Spec = np.abs(Ef)                       # Spectral Amplitude of Electric Field
phase = np.unwrap(2*np.angle(Ef))/2     # Spectral Phase, removing 2*pi jumps

# Controls
tfl = 0  # Set to 1 for Transform limited Pulse; 0 for chosen phase
cmp = 0  # Set to 1 to apply new phase to old phase, or 0 for a clean slate

#Apply Phase: Various parameters
arb = 0.0*f     # Arbitrary Phase Compensation: leave this in for later
GDD = 0.0       # Group Delay Dispersion in s/cycle (Linear Phase)
GVD = -5.7E-1   # Group Velocity Dispersion in (s/cycle)^2 (Quadratic Phase)
TOD = 0.0       # Third Order Dispersion in (s/cycle)^3 (Cubic Phase)

#Apply Phase: The Actual Compensation with controls
phi = (1-tfl)*(GDD*(f-cfreq) + 1/2*GVD*(f - cfreq)**2 + 1/6*TOD*(f - cfreq)**3 + arb + cmp*phase[0:N//2])  # Overall dispersion using a Taylor Expansion, plus an arbitrary addition for flexibility
phi = np.concatenate((phi,-np.flip(phi)))  # This now applies phase to negative frequencies too
Ec = Spec*np.exp(1j*phi)      # The spectrum with the applied phase
Et = ifft(Ec)                 # Going back to time using an ifft
Et = np.fft.ifftshift(Et)     # Putting the pulse back into the center

# Plots
plt.subplot(321)
plt.plot(t, E)                      # The Pulse in Time
plt.subplot(323)
plt.plot(f, 2.0/N * Spec[0:N//2])   # The Spectrum
plt.subplot(325)
plt.plot(f, phase[0:N//2])          # The Relative Phase Angles
plt.subplot(326)
plt.plot(f, phi[0:N//2])            # Compensated Phase Angles
plt.subplot(324)
plt.plot(f, 2.0/N * Spec[0:N//2])   # Still the Same Spectrum
plt.subplot(322)
plt.plot(t,np.real(Et))             # The Pulse with Applied Phase

plt.show()