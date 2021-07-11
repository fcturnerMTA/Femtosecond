import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# l= np.matrix([1.0,2.0,3.0,4.0,5.0,6.0,7.0])
# l= np.asarray(l)
# l= l.flatten()
 
# Il= np.matrix([0.0,0.1,0.3,0.6,1.0,1.5,2.1])
# Il= np.asarray(Il)
# Il= Il.flatten()

l= np.linspace(700, 750, 50)
Il= np.exp(-(l-725)**2/2/5**2)

w= 1e7/l
Iw=1e7*Il/w**2

wres= np.min(np.abs(w - np.roll(w, 1)))
D= np.max(w) - np.min(w)
Nfloat= np.ceil(D/wres)
N= Nfloat.astype(int) 

we= np.linspace(np.min(w), np.max(w), N)
lin= interpolate.interp1d(w, Iw, kind='linear')
Ie= lin(we)


dw= D/Nfloat
dl= np.min(np.abs(l - np.roll(l,1)))

Ite= np.trapz(Ie, dx=dw)
Itl= np.trapz(Il, dx=dl)

plt.subplot(211)
plt.plot(we,Ie)
plt.subplot(212)
plt.plot(l,Il)

print("Approximately", np.round(100/Itl*(Ite - Itl),1), "% Error between Wavenumber and Wavelength Integrals")
