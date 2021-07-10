import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

l= np.matrix([1.0,2.0,3.0,4.0,5.0,6.0,7.0])
l= np.asarray(l)
l= l.flatten()
  
Il= np.matrix([0.0,0.1,0.3,0.6,1.0,1.5,2.1])
Il= np.asarray(Il)
Il= Il.flatten()

w= 1/l
Iw=Il/w**2

res= 0.01
D= np.max(w) - np.min(w)
Nfloat= np.ceil(D/res)
N= Nfloat.astype(int) 
dw= D/Nfloat

we= np.linspace(np.min(w), np.max(w), N)
lin= interpolate.interp1d(w, Iw, kind='linear')
Ie= lin(we)

dres= 0.05
Nfloat= np.ceil(D/dres)
Nr= Nfloat.astype(int)
dwr= D/Nfloat

wr= np.linspace(np.min(we), np.max(we), Nr)
cub= interpolate.interp1d(we, Ie, kind='cubic')
Ir= cub(wr)

print(w)
print(we)
print(wr)
print(Iw)
print(Ie)
print(Ir)

plt.subplot(211)
plt.plot(we,Ie)
plt.subplot(212)
plt.plot(wr,Ir)
