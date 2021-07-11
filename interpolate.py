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

print(np.min(np.abs(w-np.roll(w,1))))
print(w)
res= np.min(np.abs(w - np.roll(w, 1)))
D= np.max(w) - np.min(w)
Nfloat= np.ceil(D/res)
N= Nfloat.astype(int) 
dw= D/Nfloat

we= np.linspace(np.min(w), np.max(w), N)
lin= interpolate.interp1d(w, Iw, kind='linear')
Ie= lin(we)

# print(w)
# print(we)
# print(wr)
# print(Iw)
# print(Ie)
# print(Ir)
# print(np.sum(Il))
# print(np.sum(Ie)*dw)
# print(np.sum(Ir)*dwr)


plt.subplot(211)
plt.plot(we,Ie)
plt.subplot(212)
plt.plot(l,Il)

print(np.trapz(Ie, dx=dw))
print(np.trapz(Il, dx=1))
