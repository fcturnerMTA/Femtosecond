# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import math
import numpy as np
import matplotlib.pyplot as plt

l= np.matrix([1.0,2.0,3.0,4.0,5.0,6.0,7.0])
#l = np.linspace(1., 7., 7)   
I= np.matrix([0.0,0.1,0.3,0.6,1.0,1.5,2.1])
w= 1/l

res= 0.01
D= np.max(w) - np.min(w)
Nfloat= np.ceil(D/res)
N = Nfloat.astype(int) 
dw = D/Nfloat

we= np.linspace(np.min(w), np.max(w), N)

#l= l.transpose()
#I= I.transpose()

print(l)
print(w)
print(I)
print(we)

#Need loop to make array
 
i= 30
d= w - we[i]
da= np.abs(d)
mda= np.min(da)
j= np.flatnonzero(da==mda)
sfloat= mda/d[0,j]
s= sfloat.astype(int) 

#print(d)
#print(da)
#print(s)

if np.min(d)==0:
    In= I[0,j]
elif s==1:
    In= (I[0,j + 1] - I[0,j])*(we[i] - w[0,j])/(w[0,j + 1] - w[0,j]) + I[0,j]
else:
    In= (I[0,j] - I[0,j - 1])*(we[i] - w[0,j - 1])/(w[0,j] - w[0,j - 1]) + I[0,j - 1]
print(we[i])
print(In)

#plt.subplot(211)
#plt.plot(l[0].transpose(),I[0].transpose())
#plt.subplot(212)
plt.plot(w[0].transpose(),I[0].transpose())