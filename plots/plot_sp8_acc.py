#!/usr/bin/env python3
from sp8_acc import *
import numpy as np
import matplotlib.pyplot as plt

lumo = 0.5;
homo = 0.6;
v = np.empty(9)
get_sp8_params(lumo,homo,4,3,v)

mc = np.empty(9)
get_sp8_monomial_coefficients(v,mc)

xv = np.linspace(0,1,500)
yv = [np.polyval(mc, x) for x in xv]
plt.plot(xv,yv)
plt.plot(lumo,np.polyval(mc, lumo),'xb')
plt.plot(homo,np.polyval(mc, homo),'xb')
plt.grid('on')
plt.show()
