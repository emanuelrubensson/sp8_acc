#!/usr/bin/env python3
from sp8_acc import *
import numpy as np
import matplotlib.pyplot as plt

def plot_sp8(homo, lumo, left, right, filename):
    v = np.empty(9)
    acc_left = True
    acc_right = True
    sp8_spec = (left,right,acc_left,acc_right)
    get_sp8_params(lumo,homo,sp8_spec,v)
    print(f'Slope <{left},{right}>: {sp8_prim(v,homo)}')
#    print(f'v: {[ "{:0.16f}".format(x) for x in v]}')
    mc = np.empty(9)
    get_sp8_monomial_coefficients(v,mc)
    xv = np.linspace(-0.1,1.1,500)
    yv = [np.polyval(mc, x) for x in xv]

    fig, ax = plt.subplots()
    ax.plot(xv,yv)
    ax.plot(lumo,np.polyval(mc, lumo),'xb')
    ax.plot(homo,np.polyval(mc, homo),'xb')
    ax.grid('on')
    ax.axis([-0.1, 1.1,-0.1,1.1])
    ax.set_title(f'<{left},{right}>')
    fig.savefig(filename, bbox_inches='tight')

lumo = 0.23;
homo = 0.25;
left = 3
right = 4
plot_sp8(homo, lumo, left, right, "plot_sp8_acc_A.pdf")
#lumo = 0.25;
#left = 1
#right = 0
#plot_sp8(homo, lumo, left, right, "plot_sp8_acc_B.pdf")
#lumo = 0.05
#left = 0
#right = 7
#plot_sp8(homo, lumo, left, right, "plot_sp8_acc_C.pdf")
