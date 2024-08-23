#!/usr/bin/env python3
import numpy as np
import sp8py
import sp8py.sp2 as sp2
import sp8py.test_utils
import matplotlib.pyplot as plt

float_type = np.float64
eps   = np.finfo(float_type).eps
n    = 4
nocc = 2
L_outer = 0.494
L_inner = 0.495
H_inner = 0.505
H_outer = 0.506

X0,D = sp8py.test_utils.get_random_X(n,nocc,L_inner,H_inner,float_type)
# Run SP2-ACC
nmin,nmax,p,alpha = sp2.get_sp2_polys(L_outer, L_inner, H_inner, H_outer)
X = X0
X,n,polys,nmul_vec,idem_err_trace,idem_err_maxabs = sp2.sp2_acc(X,nmin,nmax,p,alpha,expensive_output=1)


fig, ax = plt.subplots()
ax.semilogy(nmul_vec,idem_err_maxabs,'xb-')
ax.grid('on')
ax.set_title('SP2-ACC')
fig.savefig('plot_convergence.pdf', bbox_inches='tight')

