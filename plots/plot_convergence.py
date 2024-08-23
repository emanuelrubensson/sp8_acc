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
# Get SP2-ACC polynomials
nmin,nmax,p,alpha = sp2.get_sp2_polys(L_outer, L_inner, H_inner, H_outer)
X = X0
# Run SP2-ACC
sp2_out = sp2.sp2_acc(X,nmin,nmax,p,alpha,expensive_output=1)
# sp2_out: X,n,polys,nmul_vec,idem_err_trace,idem_err_maxabs
nmul_vec = sp2_out[3]
idem_err_maxabs = sp2_out[5]
# Run SP2-ACC again without breaking due to parameterfree stopping criterion
sp2_out_nobreak = sp2.sp2_acc(X,nmin,nmax,p,alpha,expensive_output=1,no_break=1)
nmul_vec_poststop = sp2_out_nobreak[3][len(nmul_vec)-1:]
idem_err_maxabs_poststop = sp2_out_nobreak[5][len(nmul_vec)-1:]


fig, ax = plt.subplots()
ax.semilogy(nmul_vec_poststop[0:2],idem_err_maxabs_poststop[0:2],'r-')
ax.semilogy(nmul_vec_poststop[1:],idem_err_maxabs_poststop[1:],'sr-')
ax.semilogy(nmul_vec,idem_err_maxabs,'xb-')
ax.grid('on')
ax.set_title('SP2-ACC')
fig.savefig('plot_convergence.pdf', bbox_inches='tight')

