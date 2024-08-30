#!/usr/bin/env python3
import numpy as np
import sp8py
import sp8py.sp2 as sp2
import sp8py.test_utils
import matplotlib.pyplot as plt

def plot_convergence_sp2(X, D, L_out, L_in, H_in, H_out, ax_cond, ax_idem, marker, color):
    # Get SP2-ACC polynomials
    nmin,nmax,p,alpha,gap = sp2.get_sp2_polys(L_out, L_in, H_in, H_out)
    # Run SP2-ACC
    sp2_out = sp2.sp2_acc(X,nmin,nmax,p,alpha,expensive_output=1)
    # sp2_out: X,n,polys,nmul_vec,idem_err_trace,idem_err_maxabs
    nmul_vec = sp2_out[3]
    idem_err_maxabs = sp2_out[5]
    # Run SP2-ACC again without breaking due to parameterfree stopping criterion
    sp2_out_nobreak = sp2.sp2_acc(X,nmin,nmax,p,alpha,expensive_output=1,no_break=1)
    nmul_vec_poststop = sp2_out_nobreak[3][len(nmul_vec)-1:]
    idem_err_maxabs_poststop = sp2_out_nobreak[5][len(nmul_vec)-1:]
    ax_cond.semilogy(nmul_vec,[1/x for x in gap[0:len(nmul_vec)]],marker+color+'-')
    ax_idem.semilogy(nmul_vec_poststop[0:2],idem_err_maxabs_poststop[0:2],color+'-')
    ax_idem.semilogy(nmul_vec_poststop[1:],idem_err_maxabs_poststop[1:],'s'+color+'-')
    ax_idem.semilogy(nmul_vec,idem_err_maxabs,marker+color+'-')


fig, (ax1,ax2) = plt.subplots(nrows = 2)

float_type = np.float64
n    = 4
nocc = 2
L_outer = 0.4995
L_inner = 0.4995
H_inner = 0.5005
H_outer = 0.5005
X0,D = sp8py.test_utils.get_random_X(n,nocc,L_inner,H_inner,float_type)

plot_convergence_sp2(X0, D, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'x','r')

plot_convergence_sp2(X0, D, 0.0, L_inner, H_inner, 1.0, ax1, ax2,'+','g')

X0_single = np.asanyarray(X0, order='F', dtype=np.float32)
plot_convergence_sp2(X0_single, D, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'x','r')


ax1.grid('on')
ax2.grid('on')
ax1.set_xlim( ax2.get_xlim() )
ax1.set_title('SP2-ACC')
ax2.set_xlabel('No. of multiplications')
ax1.set_ylabel('Condition number')
ax2.set_ylabel('Idempotency error')
fig.savefig('plot_convergence.pdf', bbox_inches='tight')
