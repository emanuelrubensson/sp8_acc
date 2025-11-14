#!/usr/bin/env python3
import numpy as np
import sp8py
import sp8py.sp2 as sp2
import sp8py.sp8 as sp8
import sp8py.test_utils
import matplotlib.pyplot as plt

def plot_convergence_sp2(X, D, L_out, L_in, H_in, H_out, ax_cond, ax_idem, marker, color, label=''):
    # Get SP2-ACC polynomials
    nmin,nmax,p,alpha,gap = sp2.get_sp2_polys(L_out, L_in, H_in, H_out)
    # Run SP2-ACC
    X0 = np.copy(X)
    sp2_out = sp2.sp2_acc(X0,nmin,nmax,p,alpha,expensive_output=1)
    # sp2_out: X,n,polys,nmul_vec,idem_err_trace,idem_err_maxabs
    nmul_vec = sp2_out[3]
    idem_err_maxabs = sp2_out[5]
    # Run SP2-ACC again without breaking due to parameterfree stopping criterion
    X0 = np.copy(X)
    sp2_out_nobreak = sp2.sp2_acc(X0,nmin,nmax,p,alpha,expensive_output=1,no_break=1)
    nmul_vec_poststop = sp2_out_nobreak[3][len(nmul_vec)-1:]
    idem_err_maxabs_poststop = sp2_out_nobreak[5][len(nmul_vec)-1:]
    ax_cond.semilogy(nmul_vec,[1/x for x in gap[0:len(nmul_vec)]],
                     color=color, marker=marker, linestyle='solid',
                     linewidth=1, markersize=4, label=label)
    ax_idem.semilogy(nmul_vec_poststop[0:2],idem_err_maxabs_poststop[0:2],
                     color=color, linestyle='solid', linewidth=1) # to add line between pre- and post-stop
    ax_idem.semilogy(nmul_vec_poststop[1:],idem_err_maxabs_poststop[1:],
                     color=color, marker='s', linestyle='solid',
                     linewidth=1, markersize=4) # post stop
    ax_idem.semilogy(nmul_vec,idem_err_maxabs,
                     color=color, marker=marker, linestyle='solid',
                     linewidth=1, markersize=4) # pre stop

def plot_convergence_sp8(X, D, L_out, L_in, H_in, H_out, ax_cond, ax_idem, marker, color, label='',sp2_puri=0):
    X0 = np.copy(X)
    sp8_out = sp8.sp8_acc(X0, L_out, L_in, H_in, H_out, expensive_output=1, sp2_purification=sp2_puri)
    # sp8_out: X, nmul, polys, gap, nmul_vec, idem_err_trace, idem_err_maxabs
    #          0  1     2      3    4         5               6
    gap = sp8_out[3]
    nmul_vec = sp8_out[4]
    idem_err_maxabs = sp8_out[6]
    X0 = np.copy(X)
    sp8_out_extraiter = sp8.sp8_acc(X0, L_out, L_in, H_in, H_out, expensive_output=1,fixed_niter=len(nmul_vec)+3, sp2_purification=sp2_puri)
    nmul_vec_poststop = sp8_out_extraiter[4][len(nmul_vec)-1:]
    idem_err_maxabs_poststop = sp8_out_extraiter[6][len(nmul_vec)-1:]
    ax_cond.semilogy(nmul_vec,[1/x for x in gap[0:len(nmul_vec)]],
                     color=color, marker=marker, linestyle='solid',
                     linewidth=1, markersize=4, label=label)
    ax_idem.semilogy(nmul_vec_poststop[0:2],idem_err_maxabs_poststop[0:2],
                     color=color, linestyle='solid', linewidth=1) # to add line between pre- and post-stopstop
    ax_idem.semilogy(nmul_vec_poststop[1:],idem_err_maxabs_poststop[1:],
                     color=color, marker='s', linestyle='solid',
                     linewidth=1, markersize=4) # post stop
    ax_idem.semilogy(nmul_vec,idem_err_maxabs,
                     color=color, marker=marker, linestyle='solid',
                     linewidth=1, markersize=4) # pre stop
    
fig, (ax1,ax2) = plt.subplots(nrows = 2)

float_type = np.float64
n    = 4
nocc = 2
L_inner = 0.4995
H_inner = 0.5005
L_outer = L_inner
H_outer = H_inner
X0,D = sp8py.test_utils.get_random_X(n,nocc,L_inner,H_inner,float_type)

plot_convergence_sp2(X0, D, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'x','tab:blue', 'SP2-ACC')

plot_convergence_sp2(X0, D, 0.0, L_inner, H_inner, 1.0, ax1, ax2,'+','tab:green', 'SP2')

#X0_single = np.asanyarray(X0, order='F', dtype=np.float32)
#plot_convergence_sp2(X0_single, D, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'x','r')

# plot_convergence_sp8(X0, D, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'o','black', 'SP8-ACC',sp2_puri=1)
plot_convergence_sp8(X0, D, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'o','grey', 'SP8-ACC',sp2_puri=0)

plot_convergence_sp8(X0, D, 0.0, L_inner, H_inner, 1.0, ax1, ax2,'o','tab:orange', 'SP8')


ax1.grid('on')
ax2.grid('on')
ax1.set_xlim( ax2.get_xlim() )
ax1.set_title('SP2-ACC')
ax1.legend()
ax2.set_xlabel('No. of multiplications')
ax1.set_ylabel('Condition number')
ax2.set_ylabel('Idempotency error')
fig.savefig('plot_convergence.pdf', bbox_inches='tight')
