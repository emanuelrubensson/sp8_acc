#!/usr/bin/env python3
import numpy as np
from scipy.io import mmread
from scipy.linalg import cholesky
import sp8py
import sp8py.sp2 as sp2
import sp8py.sp8 as sp8
import sp8py.sp7 as sp7
import sp8py.test_utils
import matplotlib.pyplot as plt

def check_result(X,D,Z):
    # D = Z*D_ort*Z'
    ZX   = np.matmul(Z,X) # Z*X
    ZXZT = np.matmul(ZX,Z.T)
    A = D - 2*ZXZT
    assert(np.max(np.abs(A)) < 1e-6)

def plot_convergence_sp2(X, D, Z, L_out, L_in, H_in, H_out, ax_cond, ax_idem, marker, color, linestyle, label=''):
    # Get SP2-ACC polynomials
    nmin,nmax,p,alpha,gap = sp2.get_sp2_polys(L_out, L_in, H_in, H_out)
    # Run SP2-ACC
    X0 = np.copy(X)
    sp2_out = sp2.sp2_acc(X0,nmin,nmax,p,alpha,expensive_output=1)
    # sp2_out: X,n,polys,nmul_vec,idem_err_trace,idem_err_maxabs
    check_result(sp2_out[0],D,Z)
    nmul_vec = sp2_out[3]
    idem_err_maxabs = sp2_out[5]
    # Run SP2-ACC again without breaking due to parameterfree stopping criterion
    X0 = np.copy(X)
    sp2_out_nobreak = sp2.sp2_acc(X0,nmin,nmax,p,alpha,expensive_output=1,no_break=1)
    nmul_vec_poststop = sp2_out_nobreak[3][len(nmul_vec)-1:]
    idem_err_maxabs_poststop = sp2_out_nobreak[5][len(nmul_vec)-1:]
    ax_cond.semilogy(nmul_vec,[1/x for x in gap[0:len(nmul_vec)]],
                     color=color, marker=marker, linestyle=linestyle,
                     linewidth=1, markersize=3, label=label)
    ax_idem.semilogy(nmul_vec_poststop[0:2],idem_err_maxabs_poststop[0:2],
                     color=color, linestyle='dotted', linewidth=1) # to add line between pre- and post-stop
    ax_idem.semilogy(nmul_vec_poststop[1:],idem_err_maxabs_poststop[1:],
                     color=color, marker='x', linestyle='dotted',
                     linewidth=1, markersize=4) # post stop
    ax_idem.semilogy(nmul_vec,idem_err_maxabs,
                     color=color, marker=marker, linestyle=linestyle,
                     linewidth=1, markersize=3) # pre stop

def plot_convergence_sp8(X, D, Z, L_out, L_in, H_in, H_out, ax_cond, ax_idem, marker, color, linestyle, label=''):
    X0 = np.copy(X)
    sp8_out = sp8.sp8_acc(X0, L_out, L_in, H_in, H_out, expensive_output=1)
    # sp8_out: X, nmul, polys, gap, nmul_vec, idem_err_trace, idem_err_maxabs
    #          0  1     2      3    4         5               6
    check_result(sp8_out[0],D,Z)
    gap = sp8_out[3]
    nmul_vec = sp8_out[4]
    idem_err_maxabs = sp8_out[6]
    X0 = np.copy(X)
    sp8_out_extraiter = sp8.sp8_acc(X0, L_out, L_in, H_in, H_out, expensive_output=1,fixed_niter=len(nmul_vec)+2)
    nmul_vec_poststop = sp8_out_extraiter[4][len(nmul_vec)-1:]
    idem_err_maxabs_poststop = sp8_out_extraiter[6][len(nmul_vec)-1:]
    ax_cond.semilogy(nmul_vec,[1/x for x in gap[0:len(nmul_vec)]],
                     color=color, marker=marker, linestyle=linestyle,
                     linewidth=1, markersize=4, label=label)
    ax_idem.semilogy(nmul_vec_poststop[0:2],idem_err_maxabs_poststop[0:2],
                     color=color, linestyle='dotted', linewidth=1) # to add line between pre- and post-stopstop
    ax_idem.semilogy(nmul_vec_poststop[1:],idem_err_maxabs_poststop[1:],
                     color=color, marker='x', linestyle='dotted',
                     linewidth=1, markersize=4) # post stop
    ax_idem.semilogy(nmul_vec,idem_err_maxabs,
                     color=color, marker=marker, linestyle=linestyle,
                     linewidth=1, markersize=4) # pre stop

def plot_convergence_sp7(X, D, Z, L_out, L_in, H_in, H_out, ax_cond, ax_idem, marker, color, linestyle, label=''):
    X0 = np.copy(X)
    sp7_out = sp7.sp7_acc(X0, L_out, L_in, H_in, H_out, expensive_output=1)
    # sp7_out: X, nmul, polys, gap, nmul_vec, idem_err_trace, idem_err_maxabs
    #          0  1     2      3    4         5               6
    check_result(sp7_out[0],D,Z)
    gap = sp7_out[3]
    nmul_vec = sp7_out[4]
    idem_err_maxabs = sp7_out[6]
    X0 = np.copy(X)
    sp7_out_extraiter = sp7.sp7_acc(X0, L_out, L_in, H_in, H_out, expensive_output=1,fixed_niter=len(nmul_vec)+2)
    nmul_vec_poststop = sp7_out_extraiter[4][len(nmul_vec)-1:]
    idem_err_maxabs_poststop = sp7_out_extraiter[6][len(nmul_vec)-1:]
    ax_cond.semilogy(nmul_vec,[1/x for x in gap[0:len(nmul_vec)]],
                     color=color, marker=marker, linestyle=linestyle,
                     linewidth=1, markersize=3, label=label)
    ax_idem.semilogy(nmul_vec_poststop[0:2],idem_err_maxabs_poststop[0:2],
                     color=color, linestyle='dotted', linewidth=1) # to add line between pre- and post-stopstop
    ax_idem.semilogy(nmul_vec_poststop[1:],idem_err_maxabs_poststop[1:],
                     color=color, marker='x', linestyle='dotted',
                     linewidth=1, markersize=4) # post stop
    ax_idem.semilogy(nmul_vec,idem_err_maxabs,
                     color=color, marker=marker, linestyle=linestyle,
                     linewidth=1, markersize=3) # pre stop


# float_type = np.float64
# n    = 100
# nocc = 40
# L_inner = 0.495
# H_inner = 0.505
# L_inner = 0.7143989815234746 # Aspirin
# H_inner = 0.7197035316060708 
# L_outer = L_inner
# H_outer = H_inner
# X0,D = sp8py.test_utils.get_random_X(n,nocc,L_inner,H_inner,float_type)

#### Aspirin
mtx_dirname='aspirin/aspirin_mtx/'
# INDE HOMO interval : [   -0.357278101527   -0.348181679848 ]
H_outer, H_inner = -0.357278101527, -0.348181679848
# INDE LUMO interval : [    0.035549679088    0.056525112327 ]
L_inner, L_outer = 0.035549679088, 0.056525112327
# INDE Spectrum of F: 	 [ -20.624839086600 , 51.715199740629 ]
lmin, lmax = -20.624839086600, 51.715199740629
#
F_sparse = mmread(mtx_dirname+'F_matrix_final.mtx')
F = F_sparse.todense()
D_sparse = mmread(mtx_dirname+'D_matrix_final.mtx')
D = D_sparse.todense()
S_sparse = mmread(mtx_dirname+'S_matrix.mtx')
S = S_sparse.todense()
# F_ort = Z'*F*Z
U = cholesky(S, lower=False)
Z = np.linalg.inv(U)
ZTF  = np.matmul(Z.T,F)
X0 = np.matmul(ZTF,Z)

# X = (lmax*I-X)/(lmax-lmin) = (X-lmax*I)/(lmin-lmax)
np.fill_diagonal(X0, np.diag(X0) - lmax)
np.multiply(X0,1/(lmin-lmax), out=X0)
X0 = np.asanyarray(X0, order='F', dtype=np.float64)
H_inner = (lmax-H_inner)/(lmax-lmin)
H_outer = (lmax-H_outer)/(lmax-lmin)
L_inner = (lmax-L_inner)/(lmax-lmin)
L_outer = (lmax-L_outer)/(lmax-lmin)

plt.rc('font', size=12)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, (ax1,ax2) = plt.subplots(nrows = 2)

plot_convergence_sp2(X0, D, Z, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'D','tab:green', 'solid', 'SP2-ACC')

plot_convergence_sp2(X0, D, Z, 0.0, L_inner, H_inner, 1.0, ax1, ax2,'D','tab:green', 'dashed', 'SP2')

#X0_single = np.asanyarray(X0, order='F', dtype=np.float32)
#plot_convergence_sp2(X0_single, D, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'x','r')

# plot_convergence_sp8(X0, D, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'o','black', 'SP8-ACC',sp2_puri=1)
plot_convergence_sp8(X0, D, Z, L_outer, L_inner, H_inner, H_outer, ax1, ax2,'o','tab:blue', 'solid', 'SP8-ACC')

plot_convergence_sp8(X0, D, Z, 0.0, L_inner, H_inner, 1.0, ax1, ax2,'o','tab:blue', 'dashed', 'SP8')

plot_convergence_sp7(X0, D, Z, L_outer, L_inner, H_inner, H_outer, ax1, ax2, 's', 'tab:red', 'solid', 'SP7-ACC')
plot_convergence_sp7(X0, D, Z, 0.0, L_inner, H_inner, 1.0, ax1, ax2, 's', 'tab:red', 'dashed', 'SP7')

ax1.grid('on')
ax2.grid('on')
ax1.set_xlim( ax2.get_xlim() )
#ax1.set_title('SP2-ACC')
ax1.legend()
ax2.set_xlabel('No. of multiplications')
ax1.set_ylabel('Condition number')
ax2.set_ylabel('Idempotency error')
fig.savefig('plot_convergence.pdf', bbox_inches='tight')
