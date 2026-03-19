#!/usr/bin/env python3
import numpy as np
import sp8py.test_utils
import sp8py.sp2 as sp2
import sp8py.sp8 as sp8
import sp8py.sp5 as sp5
import matplotlib.pyplot as plt
float_type = np.float64
n    = 400
nocc = 200
mu = 0.5
gap_start = -10
gap_end = -1
gap_vec = np.logspace(gap_start, gap_end, 10)
su5_gap_limit = 1e-10

sp2_acc_nmul = []
sp2_nmul     = []
sp8_acc_nmul = []
sp8_nmul     = []
sp5_acc_nmul = []
sp5_nmul     = []
sp2_acc_err = []
sp2_err     = []
sp8_acc_err = []
sp8_err     = []
sp5_acc_err = []
sp5_err     = []
for gap in gap_vec:
    homo = mu+gap/2
    lumo = mu-gap/2
    L_inner = lumo
    H_inner = homo
    L_outer = L_inner
    H_outer = H_inner
    # Construct matrix
    X,D = sp8py.test_utils.get_random_X(n,nocc,L_inner,H_inner,float_type)
    #X,D = sp8py.test_utils.get_diagonal_X(n,nocc,lumo,homo)
    ###### Run methods and save nmul
    ### SP2-ACC
    # Get SP2-ACC polynomials
    nmin,nmax,p,alpha,gapvec = sp2.get_sp2_polys(L_outer, L_inner, H_inner, H_outer)
    # Run SP2-ACC
    X0 = np.copy(X)
    sp2_out = sp2.sp2_acc(X0,nmin,nmax,p,alpha)
    sp2_acc_nmul.append(sp2_out[1])
    sp2_acc_err.append(np.max(np.abs(sp2_out[0] - D)))
    ### SP2 
    # Get SP2 polynomials
    nmin,nmax,p,alpha,gapvec = sp2.get_sp2_polys(0.0, L_inner, H_inner, 1.0)
    # Run SP2
    X0 = np.copy(X)
    sp2_out = sp2.sp2_acc(X0,nmin,nmax,p,alpha)
    sp2_nmul.append(sp2_out[1])
    sp2_err.append(np.max(np.abs(sp2_out[0] - D)))
    ### SP8-ACC
    X0 = np.copy(X)
    sp8_out = sp8.sp8_acc(X0, L_outer, L_inner, H_inner, H_outer)
    sp8_acc_nmul.append(sp8_out[1])
    sp8_acc_err.append(np.max(np.abs(sp8_out[0] - D)))    
    ### SP8
    X0 = np.copy(X)
    sp8_out = sp8.sp8_acc(X0, 0.0, L_inner, H_inner, 1.0)
    sp8_nmul.append(sp8_out[1])
    sp8_err.append(np.max(np.abs(sp8_out[0] - D)))    
    ### SP5-ACC (Su5)
    if gap >= su5_gap_limit: # Probably not needed anymore: problem with unstable fixed points solved
        X0 = np.copy(X)
        sp5_out = sp5.sp5_acc(X0, L_outer, L_inner, H_inner, H_outer)
        sp5_acc_nmul.append(sp5_out[1])
        sp5_acc_err.append(np.max(np.abs(sp5_out[0] - D)))    
    ### SP5
    X0 = np.copy(X)
    sp5_out = sp5.sp5_acc(X0, 0.0, L_inner, H_inner, 1.0, expensive_output=1)
    sp5_nmul.append(sp5_out[1])
    sp5_err.append(np.max(np.abs(sp5_out[0] - D)))    

su5_gap_vec = [gap for gap in gap_vec if gap >= su5_gap_limit]

plt.rc('font', size=12)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fig, (ax1,ax2) = plt.subplots(nrows = 2)

ax1.semilogx(gap_vec,sp2_nmul,'--D',markersize=3,color = 'tab:green', label='SP2')
ax1.semilogx(gap_vec,sp2_acc_nmul,'-D',markersize=3,color = 'tab:green', label='SP2-ACC')
ax1.semilogx(su5_gap_vec,sp5_acc_nmul,'-s',markersize=3,color = 'tab:red', label='Su5')
ax1.semilogx(gap_vec,sp8_nmul,'--o',markersize=4,color = 'tab:blue', label='SP8')
ax1.semilogx(gap_vec,sp8_acc_nmul,'-o',markersize=4,color = 'tab:blue', label='SP8-ACC')
#ax1.semilogx(gap_vec,sp5_nmul,'--s',color = 'tab:red', label='SP5')

ax2.loglog(gap_vec,sp2_err,'--D',markersize=3,color = 'tab:green', label='SP2')
ax2.loglog(gap_vec,sp2_acc_err,'-D',markersize=3,color = 'tab:green', label='SP2-ACC')
ax2.loglog(su5_gap_vec,sp5_acc_err,'-s',markersize=3,color = 'tab:red', label='Su5')
ax2.loglog(gap_vec,sp8_err,'--o',markersize=4,color = 'tab:blue', label='SP8')
ax2.loglog(gap_vec,sp8_acc_err,'-o',markersize=4,color = 'tab:blue', label='SP8-ACC')
#ax2.loglog(gap_vec,sp5_err,'--s',color = 'tab:red', label='SP5')

meps = np.finfo(np.float64).eps
ax2.loglog(gap_vec,100*meps/gap_vec,'--',markersize=4,color = 'black')


ax1.grid('on')
ax2.grid('on')
ax1.set_xlim( ax2.get_xlim() )

ax1.legend()
ax2.set_xlabel(r'Homo-lumo gap $\xi$')
ax1.set_ylabel('No. of multiplications')
ax2.set_ylabel('Max abs error')
fig.savefig('plot_gap_vs_nmul.pdf', bbox_inches='tight')



