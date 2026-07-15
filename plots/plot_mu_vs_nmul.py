#!/usr/bin/env python3
import numpy as np
import sp8py.test_utils
import sp8py.sp2 as sp2
import sp8py.sp8 as sp8
import sp8py.sp5 as sp5
import matplotlib.pyplot as plt
float_type = np.float64
n    = 400
# nocc = 200
gap = 0.005
mu_step = 0.05
mu_start = mu_step
mu_end = 1-mu_step
npoints = (round(1/mu_step) + 1) - 2
mu_vec = np.linspace(mu_start, mu_end, npoints)
sp2_acc_nmul = []
sp2_nmul     = []
sp8_acc_nmul = []
sp8_nmul     = []
sp5_acc_nmul = []
sp5_nmul     = []
for mu in mu_vec:
    print(f'---------- mu = {mu} ----------')
    homo = mu+gap/2
    lumo = mu-gap/2
    L_inner = lumo
    H_inner = homo
    L_outer = L_inner
    H_outer = H_inner
    nocc = round(n*(1.0-mu)) # nocc varied with mu 
    # Construct matrix
    X,D = sp8py.test_utils.get_random_X(n,nocc,L_inner,H_inner,float_type)
    # X,D = sp8py.test_utils.get_diagonal_X(n,nocc,lumo,homo)
    ###### Run methods and save nmul
    
    ### SP2-ACC
    print('Running SP2-ACC')
    # Get SP2-ACC polynomials
    nmin,nmax,p,alpha,gapvec = sp2.get_sp2_polys(L_outer, L_inner, H_inner, H_outer)
    # Run SP2-ACC
    X0 = np.copy(X)
    sp2_out = sp2.sp2_acc(X0,nmin,nmax,p,alpha)
    sp2_acc_nmul.append(sp2_out[1])
    print(f'  SP2-ACC error: {np.max(np.abs(sp2_out[0] - D))}')
    X2 = np.matmul(sp2_out[0],sp2_out[0])
    print(f'  SP2-ACC idem error: {np.max(np.abs(sp2_out[0] - X2))}')
    
    ### SP2 
    print('Running SP2')
    # Get SP2 polynomials
    nmin,nmax,p,alpha,gapvec = sp2.get_sp2_polys(0.0, L_inner, H_inner, 1.0)    
    # Run SP2
    X0 = np.copy(X)
    sp2_out = sp2.sp2_acc(X0,nmin,nmax,p,alpha)
    sp2_nmul.append(sp2_out[1])
    print(f'  SP2 error: {np.max(np.abs(sp2_out[0] - D))}')
    X2 = np.matmul(sp2_out[0],sp2_out[0])
    print(f'  SP2 idem error: {np.max(np.abs(sp2_out[0] - X2))}')
    
    ### SP8-ACC
    print('Running SP8-ACC')
    X0 = np.copy(X)
    sp8_out = sp8.sp8_acc(X0, L_outer, L_inner, H_inner, H_outer)
    sp8_acc_nmul.append(sp8_out[1])
    print(f'  SP8-ACC error: {np.max(np.abs(sp8_out[0] - D))}')
    X2 = np.matmul(sp8_out[0],sp8_out[0])
    print(f'  SP8-ACC idem error: {np.max(np.abs(sp8_out[0] - X2))}')
    
    ### SP8
    print('Running SP8')
    X0 = np.copy(X)
    sp8_out = sp8.sp8_acc(X0, 0.0, L_inner, H_inner, 1.0)
    sp8_nmul.append(sp8_out[1])
    print(f'  SP8 error: {np.max(np.abs(sp8_out[0] - D))}')
    X2 = np.matmul(sp8_out[0],sp8_out[0])
    print(f'  SP8 idem error: {np.max(np.abs(sp8_out[0] - X2))}')
    
    ### SP5-ACC
    print('Running Su5')
    X0 = np.copy(X)
    sp5_out = sp5.sp5_acc(X0, L_outer, L_inner, H_inner, H_outer)
    sp5_acc_nmul.append(sp5_out[1])
    print(f'  Su5 error: {np.max(np.abs(sp5_out[0] - D))}')
    X2 = np.matmul(sp5_out[0],sp5_out[0])
    print(f'  Su5 idem error: {np.max(np.abs(sp5_out[0] - X2))}')

    #    ### SP5
    #    X0 = np.copy(X)
    #    sp5_out = sp5.sp5_acc(X0, 0.0, L_inner, H_inner, 1.0, expensive_output=1)
    #    sp5_nmul.append(sp5_out[1])
    
plt.rc('font', size=12)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

plt.plot(mu_vec,sp2_nmul,'--D',markersize=3,color = 'tab:green', label='SP2')
plt.plot(mu_vec,sp2_acc_nmul,'-D',markersize=3,color = 'tab:green', label='SP2-ACC')
plt.plot(mu_vec,sp5_acc_nmul,'-s',markersize=3,color = 'tab:red', label='Su5')
plt.plot(mu_vec,sp8_nmul,'--o',markersize=4,color = 'tab:blue', label='SP8')
plt.plot(mu_vec,sp8_acc_nmul,'-o',markersize=4,color = 'tab:blue', label='SP8-ACC')
# plt.plot(mu_vec,sp5_nmul,'--s',color = 'tab:red', label='SP5')

plt.grid('on')
plt.xlabel(r'Step location $\mu$')
plt.ylabel('No. of multiplications')
plt.legend()
plt.savefig('plot_mu_vs_nmul.pdf', bbox_inches='tight')

