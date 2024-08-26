#!/usr/bin/env python3
import numpy as np
import sp8py
import sp8py.sp2 as sp2
import sp8py.test_utils

def run_test_sp2_acc(X, l_outer, l_inner, h_inner, h_outer, Dref, verbose_output):
    nmin,nmax,p,alpha,gap = sp2.get_sp2_polys(l_outer, l_inner, h_inner, h_outer)
    D,nmul,polys,nmul_vec,idem_err_trace = sp2.sp2_acc(X,nmin,nmax,p,alpha)
    norm_diff_fro = np.linalg.norm(D - Dref, ord = 'fro')
    norm_diff_max = sp8py.maxabs(D - Dref) 
    if verbose_output:
        print(f'SP2-ACC converged using {nmul} multiplications.')
        print(f'||D-Dref||_F   = {norm_diff_fro}')
        print(f'||D-Dref||_max = {norm_diff_max}')
    assert( norm_diff_max < np.sqrt(eps) )

verbose_output = False
float_type = np.float64
eps   = np.finfo(float_type).eps
n    = 4
nocc = 2
lumo = 0.49
homo = 0.51
X, Dref = sp8py.test_utils.get_random_X(n,nocc,lumo,homo,float_type)

if verbose_output:
    print('Running SP2-ACC')
run_test_sp2_acc(X, lumo, lumo, homo, homo, Dref, verbose_output)

if verbose_output:
    print('Running SP2-ACC without acceleration')
run_test_sp2_acc(X, 0,    lumo, homo, 1.0,  Dref, verbose_output)

# Make sure method does not fail with trivial input
X = np.array([[1.0,0],[0,1.0]])
X = np.asanyarray(X, order='F', dtype=float_type)
nmin,nmax,p,alpha,gap = sp2.get_sp2_polys(0.0, 0.0, 1.0, 1.0)
D,nmul,polys,nmul_vec,idem_err_trace = sp2.sp2_acc(X,nmin,nmax,p,alpha)
