import numpy as np
import sp8py
import matplotlib.pyplot as plt

def sp2_acc(X,nmin,nmax,p,alpha):
    """
    X: input matrix with all eigenvalues in [0,1]
    """
    nmul = 0
    Csp2 = 4.5    
    I = np.eye(X.shape[0],order='F',dtype=type(X[0][0]))
    # check properties with I.flags
    X2 = np.empty_like(X)
    nmul_vec = [nmul]
    idem_err_trace = [sp8py.trace_XmX2(X)]
    polys = []
    for i in range(nmax):
        a = alpha[i]
        # Using own multiply to make sure same routine is used throughout. 
        # (rather than np.matmul(X,X,out=X2))
        sp8py.matmul(X,X,X2)
        nmul += 1
        if p[i]:
            X = ((1-a)**2)*I + 2*(1-a)*a*X + (a**2)*X2
            polys.append([a**2, 2*(1-a)*a, (1-a)**2])
        else:
            X = 2*a*X - (a**2)*X2
            polys.append([-a**2, 2*a, 0])
        nmul_vec.append(nmul)
        idem_err_trace.append(sp8py.trace_XmX2(X))
        if idem_err_trace[-1] <= 0:
            break
        if i >= nmin and p[i] != p[i-1] and idem_err_trace[-1] > Csp2*idem_err_trace[-3]**2:
            break
    return X,nmul,polys,nmul_vec,idem_err_trace

def get_sp2_polys(L_outer, L_inner, H_inner, H_outer):
    """ 
    Implements Alg. 4 from http://dx.doi.org/10.1021/acs.jctc.6b00626
    L_outer, L_inner: outer and inner lumo eigenvalue bounds
    H_inner, H_outer: inner and outer homo eigenvalue bounds
    """
    eps   = np.finfo(type(L_inner)).eps
    l_lo  = L_outer       # gamma
    l_up  = L_inner
    h_up  = 1 - H_inner   # beta
    h_lo  = 1 - H_outer
    delta = 0.01
    p     = []
    alpha = []
    i = 0
    while i < 2 or h_up-h_up**2 > eps*eps or l_up-l_up**2 > eps*eps:
        i += 1
        if h_lo < delta and l_lo < delta:
            h_lo = 0
            l_lo = 0
            nmin = max(i,2)
            delta = 0
        if l_up > h_up or (i%2 != 0 and l_up == h_up):
            # X = X**2 (but with scaling)
            p.append(True)
            a = 2/(2-l_lo)
            l_lo = ((1-a)+a*l_lo)**2
            l_up = ((1-a)+a*l_up)**2
            h_lo = 2*a*h_lo - (a*h_lo)**2
            h_up = 2*a*h_up - (a*h_up)**2
        else:
            # X = 2*X - X**2 (but with scaling)
            p.append(False)
            a = 2/(2-h_lo)
            h_lo = ((1-a)+a*h_lo)**2
            h_up = ((1-a)+a*h_up)**2
            l_lo = 2*a*l_lo - (a*l_lo)**2
            l_up = 2*a*l_up - (a*l_up)**2
        alpha.append(a)
    nmax = i
    return (nmin,nmax,p,alpha)
