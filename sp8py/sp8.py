import numpy as np
import sp8py
#import sp2
import matplotlib.pyplot as plt

def poly8_simple(mc,X):
    I = np.eye(X.shape[0],order='F',dtype=type(X[0][0]))
    Y = mc[8]*I
    Z = X
    Y += mc[7]*Z
    for i in [6,5,4,3,2,1,0]:
        Z = np.matmul(Z,X)
        Y += mc[i]*Z
    Y = np.asanyarray(Y, order='F', dtype=type(X[0][0]))
    return Y

def sp8_acc(X, L_outer, L_inner, H_inner, H_outer, expensive_output=0):
    Csp8 = 85
    qsp8 = 4
    Csp8sp2 = 275
    qsp8sp2 = 5    
    precond_limit = 0.01
    X2 = np.empty_like(X)
    M2 = np.empty_like(X)
    M3 = np.empty_like(X)
    v  = np.empty(9)
    mc = np.empty(9)
    idem_err_trace = [sp8py.trace_XmX2(X)]
    idem_err_maxabs = []    
    precond_phase = (L_inner > precond_limit) or (H_inner < 1-precond_limit)
    polys = []
    print(idem_err_trace[-1])

    # Preconditioning phase
    print('Preconditioning phase')
    nmul = 0
    nmul_vec = [nmul]
    while precond_phase:
        sp8py.get_sp8_params_max_gap(L_outer, L_inner, H_inner, H_outer, v)
        sp8py.get_sp8_monomial_coefficients(v, mc)
        # Evaluate SP8 for X
        sp8py.matmul(X, X, X2)
        if expensive_output:
            idem_err_maxabs.append(sp8py.maxabs(X-X2))
        nmul += 1
        sp8py.poly_8_eval(mc, X, X2, M2, M3)        
        # X = poly8_simple(mc,X)
        nmul += 2
        polys.append(list(mc))
        # Update homo lumo bounds accordingly
        L_outer = np.polyval(mc, L_outer)
        L_inner = np.polyval(mc, L_inner)
        H_inner = np.polyval(mc, H_inner)
        H_outer = np.polyval(mc, H_outer)
        L_outer = max(L_outer,0.0)
        H_outer = min(H_outer,1.0)
        # Compute Tr[X-X**2] (without computing X**2)
        nmul_vec.append(nmul)
        idem_err_trace.append( sp8py.trace_XmX2(X) )
        precond_phase = ((L_inner > precond_limit) or
                         (H_inner < 1-precond_limit))
        print(idem_err_trace[-1],L_inner,H_inner)
    # Purification phase
    print('Purification phase')

# ALTERNATIVE TO RUN SP2 IN PURIFICATION PHASE
#    nmin,nmax,p,alpha = sp2.get_sp2_polys(L_outer, L_inner, H_inner, H_outer)
#    X,npuri,polys_puri = sp2.sp2_acc(X,nmin,nmax,p,alpha)
#    n+=npuri
#    return X,n
    
    sp8py.matmul(X, X, X2)
    if expensive_output:
        idem_err_maxabs.append(sp8py.maxabs(X-X2))
    nmul += 1
    while 1:
        # Need to do at least 1 iteration in purification phase
        if (L_inner + H_inner > 1):
            left  = 4
            right = 3
        else:
            left  = 3
            right = 4
        info = sp8py.get_sp8_params(0,1,(left, right,False,False), v)
        assert(info==0)
        sp8py.get_sp8_monomial_coefficients(v, mc)
        # Evaluate SP8 for X
        sp8py.poly_8_eval(mc, X, X2, M2, M3)
        # X = poly8_simple(mc,X)
        nmul += 2
        polys.append(list(mc))
        # Update inner homo lumo bounds accordingly
        L_inner = np.polyval(mc, L_inner)
        H_inner = np.polyval(mc, H_inner)
        # Compute Tr[X-X**2] (without computing X**2)
        nmul_vec.append(nmul)
        idem_err_trace.append( sp8py.trace_XmX2(X) )
        print(f'IDEM: {idem_err_trace[-1]}')
        if idem_err_trace[-1] <= 0:
            break
        # First convergence check, error reduction with applied 8th-degree
        print(f'A-test: {idem_err_trace[-1]} >? {Csp8*idem_err_trace[-2]**qsp8} ---- {idem_err_trace[-2]}')
        if idem_err_trace[-1] > Csp8*idem_err_trace[-2]**qsp8:
            print(f'A: {idem_err_trace[-1]} > {Csp8*idem_err_trace[-2]**qsp8} ---- {idem_err_trace[-2]}')
            break
        # Attempt with only one more multiply
        sp8py.matmul(X, X, X2)
        if expensive_output:
            idem_err_maxabs.append(sp8py.maxabs(X-X2))
        nmul += 1
        # if (L_inner + H_inner < 1):
        if left == 4:    # attempt with 2x - x**2
            # M2 = 2*X-X2
            np.multiply(2.0,X,out=M3)
            np.subtract(M3,X2,out=M2)            
            poly = [-1, 2, 0]
        else:            # attempt with x**2
            # M2 = X2
            M2 = np.copy(X2)
            poly = [1, 0, 0]
        M2_idem_err_trace = sp8py.trace_XmX2(M2)
        # Second convergence check,
        # error reduction with applied 8th-degree + 2nd degree
        print(f'B-test: {M2_idem_err_trace} >? {Csp8sp2*idem_err_trace[-2]**qsp8sp2} ---- {idem_err_trace[-2]}')
        if M2_idem_err_trace <= 0 or M2_idem_err_trace > Csp8sp2*idem_err_trace[-2]**qsp8sp2:
            print(f'B: {M2_idem_err_trace} > {Csp8sp2*idem_err_trace[-2]**qsp8sp2} ---- {idem_err_trace[-2]}')
            X = np.copy(M2)
            nmul_vec.append(nmul)
            idem_err_trace.append(M2_idem_err_trace)
            polys.append(poly)            
            break
        # continue
        print(idem_err_trace[-1],L_inner,H_inner)
    if expensive_output:
        sp8py.matmul(X,X,X2)
        idem_err_maxabs.append(sp8py.maxabs(X-X2))
        return X,nmul,polys,nmul_vec,idem_err_trace,idem_err_maxabs
    return X,nmul,polys,nmul_vec,idem_err_trace
