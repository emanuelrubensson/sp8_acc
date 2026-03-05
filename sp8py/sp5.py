import numpy as np
import sp8py
from scipy.optimize import root

def poly5_simple(mc,X):
    I = np.eye(X.shape[0],order='F',dtype=type(X[0][0]))
    Y = mc[5]*I
    Z = X
    Y += mc[4]*Z
    for i in [3,2,1,0]:
        Z = np.matmul(Z,X)
        Y += mc[i]*Z
    Y = np.asanyarray(Y, order='F', dtype=type(X[0][0]))
    return Y

def sp5_polynomial_coeffs(L):
    def polynomial_coefficients(s1, s2):
        A = 1.0 + s1 - s1**2 + s2 - s2**2
        B = s1**2 - s1 + s2**2 - s2
        C = s1 * s2 * (1.0 - s1) * (1.0 - s2)
        D = 1.0/30.0 + B/6.0 + C   # normalization from p(1)=1
        c = 1.0 / D
        a5 = c/5.0
        a4 = -c/2.0
        a3 = c*A/3.0
        a2 = c*B/2.0
        a1 = c*C
        a0 = 0.0
        return np.array([a5, a4, a3, a2, a1, a0])    
    def equations(s):
        """
        s: numpy array with s[0]=s1, s[1]=s2
        Output: numpy array [eq1, eq2] where
            eq1 : I(s1)
            eq2 : I(L) - I(s2)
        with I(x)=\int_0^x (t-s1)(t-s2)(t-(1-s1))(t-(1-s2)) dt.
        """
        s1 = s[0]
        s2 = s[1]
        # Coefficients for F(t)=t^4 - 2 t^3 + A t^2 + B t + C
        A = 1.0 + s1 - s1**2 + s2 - s2**2
        B = s1**2 - s1 + s2**2 - s2
        C = s1 * s2 * (1.0 - s1) * (1.0 - s2)
        def I(x):
            return (x**5) / 5.0 - (x**4) / 2.0 + (A / 3.0) * (x**3) + (B / 2.0) * (x**2) + C * x
        eq1 = I(s1)          # p(s1)=0  <=> I(s1)=0 since p(x)=c I(x) and c≠0
        eq2 = I(L) - I(s2)   # p(s2)=p(L) <=> I(s2)=I(L)
        return np.array([eq1, eq2])
    initial_guess = np.array([0.7*L, 0.2*L])
    sol = root(equations, initial_guess, method="hybr")
    if not sol.success:
        raise RuntimeError("Nonlinear solve failed: " + sol.message)
    return polynomial_coefficients(sol.x[0],sol.x[1])



def sp5_acc(X, L_outer, L_inner, H_inner, H_outer, expensive_output=0, fixed_niter=0):
    Csp5 = 16
    qsp5 = 3
    precond_limit = 0.01
    X2 = np.empty_like(X)
    mu = (L_inner+H_inner)/2
    polys = []
    nmul = 0
    nmul_vec = [nmul]
    gap = [H_inner-L_inner]
    idem_err_trace = [sp8py.trace_XmX2(X)]
    idem_err_maxabs = []    
    if expensive_output:
        sp8py.matmul(X, X, X2)
        idem_err_maxabs.append(np.max(np.abs(X-X2)))
    if mu > 0.5:        
        X = np.multiply(X,0.5/mu, out=X) # X = (0.5/mu)*X
        L_inner = (0.5/mu)*L_inner
        H_inner = (0.5/mu)*H_inner
        L_outer = (0.5/mu)*L_outer
        H_outer = (0.5/mu)*H_outer
        polys.append([0.5/mu, 0])
    elif mu < 0.5:
        np.fill_diagonal(X, np.diag(X) - 1) # X = X-I
        X = np.multiply(X,1/(2*(1-mu)), out=X) # X = (1/(2*(1-mu)))*X
        np.fill_diagonal(X, np.diag(X) + 1) # X = X+I
        L_inner = (L_inner-1)/(2*(1-mu)) + 1
        H_inner = (H_inner-1)/(2*(1-mu)) + 1
        L_outer = (L_outer-1)/(2*(1-mu)) + 1
        H_outer = (H_outer-1)/(2*(1-mu)) + 1
        polys.append([1/(2*(1-mu)), 1-1/(2*(1-mu))])
    else:
        polys.append([1,0])
    nmul_vec.append(nmul)
    gap.append(H_inner-L_inner)
    idem_err_trace.append( sp8py.trace_XmX2(X) )
    if expensive_output:
        sp8py.matmul(X, X, X2)
        idem_err_maxabs.append(np.max(np.abs(X-X2)))
    while L_outer > precond_limit:
        mc5 = sp5_polynomial_coeffs(L_outer)
        mc5_prim = np.polyder(mc5)
        print('SP5 slope at 0.5: ',np.polyval(mc5_prim,0.5))
        X = poly5_simple(mc5,X) # FIXME: USE SCHEME USING ONLY 3 MULTS
        polys.append(mc5)
        nmul += 3
        L_inner = np.polyval(mc5, L_inner)
        H_inner = np.polyval(mc5, H_inner)
        L_outer = np.polyval(mc5, L_outer)
        H_outer = np.polyval(mc5, H_outer)
        nmul_vec.append(nmul)
        gap.append(H_inner-L_inner)
        idem_err_trace.append( sp8py.trace_XmX2(X) )
        if expensive_output:
            sp8py.matmul(X, X, X2)
            idem_err_maxabs.append(np.max(np.abs(X-X2)))

    mc5 = [6.0, -15.0, 10.0, 0., 0., 0.]
    while 1:
        X = poly5_simple(mc5,X) # FIXME: USE SCHEME USING ONLY 3 MULTS
        polys.append(mc5)
        nmul += 3
        L_inner = np.polyval(mc5, L_inner)
        H_inner = np.polyval(mc5, H_inner)
        L_outer = np.polyval(mc5, L_outer)
        H_outer = np.polyval(mc5, H_outer)
        nmul_vec.append(nmul)
        gap.append(H_inner-L_inner)
        idem_err_trace.append( sp8py.trace_XmX2(X) )
        if expensive_output:
            sp8py.matmul(X, X, X2)
            idem_err_maxabs.append(np.max(np.abs(X-X2)))
        if fixed_niter:
            if len(nmul_vec) >= fixed_niter:
                break
            continue
        if idem_err_trace[-1] <= 0:
            break
        if idem_err_trace[-1] > Csp5*idem_err_trace[-2]**qsp5:
            break
    return X,nmul,polys,gap,nmul_vec,idem_err_trace,idem_err_maxabs
