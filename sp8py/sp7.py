import numpy as np
import sp8py
from scipy.optimize import root
import matplotlib.pyplot as plt

# This code was written only for the purpose of comparing performance
# of sp7 schemes with sp8.

def poly7_simple(mc,X):
    I = np.eye(X.shape[0],order='F',dtype=type(X[0][0]))
    Y = mc[7]*I
    Z = X
    Y += mc[6]*Z
    for i in [5,4,3,2,1,0]:
        Z = np.matmul(Z,X)
        Y += mc[i]*Z
    Y = np.asanyarray(Y, order='F', dtype=type(X[0][0]))
    return Y

def sp7_polynomial_coeffs(L, initial_guess=None):
    """
    Compute monomial coefficients of the 7th-degree polynomial p(x)
    defined by:
      p'(x) = c * prod_{k=1}^3 (x-r_k)(x-(1-r_k))
    with constraints:
      p(0)=p(L), p(r1)=0, p(r2)=p(L), p(r3)=0, p(0.5)=0.5

    Parameters
    ----------
    L : float
        Must satisfy 0 < L < 0.5
    initial_guess : array_like, optional
        Initial guess for (s1, s2, s3). If None, a heuristic guess is used.

    Returns
    -------
    coeffs : ndarray
        Monomial coefficients of p(x), highest degree first.
    params : dict
        Dictionary with keys 'r1','r2','r3','c'
    """

    a = 0.5
    b = 0.5 - L
    if not (0 < b < a):
        raise ValueError("Require 0 < L < 0.5")

    # Default initial guess
    if initial_guess is None:
        initial_guess = np.array([
            b + 0.25*(a-b),
            b + 0.55*(a-b),
            b + 0.85*(a-b)
        ])

    def equations(s):
        s1, s2, s3 = s
        E1 = s1**2 + s2**2 + s3**2
        E2 = s1**2*s2**2 + s1**2*s3**2 + s2**2*s3**2
        E3 = s1**2*s2**2*s3**2

        def H(t):
            return (t**7)/7 - E1*(t**5)/5 + E2*(t**3)/3 - E3*t

        return np.array([
            H(s1) - H(s3),
            H(s2) - H(b),
            H(a)  - H(b)
        ])

    sol = root(equations, initial_guess, method="hybr")
    if not sol.success:
        raise RuntimeError("Nonlinear solve failed: " + sol.message)

    s1, s2, s3 = sol.x

    # Recover r's
    r1, r2, r3 = 0.5 - s1, 0.5 - s2, 0.5 - s3

    # Build polynomial
    E1 = s1**2 + s2**2 + s3**2
    E2 = s1**2*s2**2 + s1**2*s3**2 + s2**2*s3**2
    E3 = s1**2*s2**2*s3**2

    def H(t):
        return (t**7)/7 - E1*(t**5)/5 + E2*(t**3)/3 - E3*t

    c = float(1.0 / (2*H(s1)))

    # p(x) = 1/2 + c * H(x-1/2)
    # Expand into monomials
    y = np.poly1d([1, -0.5])
    Hpoly = (1/7)*y**7 - (E1/5)*y**5 + (E2/3)*y**3 - E3*y
    ppoly = c*Hpoly + 0.5
#    print(ppoly,r1,r2,r3,c)
    
    return ppoly.coeffs, {
        "r1": r1,
        "r2": r2,
        "r3": r3,
        "c":  c
    }

def sp7_acc(X, L_outer, L_inner, H_inner, H_outer, expensive_output=0, fixed_niter=0):
    Csp7 = 64
    qsp7 = 4
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
        idem_err_maxabs.append(sp8py.maxabs(X-X2))
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
        idem_err_maxabs.append(sp8py.maxabs(X-X2))
    while L_outer > precond_limit:
        mc7,params = sp7_polynomial_coeffs(L_outer)
        mc7_prim = np.polyder(mc7)
        print('SP7 slope at 0.5: ',np.polyval(mc7_prim,0.5))
        X = poly7_simple(mc7,X) # FIXME: USE SCHEME USING ONLY 4 MULTS
        polys.append(mc7)
        nmul += 4
#        xx = np.linspace(0,1,100)
#        plt.plot(xx,np.polyval(mc7, xx))
#        plt.show()
        L_inner = np.polyval(mc7, L_inner)
        H_inner = np.polyval(mc7, H_inner)
        L_outer = np.polyval(mc7, L_outer)
        H_outer = np.polyval(mc7, H_outer)
        nmul_vec.append(nmul)
        gap.append(H_inner-L_inner)
        idem_err_trace.append( sp8py.trace_XmX2(X) )
        if expensive_output:
            sp8py.matmul(X, X, X2)
            idem_err_maxabs.append(sp8py.maxabs(X-X2))

    mc7 = [-20., 70., -84., 35., 0., 0., 0., 0.]
    while 1:
        X = poly7_simple(mc7,X) # FIXME: USE SCHEME USING ONLY 4 MULTS
        polys.append(mc7)
        nmul += 4
        L_inner = np.polyval(mc7, L_inner)
        H_inner = np.polyval(mc7, H_inner)
        L_outer = np.polyval(mc7, L_outer)
        H_outer = np.polyval(mc7, H_outer)
        nmul_vec.append(nmul)
        gap.append(H_inner-L_inner)
        idem_err_trace.append( sp8py.trace_XmX2(X) )
        if expensive_output:
            sp8py.matmul(X, X, X2)
            idem_err_maxabs.append(sp8py.maxabs(X-X2))
        if fixed_niter:
            if len(nmul_vec) >= fixed_niter:
                break
            continue
        if idem_err_trace[-1] <= 0:
            break
        if idem_err_trace[-1] > Csp7*idem_err_trace[-2]**qsp7:
            break
    return X,nmul,polys,gap,nmul_vec,idem_err_trace,idem_err_maxabs
