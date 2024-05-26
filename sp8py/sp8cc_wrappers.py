import ctypes
import glob
import os
from numpy.ctypeslib import ndpointer
import numpy as np

lib_search_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               '..','**','sp8cc*.so'))
libfile = glob.glob(lib_search_path, recursive=True)[0]
sp8_cwrappers_lib = ctypes.cdll.LoadLibrary(libfile)

def none_is_another(lst):
    for ind, el in enumerate(lst):
        for ind2 in range(ind+1,len(lst)):
            if el is lst[ind2]:
                return False
    return True

get_sp8_params_no_acc = sp8_cwrappers_lib.get_sp8_params_no_acc
get_sp8_params_no_acc.restype = None
get_sp8_params_no_acc.argtypes = [ctypes.c_int,
                                  ctypes.c_int,
                                  ndpointer(ctypes.c_double, ndim=1, shape=(9),flags="C_CONTIGUOUS")]

get_sp8_params_helper = sp8_cwrappers_lib.get_sp8_params
get_sp8_params_helper.restype = ctypes.c_int
get_sp8_params_helper.argtypes = [ctypes.c_double,
                                  ctypes.c_double,
                                  ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_bool,
                                  ctypes.c_bool,
                                  ndpointer(ctypes.c_double, ndim=1, shape=(9),flags="C_CONTIGUOUS")] 
def get_sp8_params(lumo,homo,sp8_spec,v):
    left      = sp8_spec[0]
    right     = sp8_spec[1]
    acc_left  = sp8_spec[2]
    acc_right = sp8_spec[3]
    return get_sp8_params_helper(lumo,homo,left,right,acc_left,acc_right,v)
    
get_sp8_params_max_slope = sp8_cwrappers_lib.get_sp8_params_max_slope
get_sp8_params_max_slope.restype = None
get_sp8_params_max_slope.argtypes = [ctypes.c_double,
                                     ctypes.c_double,
                                     ndpointer(ctypes.c_double, ndim=1, shape=(9),flags="C_CONTIGUOUS")]

get_sp8_params_max_gap = sp8_cwrappers_lib.get_sp8_params_max_gap
get_sp8_params_max_gap.restype = None
get_sp8_params_max_gap.argtypes = [ctypes.c_double,
                                   ctypes.c_double,
                                   ctypes.c_double,
                                   ctypes.c_double,
                                   ndpointer(ctypes.c_double, ndim=1, shape=(9),flags="C_CONTIGUOUS")]

get_sp8_monomial_coefficients = sp8_cwrappers_lib.get_sp8_monomial_coefficients
get_sp8_monomial_coefficients.restype = None
get_sp8_monomial_coefficients.argtypes = [ndpointer(ctypes.c_double, ndim=1, shape=(9), flags="C_CONTIGUOUS"),
                                          ndpointer(ctypes.c_double, ndim=1, shape=(9), flags="C_CONTIGUOUS")]

sp8_prim = sp8_cwrappers_lib.sp8_prim
sp8_prim.restype = ctypes.c_double
sp8_prim.argtypes = [ndpointer(ctypes.c_double, ndim=1, shape=(9), flags="C_CONTIGUOUS"),ctypes.c_double]

sp8_eval = sp8_cwrappers_lib.sp8_eval
sp8_eval.restype = ctypes.c_double
sp8_eval.argtypes = [ndpointer(ctypes.c_double, ndim=1, shape=(9), flags="C_CONTIGUOUS"),ctypes.c_double]

matmul_single_helper = sp8_cwrappers_lib.matmul_single
matmul_single_helper.restype = None
matmul_single_helper.argtypes = [ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                 ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                 ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                 ctypes.c_int]

matmul_double_helper = sp8_cwrappers_lib.matmul_double
matmul_double_helper.restype = None
matmul_double_helper.argtypes = [ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                 ctypes.c_int]
def matmul(A,B,C):
    assert(np.ndim(A)==np.ndim(B)==np.ndim(C)==2)
    mA,nA = np.shape(A)
    mB,nB = np.shape(B)
    mC,nC = np.shape(C)    
    assert(mA == nA == mB == nB == mC == nC)
    if isinstance(A[0][0], float): # OBS: "np.float64" is a "float" but not the other way around
        matmul_double_helper(A,B,C,nA)
    if isinstance(A[0][0], np.float32):
        matmul_single_helper(A,B,C,nA)
    return

sastre_poly_8_eval_single_helper = sp8_cwrappers_lib.sastre_poly_8_eval_single
sastre_poly_8_eval_single_helper.restype = None
sastre_poly_8_eval_single_helper.argtypes = [ndpointer(ctypes.c_float, ndim=1, shape=(9), flags="C_CONTIGUOUS"),
                                             ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                             ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                             ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                             ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                             ctypes.c_int]
sastre_poly_8_eval_double_helper = sp8_cwrappers_lib.sastre_poly_8_eval_double
sastre_poly_8_eval_double_helper.restype = None
sastre_poly_8_eval_double_helper.argtypes = [ndpointer(ctypes.c_double, ndim=1, shape=(9), flags="C_CONTIGUOUS"),
                                             ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                             ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                             ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                             ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                             ctypes.c_int]
def sastre_poly_8_eval(mc, X, X2, M2, M3):
    assert(np.ndim(X)==np.ndim(X2)==np.ndim(M2)==np.ndim(M3)==2)
    mX,nX   = np.shape(X)
    mX2,nX2 = np.shape(X2)
    mM2,nM2 = np.shape(M2)
    mM3,nM3 = np.shape(M3)
    assert(mX == nX == mX2 == nX2 == mM2 == nM2 == mM3 == nM3)
    if isinstance(X[0][0], float):
        sastre_poly_8_eval_double_helper(mc, X, X2, M2, M3, nX)
        return
    if isinstance(X[0][0], np.float32):
        mc = np.asanyarray(mc, dtype=np.float32)
        sastre_poly_8_eval_single_helper(mc, X, X2, M2, M3, nX)
        return

poly_8_eval_single_helper = sp8_cwrappers_lib.poly_8_eval_single
poly_8_eval_single_helper.restype = None
poly_8_eval_single_helper.argtypes = [ndpointer(ctypes.c_float, ndim=1, shape=(9), flags="C_CONTIGUOUS"),
                                      ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                      ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                      ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                      ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                      ctypes.c_int]
poly_8_eval_double_helper = sp8_cwrappers_lib.poly_8_eval_double
poly_8_eval_double_helper.restype = None
poly_8_eval_double_helper.argtypes = [ndpointer(ctypes.c_double, ndim=1, shape=(9), flags="C_CONTIGUOUS"),
                                      ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                      ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                      ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                      ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                      ctypes.c_int]
def poly_8_eval(mc, X, X2, M2, M3):
    assert( none_is_another( (X, X2, M2, M3) ) )
    assert(np.ndim(X)==np.ndim(X2)==np.ndim(M2)==np.ndim(M3)==2)
    mX,nX   = np.shape(X)
    mX2,nX2 = np.shape(X2)
    mM2,nM2 = np.shape(M2)
    mM3,nM3 = np.shape(M3)
    assert(mX == nX == mX2 == nX2 == mM2 == nM2 == mM3 == nM3)
    if isinstance(X[0][0], float):
        poly_8_eval_double_helper(mc, X, X2, M2, M3, nX)
        return
    if isinstance(X[0][0], np.float32):
        mc = np.asanyarray(mc, dtype=np.float32)
        poly_8_eval_single_helper(mc, X, X2, M2, M3, nX)
        return

trace_XmX2_single_helper = sp8_cwrappers_lib.trace_XmX2_single
trace_XmX2_single_helper.restype = ctypes.c_float
trace_XmX2_single_helper.argtypes = [ndpointer(ctypes.c_float, ndim=2, flags="F_CONTIGUOUS"),
                                     ctypes.c_int]
trace_XmX2_double_helper = sp8_cwrappers_lib.trace_XmX2_double
trace_XmX2_double_helper.restype = ctypes.c_double
trace_XmX2_double_helper.argtypes = [ndpointer(ctypes.c_double, ndim=2, flags="F_CONTIGUOUS"),
                                     ctypes.c_int]
def trace_XmX2(X):
    assert(np.ndim(X) == 2)
    mX,nX   = np.shape(X)
    assert(mX == nX)
    if isinstance(X[0][0], float):
        return trace_XmX2_double_helper(X,nX)
    if isinstance(X[0][0], np.float32):
        return trace_XmX2_single_helper(X,nX)
