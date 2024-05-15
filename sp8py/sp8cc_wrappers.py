import ctypes
import glob
import os
from numpy.ctypeslib import ndpointer
import numpy as np

lib_search_path = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                               '..','**','sp8cc*.so'))
libfile = glob.glob(lib_search_path, recursive=True)[0]
sp8_cwrappers_lib = ctypes.cdll.LoadLibrary(libfile)

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
        assert(isinstance(B[0][0], float));
        assert(isinstance(C[0][0], float));
        matmul_double_helper(A,B,C,nA)
    if isinstance(A[0][0], np.float32):
        assert(isinstance(B[0][0], np.float32));
        assert(isinstance(C[0][0], np.float32));
        matmul_single_helper(A,B,C,nA)
    return
