import ctypes
import glob
import os
from numpy.ctypeslib import ndpointer

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

sp8 = sp8_cwrappers_lib.sp8_eval
sp8.restype = ctypes.c_double
sp8.argtypes = [ndpointer(ctypes.c_double, ndim=1, shape=(9), flags="C_CONTIGUOUS"),ctypes.c_double]
