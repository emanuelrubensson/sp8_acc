import ctypes
from numpy.ctypeslib import ndpointer
sp8_cwrappers_lib = ctypes.cdll.LoadLibrary("../source/sp8_cwrappers.so")

get_sp8_params = sp8_cwrappers_lib.get_sp8_params
get_sp8_params.restype = None
get_sp8_params.argtypes = [ctypes.c_double,
                           ctypes.c_double,
                           ctypes.c_int,
                           ctypes.c_int,
                           ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

get_sp8_monomial_coefficients = sp8_cwrappers_lib.get_sp8_monomial_coefficients
get_sp8_monomial_coefficients.restype = None
get_sp8_monomial_coefficients.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                          ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]
