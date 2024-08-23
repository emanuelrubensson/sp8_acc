import numpy as np
def maxabs(A):
    return np.linalg.norm(A.reshape(A.size),ord=np.inf)
