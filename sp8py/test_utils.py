import numpy as np
def get_random_X(n,nocc,lumo,homo,float_type = float):
    """
    Constructs a symmetric matrix with 
    all eigenvalues in [0,1],
    n-nocc eigenvalues in [0, lumo] with equidistant spacing,
    nocc eigenvalues in [homo, 1] with equidistant spacing.
    Eigenvectors are taken from Q of a QR-factorization of a random matrix.
    Returns also the projection matrix (density matrix) for projection
    onto the subspace corresponding to the nocc largest eigenvalues. 
    """
    T   = np.random.randn(n,n)
    Q,R = np.linalg.qr(T)
    ev = np.concatenate((np.linspace(0,lumo,n-nocc),
                         np.linspace(homo,1,nocc)))
    X = np.matmul(Q,np.diag(ev))
    X = np.matmul(X,np.transpose(Q))
    X = np.asanyarray(X, order='F', dtype=float_type)
    D = np.matmul(Q[:,-nocc:],np.transpose(Q[:,-nocc:]))
    return X,D
