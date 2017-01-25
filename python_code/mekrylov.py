import scipy.sparse as sparse
import matplotlib.pyplot as plt
#import scipy.io as io
import numpy as np 
import scipy.sparse.linalg as spla
import pyamg
from math import sqrt, atan, cos, sin, pi, atan2
from numpy.linalg import norm
#from scipy.io import mmwrite
from nutils import *
from numpy.linalg import solve
from scipy.linalg.blas import get_blas_funcs
#from mskrylov import opt_tau_anal
from plot_misc import opt_tau_anal
import time

def plot_rel_convergence(resvec):
    it = len(resvec)
    x_as = np.linspace(0,it,it)
    
    with plot.PyPlot( 'conv_megmres', figsize=(10,10)) as plt:
        plt.semilogy(x_as, resvec[:]/resvec[0])
        plt.title('Convergence of global GMRES')       
        plt.xlabel('Number of operator applications')
        plt.ylabel('Relative residual norm')
        plt.ylim((1e-8,1))
        plt.grid()
        
        
def megmres(A, B, m=200, X0=None, tol=1e-8, maxit=None, M1=None, callback=None):
    size = B.shape
    if maxit is None:
        maxit = 2*np.prod(size)
    if M1 is None:
        # No preconditioner class
        class __NoPrecond__(object):
            def solve(self,_X_): return _X_
        M1 = __NoPrecond__()
    if X0 is None:
        X0 = np.zeros(size, dtype = complex)
    
    X    = np.array(X0)
    bnrm = norm(B)
    info = 1

    # Check for zero rhs:    
    if bnrm == 0.0:
        # Solution is null-vector
        info = 0
        return np.zeros(size), info
    # Compute initial residual:
    R = B - A.dot(X)
    rnrm = norm(R)
    # Relative tolerance
    tolb = tol*bnrm
    if callback is not None:    
        callback(rnrm)

    if rnrm < tolb:
        # Initial guess is a good enough solution
        info = 0
        return X, info
    
    # Initialization
    rotmat = get_blas_funcs('rotg', dtype=np.complex128) # call to ZROTG
    V  = [np.zeros(size, dtype=complex) for i in range(0, m+1)]
    H  = np.zeros((m+1, m), dtype=complex)
    cs = np.zeros(m+1, dtype=np.float64)
    cs_tmp = np.zeros(1, dtype=np.complex128)
    sn = np.zeros(m+1, dtype=np.complex128)
    e1 = np.zeros(m+1, dtype=complex)
    e1[0] = 1.
    for _iter in range(0, maxit):
        # Begin iteration
        V[0] = R/rnrm
        s = rnrm*e1
        for i in range(0, m):
            # Construct orthonormal basis
            # using Gram-Schmidt
            W = A.dot(M1.solve(V[i]))
            
            for k in range(0, i+1):
                H[k, i] = np.vdot(V[k],W)
                W = W - H[k, i]*V[k]
            H[i+1, i] = norm(W)
            V[i+1]    = W/H[i+1, i]
            
            for k in range(0, i):
                # Apply Givens rotation
                temp      = cs[k]*H[k, i] + sn[k]*H[k+1, i]
                H[k+1, i] = -np.conj(sn[k])*H[k, i] + cs[k]*H[k+1, i]
                H[k, i]   = temp
                
            cs_tmp, sn[i] = rotmat(H[i, i], H[i+1, i])
            cs[i]   = cs_tmp.real # BUGFIX: BLAS wrapper out params
            temp    = cs[i]*s[i]
            s[i+1]  = -np.conj(sn[i])*s[i]
            s[i]    = temp
            H[i, i] = cs[i]*H[i, i] + sn[i]*H[i+1, i]
            H[i+1, i] = 0.0
             
            rnrm = abs(s[i+1])
            if callback is not None:    
                callback(rnrm)
             
            if rnrm < tolb:
                y = solve(H[:i, :i],  s[:i])
                Xtmp = np.zeros(size, dtype=complex)
                for k in range(0, i):
                    Xtmp += y[k]*V[k]
                X += M1.solve(Xtmp)
                info = 0
                return X, info

        y = solve(H[:m, :m],  s[:m])
        Xtmp = np.zeros(size, dtype=complex)
        for k in range(0, k):
            Xtmp += y[k]*V[k]
        X += M1.solve(Xtmp)
        R = B - A.dot(X)
        rnrm = norm(R)
        if callback is not None:    
            callback(rnrm)
        if rnrm < tolb:
            info = 0
            break
        
    return X, info



def me_driver(K, C, M, b, freq, tau, damping, tol, maxit, iLU, fill_factor, plot_resnrm):
    
    class vG_op:
        def __init__(self, K, C, M, Om):
            self.K  = K
            self.C  = C
            self.M  = M
            self.Om = Om
            self.type = complex
        def dot(self, X):
            return self.K.dot(X) + 1j*( self.C.dot( ((self.Om).dot(X.T)).T ) ) - self.M.dot( ((self.Om**2).dot(X.T)).T )

    class precon:
        def __init__(self, K, C, M, tau, timing=False):
            P       = K+1j*tau*C-tau**2*M
            t0 = time.time()
            self.P  = spla.splu(P.tocsc())
            te = time.time()
            if timing:
                print('LU decomposition:'+str(te-t0))
        def solve(self, X):
            return self.P.solve(X)
        
    class precon_ilu:
        def __init__(self, K, C, M, tau, fill_factor=1.0, timing=False):
            P      = K+1j*tau*C-tau**2*M
            t0 = time.time()
            self.P = spla.spilu( P.tocsc(), fill_factor=fill_factor)
            te = time.time()
            if timing:
                print('iLU({}) decomposition:'.format(fill_factor)+str(te-t0))
        def solve(self, X):
            return self.P.solve(X) 
        
    class precon_amg:
        def __init__(self, K, C, M, tau, tol=1e-1, timing=False):
            P        = K+1j*tau*C-tau**2*M
            t0 = time.time()
            self.P   = pyamg.smoothed_aggregation_solver(P.tocsr(), max_levels=8, max_coarse=1, keep=True)
            te = time.time()
            self.tol = tol
            if timing:
                print(self.P)
                print('AMG setup:'+str(te-t0))
        def solve(self, X):
            Y = np.empty(X.shape, dtype=complex)
            for k in range(X.shape[1]):
                Y[:,k] = self.P.solve(X[:,k], maxiter=10, cycle='V', accel=None) 
            return Y
        
    class convergenceHistory:
        def __init__(self, plot_resnrm=False):
            self.resvec = []
            self.plot_resnrm = plot_resnrm
        def callback(self, _rnrm_):
            self.resvec.append(_rnrm_)
            if self.plot_resnrm:
                print(str(len(self.resvec))+' - '+str(_rnrm_))


    om  = 2.0*pi*freq*(1.0-1j*damping)
    Om  = sparse.diags(om,0)
    tau = tau*max(om.real)
    if tau.real<0.0:
       tau = opt_tau_anal(0.0,min(om.real),max(om.real))
    B = (b*np.ones((len(om),1))).T   
    A = vG_op(K, C, M, Om)
    
    if not iLU:
        P = precon(K, C, M, tau, timing=True)
    else:
        P = precon_ilu(K, C, M, tau, fill_factor=fill_factor, timing=True)
        #P = precon_amg(K, C, M, tau, tol=1e-1, timing=True)   
    
    X0      = np.zeros(B.shape, dtype=complex)
    bnrm    = np.linalg.norm(B)
    res     = convergenceHistory(plot_resnrm=plot_resnrm)
        
    X, info = megmres(A, B, X0=X0, tol=tol, maxit=maxit, M1=P, callback=res.callback)
    plot_rel_convergence(res.resvec)
       
    return X.T, len(res.resvec)
