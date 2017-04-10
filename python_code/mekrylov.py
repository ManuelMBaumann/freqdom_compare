import scipy.sparse as sparse
import matplotlib.pyplot as plt
#import scipy.io as io
import numpy as np 
import scipy.sparse.linalg as spla
#import pyamg
from math import sqrt, atan, cos, sin, pi, atan2
from numpy.linalg import norm
#from scipy.io import mmwrite
from nutils import *
from numpy.linalg import solve
from scipy.linalg.blas import get_blas_funcs
from plot_misc import *
import time
import cmath

class convergenceHistory:
    def __init__(self, plot_resnrm=True):
        self.resvec = []
        self.plot_resnrm = plot_resnrm
    def callback(self, _rnrm_):
        self.resvec.append(_rnrm_)
        if self.plot_resnrm:
            print(str(len(self.resvec))+' - '+str(_rnrm_))

class __NoPrecond__(object):
    def solve(self,_X_): return _X_
   

def vectorize_me(A, P=None):
    # simplified for right operators being diag matrices
    Pflag = 0
    if P==None:
        P = __NoPrecond__()
        Pflag = 1
    
    N   = A.K.shape[0]
    Nom = A.Om.shape[0]
    Eij = np.zeros((N,Nom), dtype=complex)
    A_blk = np.zeros((N*Nom,N*Nom), dtype=complex)  
    for i in range(N):
        for j in range(Nom):
            Eij[i,j] = 1.0
            A_blk[j*N:(j+1)*N,i+j*N] = A.dot(P.solve(Eij))[:,j]
            Eij[i,j] = 0.0
            
    with plot.PyPlot( 'blk_eigs', figsize=(10,10)) as plt:
        vals = np.linalg.eigvals(A_blk)
        plt.plot(vals.real, vals.imag, 'bx', markersize=4)
        plt.axhline(linewidth=0.5, color='k')
        plt.axvline(linewidth=0.5, color='k')
        plt.axis('equal')
    if Pflag==1:
        with plot.PyPlot( 'blk_spy', ndigits=0 ) as plt:    
            plt.spy( A_blk, markersize=0.8, precision=0.05)


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


def me_driver(K, C, M, b, freq, tau, damping, tol, maxit, iLU, rot, fill_factor, plot_resnrm):
                         
    class vG_op:
        def __init__(self, K, C, M, Om, P):
            self.K  = K
            self.C  = C
            self.M  = M
            self.Om = Om
            self.P  = P
            self.type = complex
        def dot(self, X):
            X =  self.P.solve(X)
            return self.K.dot(X) + 1j*( self.C.dot( ((self.Om).dot(X.T)).T ) ) - self.M.dot( ((self.Om**2).dot(X.T)).T )
        def resub(self, X):
            return self.P.solve(X)
        
    class precon:
        def __init__(self, K, C, M, tau, timing=False):
            P       = K+1j*tau*C-tau**2*M
            t0 = time.time()
            self.P  = spla.splu(P.tocsc())
            te = time.time()
            if timing:
                print('LU decomposition:'+str(te-t0))
                print('tau = '+str(tau))
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
                print('tau = '+str(tau))
        def solve(self, X):
            return self.P.solve(X) 
        
    class rot_precon:
        def __init__(self, eta, tau):
            c1   = (0-np.conj(tau))/(tau-np.conj(tau)) - eta[0]
            phi1 = cmath.polar(c1)[1]
            rot  = np.ones((len(eta),), dtype=complex)
            for k in range(0,len(eta)):
                ck     = (0-np.conj(tau))/(tau-np.conj(tau)) - eta[k]
                phik   = cmath.polar(ck)[1]
                rot[k] = np.exp(-1j*(phik-phi1))
            self.R  = sparse.diags(rot,0)        
            self.IE = sparse.identity(len(eta)) - sparse.diags(eta,0)
        def solve(self, X):
            return (self.IE.dot(self.R.dot(X.T))).T
        
         
    om  = 2.0*pi*freq*(1.0-1j*damping)
    Om  = sparse.diags(om,0)
    
    if tau.real<0.0:
        tau2 = opt_tau_anal( 2.0*damping/(1.0-damping**2), (1.0-damping**2)*min((2.0*pi*freq)**2), (1.0-damping**2)*max((2.0*pi*freq)**2) )
        tau  = np.sqrt(tau2)
    else:
        tau  = tau*max(om.real)
        tau2 = tau**2
        
    eta = om**2/(om**2-tau2)
    
    # Define preconditioners
    if rot:
        P1 = rot_precon( eta, tau2 )
    else:
        P1 = __NoPrecond__()
    if not iLU:
        P2 = precon( K, C, M, tau, timing=True )
    else:
        P2 = precon_ilu( K, C, M, tau, fill_factor=fill_factor, timing=True )
    
    
    A = vG_op(K, C, M, Om, P1)
    B = (b*np.ones((len(om),1))).T   
    
    #vectorize_me(vG_op(K, C, M, Om, __NoPrecond__()), P2)
    #vectorize_me(A, P2)
    
    # Run global GMRES
    X0      = np.zeros(B.shape, dtype=complex)
    res     = convergenceHistory(plot_resnrm=plot_resnrm)
    X, info = megmres(A, B, X0=X0, tol=tol, maxit=maxit, M1=P2, callback=res.callback)
    X = A.resub(X)
        
    plot_meconvergence(res.resvec)
    I  = sparse.identity(M.shape[0])
    AA = sparse.bmat([[1j*C,K],[I,None]])
    BB = sparse.bmat([[M,None],[None,I]])
    plot_circles_on_circle( AA, BB, om**2, tau**2, damping)
    if rot:
        plot_circles_on_circle( AA, BB, om**2, tau2, damping, rot=P1.R.todense() )
    
    return X.T, len(res.resvec)
