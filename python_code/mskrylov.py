import scipy.sparse as sparse
import matplotlib.pyplot as plt
import scipy.io as io
from math import sqrt, atan, cos, sin, pi, atan2
import numpy as np 
from nutils import *
from plot_misc import *
import scipy.sparse.linalg as spla
import time
import pyamg

class convergenceHistory:
    def __init__(self, plot_resnrm=False):
        self.resvec = []
        self.plot_resnrm = plot_resnrm
    def callback(self, _rnrm_):
        self.resvec.append(_rnrm_)
        if self.plot_resnrm:
            print(str(len(self.resvec))+' - '+str(_rnrm_))

def pascal(n):
    if n==0:
        return [1]
    else:
        N = pascal(n-1)
        return [1] + [N[i] + N[i+1] for i in range(n-1)] + [1]
    
def polypre( dg_pp, eta, tau ):
    Nom    = len(eta)
    center = -np.conj(tau)/(tau-np.conj(tau))
    radius = abs(tau/(tau-np.conj(tau))) 
    omega  = 1.0/center
    #omega  = (1.0/center - 1.0/(2*radius)) 

    # convert to a monic representation:
    coeffs = np.zeros((dg_pp+1,dg_pp+1), dtype=complex)
    for ii in range(dg_pp+1):
        binom_coeffs = pascal(ii)
        for jj in range(ii+1):
            binom_coeffs[jj] = binom_coeffs[jj]*(-1.0)**jj*omega**jj
        coeffs[ii,0:ii+1] = binom_coeffs
    gamma = np.sum(coeffs, axis=0)
    
    # Now compute the parameters of the shifted systems:
    eta_hat = np.zeros((Nom,), dtype=complex)
    gamma_k = np.zeros((dg_pp+1,Nom), dtype=complex)
    for i in range(Nom):
        gamma_k[-1,i] = gamma[-1] 
        for j in range(dg_pp-1,-1,-1):
            gamma_k[j,i] = gamma[j] + eta[i]*gamma_k[j+1,i]
        eta_hat[i] = eta[i]*gamma_k[0,i]
        
    return gamma, gamma_k, eta_hat


def msgmres( A, Pl, Pr, b, omega, tau, tol, maxit, callback=None ):
    
    n           = len(b)
    Nom         = len(omega)
    x           = np.zeros((n,Nom), dtype=complex)
    resnrm      = np.zeros((maxit+1,Nom))
    V           = np.zeros((n,maxit+1), dtype=complex)
    H           = np.zeros((maxit+1,maxit), dtype=complex)
    
    beta        = np.linalg.norm(b)
    resnrm[0,:] = beta
    v           = (1.0/beta)*b
    V[:,0]      = v
    
    for k in range(maxit):
        
        v = Pl.solve(A.dot(Pr.solve(v)))
            
        # Gram Schmidt
        for j in range(k+1):
            H[j,k] = np.vdot(V[:,j],v)
            v      = v - H[j,k]*V[:,j]

        H[k+1,k] = np.linalg.norm(v)
        v        = (1.0/H[k+1,k])*v
        V[:,k+1] = v
        
        # Solve shifted Hessenberg systems
        Hr    = H[:k+2,:k+1]       
        e1    = np.zeros(k+2)
        e1[0] = beta 
        Ir    = np.eye(k+2,k+1)
        y     = np.zeros((k+1,Nom), dtype=complex)
        for ii in range(Nom):
            Hr_s           = Hr-omega[ii]*Ir           
            y[:,ii]        = np.linalg.lstsq(Hr_s,e1)[0] # y = min || Hr_s*y - e1 ||
            resnrm[k+1,ii] = np.linalg.norm(e1 - np.dot(Hr_s,y[:,ii]))  
            
        normr = np.max(resnrm[k+1,:])
        if callback is not None:
            callback(normr)
        if normr<tol*beta:
            break
        
    x      = np.dot(V[:,:k+1],y) 
    it     = k+1
    resnrm = resnrm[:it+1,:]
        
    return x, it, resnrm

 
def poly_driver( K, C, M, rhs, freq, tau, damping, tol, maxit, dg_pp, iLU, fill_factor, plot_resnrm ):
    
    Nom = len(freq)
    n   = len(rhs)
    trace_C = np.sum(C.diagonal())
    
    if trace_C>0:
        
        class make_sys:
            def __init__(self, K, C, M):
                n = K.shape[0]
                I = sparse.identity(n)
                A = sparse.bmat([[1j*C,K],[I,None]])
                self.A = A
                #self.K = K
                #self.C = C
                #self.M = M
                self.type = complex
            def dot(self, x):
                #n = len(x)
                #return np.concatenate( (1j*self.C*x[:n//2]+self.K*x[n//2:],x[:n//2]) )
                return self.A*x
            
        class shift_precon:
            def __init__(self, K, C, M, tau, eta, timing=True):
                self.C = C
                self.M = M
                self.tau = tau
                self.eta = eta
                P = K+1j*tau*C-tau**2*M
                t0 = time.time()
                self.P  = spla.splu(P.tocsc())
                te = time.time()
                if timing:
                    print('LU decomposition:'+str(te-t0)) 
                    print('tau = '+str(tau))
            def solve(self, x):
                n  = len(x)
                x1 = x[n//2:]
                x2 = x[:n//2] + (-1j*self.C*x[n//2:] + self.tau*self.M*x[n//2:])
                x2 = self.P.solve(x2)
                x1 = x1 + self.tau*x2
                return np.concatenate((x1,x2))
            def resub(self, X):
                N   = X.shape[0]
                Nom = X.shape[1]
                for i in range(Nom):
                    X[:,i] = self.solve(X[:,i])
                    X[:,i] = (1.0-self.eta[i])*X[:,i]
                return X[N//2:,:]
            
        class shift_precon_ilu:
            def __init__(self, K, C, M, tau, eta, fill_factor=1.0, timing=True):
                self.C = C
                self.M = M
                self.tau = tau
                self.eta = eta
                self.P = K+1j*tau*C-tau**2*M
                t0 = time.time()
                Pinv  = spla.spilu(self.P.tocsc(), fill_factor=fill_factor)
                Pinv_x = lambda x: Pinv.solve(x)
                self.Pinv = spla.LinearOperator(self.P.shape, Pinv_x)
                te = time.time()
                if timing:
                    print('iLU({}) decomposition:'.format(fill_factor)+str(te-t0))
                    print('tau = '+str(tau))
            def solve(self, x):
                n  = len(x)
                x1 = x[n//2:]
                x2 = x[:n//2] + (-1j*self.C*x[n//2:] + self.tau*self.M*x[n//2:])
                x2, info = spla.gmres(self.P, x2, tol=1e-16, restart=200, maxiter=600, M=self.Pinv)
                assert(info==0) # apply shift-and-invert exactly
                x1 = x1 + self.tau*x2
                return np.concatenate((x1,x2))
            def resub(self, X):
                N   = X.shape[0]
                Nom = X.shape[1]
                for i in range(Nom):
                    X[:,i] = self.solve(X[:,i])
                    X[:,i] = (1.0-self.eta[i])*X[:,i]
                return X[N//2:,:]
            
        #class shift_precon_amg:
            #def __init__(self, K, C, M, tau, tol=1e-1, timing=True):
                #self.C = C
                #self.M = M
                #self.tau = tau
                #self.eta = eta
                #self.P = K+1j*tau*C-tau**2*M
                #t0 = time.time()
                #self.settings = {"max_levels": 8,
                                 #"max_coarse": 1,
                                 #"keep": False,
                                 #"presmoother": ("gauss_seidel",
                                         #{"sweep": "symmetric"}),
                                  #"postsmoother": ("gauss_seidel",
                                         #{"sweep": "symmetric"})}
                #Pinv = pyamg.smoothed_aggregation_solver(self.P.tocsr(), **self.settings)
                #Pinv_x = lambda x: Pinv.solve(x, tol=1e-2, maxiter=100, cycle='V', accel=None) 
                #self.Pinv = spla.LinearOperator(self.P.shape, Pinv_x)
                #te = time.time()
                #if timing:
                    #print(Pinv)
                    #print('AMG setup:'+str(te-t0))
            #def solve(self, x):
                #n  = len(x)
                #x1 = x[n//2:]
                #x2 = x[:n//2] + (-1j*self.C*x[n//2:] + self.tau*self.M*x[n//2:])
                #x2, info = spla.gmres(self.P, x2, tol=1e-16, restart=200, maxiter=600, M=self.Pinv)
                #assert(info==0) # apply shift-and-invert exactly
                #x1 = x1 + self.tau*x2
                #return np.concatenate((x1,x2))
            #def resub(self, X):
                #N   = X.shape[0]
                #Nom = X.shape[1]
                #for i in range(Nom):
                    #X[:,i] = self.solve(X[:,i])
                    #X[:,i] = (1.0-self.eta[i])*X[:,i]
                #return X[N//2:,:]

        class poly_precon:
            def __init__(self, A, P, dg_pp, tau, eta):
                self.A = A
                self.P = P
                self.dg_pp = dg_pp
                if dg_pp>0:
                    self.gamma, self.gamma_k, self.eta_k = polypre( dg_pp, eta, tau )
                else:
                    self.eta_k = eta
            def solve(self, x):
                if self.dg_pp>0:
                    n = len(self.gamma)
                    y = self.gamma[-1]*x
                    for i in range(n-2,-1,-1):
                        y = self.A.dot(self.P.solve(y)) + self.gamma[i]*x
                    return y
                else:
                    return x
            def resub(self, X):
                n   = X.shape[0]
                Nom = X.shape[1]
                Y   = np.zeros(X.shape, dtype=complex)
                if self.dg_pp>0:
                    for j in range(Nom):
                        Y[:,j] = self.gamma_k[-1,j]*X[:,j]
                        for i in range(self.dg_pp-1,-1,-1):
                            Y[:,j] = self.A.dot(self.P.solve(Y[:,j])) + self.gamma_k[i,j]*X[:,j]
                    return Y
                else:
                    return X
                
                
        om  = 2.0*pi*freq*(1.0-1j*damping)
        tau = tau*max(om.real)
        if tau.real<0.0:
            tau = opt_tau_anal(damping, min(om.real), max(om.real))
        eta = om/(om-tau) 
        
        A  = make_sys( K, C, M )
        if not iLU:
            Pr = shift_precon( K, C, M, tau, eta )
        else:
            Pr = shift_precon_ilu( K, C, M, tau, eta, fill_factor=fill_factor )
            #Pr = shift_precon_amg( K, C, M, tau, eta )
            
        Pl = poly_precon( A, Pr, dg_pp, tau, eta )
        b  = np.concatenate( (rhs,np.zeros((n,))) )

    else:
        
        class make_sys:
            def __init__(self, K, C, M):
                self.A = K
                self.type = complex
            def dot(self, x):
                return self.A*x
            
        class shift_precon:
            def __init__(self, K, C, M, tau, eta, timing=True):
                self.eta = eta
                P       = K-tau*M
                t0 = time.time()
                self.P  = spla.splu(P.tocsc())
                te = time.time()
                if timing:
                    print('LU decomposition:'+str(te-t0))
                    print('tau = '+str(tau))
            def solve(self, x):
                return self.P.solve(x)
            def resub(self, X):
                Nom = X.shape[1]
                X = self.solve(X)
                for i in range(Nom):
                    X[:,i] = (1.0-self.eta[i])*X[:,i]
                return X
            
        class shift_precon_ilu:
            def __init__(self, K, C, M, tau, eta, fill_factor=1.0, timing=True):
                self.eta = eta
                self.P = K-tau*M
                t0 = time.time()
                Pinv   = spla.spilu(self.P.tocsc(), fill_factor=fill_factor)
                Pinv_x = lambda x: Pinv.solve(x)
                self.Pinv = spla.LinearOperator(self.P.shape, Pinv_x)
                te = time.time()
                if timing:
                    print('iLU({}) decomposition:'.format(fill_factor)+str(te-t0))  
                    print('tau = '+str(tau))
            def solve(self, x):
                y, info = spla.gmres(self.P, x, tol=1e-16, restart=200, maxiter=600, M=self.Pinv)
                assert(info==0) # apply shift-and-invert exactly
                return y
            def resub(self, X):
                Nom = X.shape[1]
                X = self.solve(X)
                for i in range(Nom):
                    X[:,i] = (1.0-self.eta[i])*X[:,i]
                return X
            
        class poly_precon:
            def __init__(self, A, P, dg_pp, tau, eta):
                self.A = A
                self.P = P
                self.dg_pp = dg_pp
                if dg_pp>0:
                    self.gamma, self.gamma_k, self.eta_k = polypre( dg_pp, eta, tau )
                else:
                    self.eta_k = eta
            def solve(self, x):
                if self.dg_pp>0:
                    n = len(self.gamma)
                    y = self.gamma[-1]*x
                    for i in range(n-2,-1,-1):
                        y = self.A.dot(self.P.solve(y)) + self.gamma[i]*x
                    return y
                else:
                    return x
            def resub(self, X):
                Nom = X.shape[1]
                Y   = np.zeros(X.shape, dtype=complex)
                if self.dg_pp>0:
                    for j in range(Nom):
                        Y[:,j] = self.gamma_k[-1,j]*X[:,j]
                        for i in range(self.dg_pp-1,-1,-1):
                            Y[:,j] = self.A.dot(self.P.solve(Y[:,j])) + self.gamma_k[i,j]*X[:,j]
                    return Y
                else:
                    return X
            
        om  = (2.0*pi*freq*(1.0-1j*damping))**2
        tau = tau*max(om.real)
        if tau.real<0.0:
            #tau = opt_tau_anal(damping, min(om.real), max(om.real) )
            tau = opt_tau_anal(2.0*damping/(1.0-damping**2), min(om.real), max(om.real))
        eta = om/(om-tau)
        
        A  = make_sys( K, C, M )
        if not iLU:
            Pr = shift_precon( K, C, M, tau, eta )
        else:
            Pr = shift_precon_ilu( K, C, M, tau, eta, fill_factor=fill_factor )
        Pl = poly_precon( A, Pr, dg_pp, tau, eta )
        b  = rhs
        

    
    res     = convergenceHistory(plot_resnrm=plot_resnrm)
    y, it, resnrm = msgmres( A, Pl, Pr, b, Pl.eta_k, tau, tol, maxit, callback=res.callback )
    x = Pr.resub( Pl.resub(y) )
    plot_msconvergence(resnrm)
    
    return x.T, it
  
  
def msfom( A, P, b, omega, tau, tol, maxit ):   
    n           = len(b)
    Nom         = len(omega)
    collfac     = np.zeros((Nom,), dtype=complex) # co-linearity factors
    x           = np.zeros((n,Nom), dtype=complex)
    V           = np.zeros((n,maxit+1), dtype=complex)
    H           = np.zeros((maxit+1,maxit), dtype=complex)
    
    beta        = np.linalg.norm(b)
    v           = (1.0/beta)*b
    V[:,0]      = v
        
    for k in range(maxit):
        
        v = A.dot(P.solve(v))
            
        # Gram Schmidt
        for j in range(k+1):
            H[j,k] = np.vdot(V[:,j],v)
            v      = v - H[j,k]*V[:,j]

        H[k+1,k] = np.linalg.norm(v)
        v        = (1.0/H[k+1,k])*v
        V[:,k+1] = v
        
        # Solve shifted Hessenberg systems
        Hk    = H[:k+1,:k+1]       
        e1    = np.zeros(k+1)
        e1[0] = beta 
        Ik    = np.eye(k+1)
        y     = np.zeros((k+1,Nom), dtype=complex)
        for ii in range(Nom-1,-1,-1):
            Hr_s       = Hk-omega[ii]*Ik
            y[:,ii]    = np.linalg.solve(Hr_s,e1)
            collfac[ii] = y[-1,ii]/y[-1,-1]
            
        normr = np.max( abs(collfac)*np.linalg.norm(-H[k+1,k]*y[-1,-1]*v) )
        if normr<tol*beta:
            break
        
    x      = np.dot(V[:,:k+1],y) 
    it     = k+1
        
    return x, collfac, it


def fmsgmres( A, P, b, omega, tau, tol_i, tol_o, maxit_i, maxit_o, callback=None ):
    n        = len(b)
    Nom      = len(omega)
    x        = np.zeros((n,Nom), dtype=complex)
    resnrm   = np.zeros((maxit_o+1,Nom))
    V        = np.zeros((n,maxit_o+1), dtype=complex)
    H        = np.zeros((maxit_o+1,maxit_o), dtype=complex)
    Z        = np.zeros((n,maxit_o,Nom), dtype=complex)
    collfacs = np.zeros((maxit_o,Nom), dtype=complex)
    
    beta        = np.linalg.norm(b)
    resnrm[0,:] = beta
    v           = (1.0/beta)*b
    V[:,0]      = v
    
    for k in range(maxit_o):
        
        Z[:,k,:], collfacs[k,:], it_i = msfom(A, P, V[:,k], omega, tau, tol_i, maxit_i)
        print( 'Outer #'+str(k)+' No. inner iters:'+str(it_i) )
        v = A.dot(P.solve(Z[:,k,-1]))
            
        # Gram-Schmidt
        for j in range(k+1):
            H[j,k] = np.vdot(V[:,j],v)
            v      = v - H[j,k]*V[:,j]
        H[k+1,k] = np.linalg.norm(v)
        v        = (1.0/H[k+1,k])*v
        V[:,k+1] = v
        
        # Solve shifted Hessenberg systems
        Hk    = H[:k+2,:k+1]       
        e1    = np.zeros(k+2)
        e1[0] = beta 
        Ik    = np.eye(k+2,k+1)
        y     = np.zeros((k+1,Nom), dtype=complex)
        for ii in range(Nom):
            Gamma          = np.diag(collfacs[:k+1,ii])
            Hk_s           = np.dot(Hk-Ik,Gamma)+Ik
            y[:,ii]        = np.linalg.lstsq(Hk_s,e1)[0] # y = min ||Hr_s*y - e1||
            resnrm[k+1,ii] = np.linalg.norm(e1 - np.dot(Hk_s,y[:,ii]))  
         
        normr = np.max(resnrm[k+1,:])
        if callback is not None:
            callback(normr)
        if normr<tol_o*beta:
            break
    
    for ii in range(Nom):
        x[:,ii] = np.dot(Z[:,:k+1,ii],y[:,ii]) 
    it     = k+1
    resnrm = resnrm[:it+1,:]

    return x, it, resnrm


def nested_driver( K, C, M, rhs, freq, tau, damping, tol_i, tol_o, maxit_i, maxit_o, iLU, fill_factor, plot_resnrm ):
    
    Nom = len(freq)
    n   = len(rhs)
    trace_C = np.sum(C.diagonal())
    
    if trace_C>0:
        
        class make_sys:
            def __init__(self, K, C, M, base):
                n = K.shape[0]
                I = sparse.identity(n)
                A = sparse.bmat([[1j*C,K],[I,None]])
                B = sparse.bmat([[M,None],[None,I]])
                self.A = A-base*B
                self.type = complex
            def dot(self, x):
                return self.A*x
            
        class shift_precon:
            def __init__(self, K, C, M, tau, eta, timing=True):
                self.C = C
                self.M = M
                self.tau = tau
                self.eta = eta
                P       = K+1j*tau*C-tau**2*M
                t0 = time.time()
                self.P  = spla.splu(P.tocsc())
                te = time.time()
                if timing:
                    print('LU decomposition:'+str(te-t0))
                    print('tau = '+str(tau))
            def solve(self, x):
                n  = len(x)
                x1 = x[n//2:]
                x2 = x[:n//2] + (-1j*self.C*x[n//2:] + self.tau*self.M*x[n//2:])
                x2 = self.P.solve(x2)
                x1 = x1 + self.tau*x2
                return np.concatenate((x1,x2))
            def resub(self, X):
                n   = X.shape[0]
                Nom = X.shape[1]
                for i in range(Nom):
                    X[:,i] = self.solve(X[:,i])
                    X[:,i] = (1.0-self.eta[i])*X[:,i]
                return X[n//2:,:]
            
        class shift_precon_ilu:
            def __init__(self, K, C, M, tau, eta, fill_factor=1.0, timing=True):
                self.C = C
                self.M = M
                self.tau = tau
                self.eta = eta
                self.P = K+1j*tau*C-tau**2*M
                t0 = time.time()
                Pinv   = spla.spilu(self.P.tocsc(), fill_factor=fill_factor)
                Pinv_x = lambda x: Pinv.solve(x)
                self.Pinv = spla.LinearOperator(self.P.shape, Pinv_x)
                te = time.time()
                if timing:
                    print('iLU({}) decomposition:'.format(fill_factor)+str(te-t0))   
                    print('tau = '+str(tau))
            def solve(self, x):
                n  = len(x)
                x1 = x[n//2:]
                x2 = x[:n//2] + (-1j*self.C*x[n//2:] + self.tau*self.M*x[n//2:])
                x2, info = spla.gmres(self.P, x2, tol=1e-16, restart=200, maxiter=600, M=self.Pinv)
                assert(info==0) # apply shift-and-invert exactly
                x1 = x1 + self.tau*x2
                return np.concatenate((x1,x2))
            def resub(self, X):
                n   = X.shape[0]
                Nom = X.shape[1]
                for i in range(Nom):
                    X[:,i] = self.solve(X[:,i])
                    X[:,i] = (1.0-self.eta[i])*X[:,i]
                return X[n//2:,:]
            
       
        om   = 2.0*pi*freq*(1.0-1j*damping)
        base = om[-1]
        tau = tau*max(om.real)
        if tau.real<0.0:
            tau = opt_tau_anal(damping, min(om.real), max(om.real) )
        om   = om - base
        eta  = om/(om-(tau-base))
        
        A = make_sys( K, C, M, base )
        if not iLU:
            P = shift_precon( K, C, M, tau, eta )
        else:
            P = shift_precon_ilu( K, C, M, tau, eta, fill_factor=fill_factor )
        b = np.concatenate( (rhs,np.zeros((n,))) )
        
    else:
        
        class make_sys:
            def __init__(self, K, C, M, base):
                self.A = K-base*M
                self.type = complex
            def dot(self, x):
                return self.A*x
            
        class shift_precon:
            def __init__(self, K, C, M, tau, eta, timing=True):
                self.eta = eta
                P       = K-tau*M
                t0 = time.time()
                self.P  = spla.splu(P.tocsc())
                te = time.time()
                if timing:
                    print('LU decomposition:'+str(te-t0))
                    print('tau = '+str(tau))
            def solve(self, x):
                return self.P.solve(x)
            def resub(self, X):
                Nom = X.shape[1]
                X = self.solve(X)
                for i in range(Nom):
                    X[:,i] = (1.0-self.eta[i])*X[:,i]
                return X
            
        class shift_precon_ilu:
            def __init__(self, K, C, M, tau, eta, fill_factor=1.0, timing=True):
                self.eta = eta
                self.P = K-tau*M
                t0 = time.time()
                Pinv   = spla.spilu(self.P.tocsc(), fill_factor=fill_factor)
                Pinv_x = lambda x: Pinv.solve(x)
                self.Pinv = spla.LinearOperator(self.P.shape, Pinv_x)
                te = time.time()
                if timing:
                    print('iLU({}) decomposition:'.format(fill_factor)+str(te-t0))   
                    print('tau = '+str(tau))
            def solve(self, x):
                y, info = spla.gmres(self.P, x, tol=1e-16, restart=200, maxiter=600, M=self.Pinv)
                assert(info==0) # apply shift-and-invert exactly
                return y
            def resub(self, X):
                Nom = X.shape[1]
                X = self.solve(X)
                for i in range(Nom):
                    X[:,i] = (1.0-self.eta[i])*X[:,i]
                return X
        
        
        om   = (2.0*pi*freq*(1.0-1j*damping))**2
        base = om[-1]
        tau  = tau*max(om.real)
        if tau.real<0.0:
            #tau = opt_tau_anal(damping, min(om.real), max(om.real) )
            tau = opt_tau_anal(2.0*damping/(1.0-damping**2), min(om.real), max(om.real))
        om  = om - base
        eta = om/(om-(tau-base))
        
        A  = make_sys( K, C, M, base )
        if not iLU:
            P  = shift_precon( K, C, M, tau, eta )
        else:
            P  = shift_precon_ilu( K, C, M, tau, eta, fill_factor=fill_factor )
        b  = rhs
       
       
      
    res     = convergenceHistory(plot_resnrm=plot_resnrm)
    y, it, resnrm = fmsgmres( A, P, b, eta, tau, tol_i, tol_o, maxit_i, maxit_o, callback=res.callback )
    x = P.resub( y )

    plot_msconvergence(resnrm)

    return x.T, it