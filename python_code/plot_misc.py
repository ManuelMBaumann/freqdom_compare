import scipy.sparse as sparse
import matplotlib.pyplot as plt
import scipy.io as io
from math import sqrt, atan, cos, sin, pi, atan2
import numpy as np 
from nutils import *

cc  = list('gbcmy')
  
  
def opt_tau_anal(e,w,W):    
    r        = sqrt(w*W*(1.0+e**2))
    th       = atan(-sqrt( (e**2*(W+w)**2+(W-w)**2) /(4.0*w*W) ))
    return r*cos(th) + 1j*(r*sin(th))

def plot_circles_on_circle(A, B, om, tau, dd, plot_spec=False, rot=False):
    NOP = 100
    th  = np.linspace(0.0,2.0*pi,NOP)
    Nom = len(om)
    
    col = list('r')
    j = -1
    for k in range(1,Nom-1):
        j=j+1
        if (j>4):
            j=0
        col.append(cc[j])
    col.append('r')

    eta = om/(om-tau)
    C   = 0.0 + 1j*( (dd*abs(tau)**2)/(2.0*tau.imag*(tau.imag+dd*tau.real)) )
    R   = sqrt( abs(tau)**2*(dd**2+1.0)/(4.0*(tau.imag+dd*tau.real)**2) )
    X   = R*np.cos(th)+C.real
    Y   = R*np.sin(th)+C.imag
    
    with plot.PyPlot( 'circles', figsize=(10,10)) as plt:
        plt.plot(X, Y, 'k')
        plt.plot(C.real, C.imag, 'kx', markersize=10)
        
        for k in range(0,Nom):
            ck = -np.conj(tau)/(tau-np.conj(tau)) - eta[k]
            r  = abs(tau/(tau-np.conj(tau)))
            x = r*np.cos(th)+ck.real
            y = r*np.sin(th)+ck.imag
            if rot is not False:
                tmp = x + 1j*y
                tmp = tmp*rot[k,k]
                ck  = ck*rot[k,k]
                plt.plot(tmp.real, tmp.imag, col[k]+'--')
                plt.plot(ck.real, ck.imag, col[k]+'x', markersize=10)
            else:
                plt.plot(x, y, col[k]+'--')
                plt.plot(ck.real, ck.imag, col[k]+'x', markersize=10)
            
            if plot_spec:
                n = A.shape[0]
                I = sparse.identity(n).tocsc()
                P = (A - tau*B).tocsc()
                Pinv = sparse.linalg.inv(P)
                vals, vecs  = sparse.linalg.eigs(A.tocsc()*Pinv.tocsc()-eta[k]*I,k=n-2)
                plt.plot(vals.real, vals.imag, col[k]+'x', markersize=4)
            
        plt.axhline(linewidth=0.5, color='k')
        plt.axvline(linewidth=0.5, color='k')
        plt.axis('equal')           

def plot_msconvergence(resvec):
    Nom = resvec.shape[1]
    it  = resvec.shape[0]

    col = list('r')
    j = -1
    for k in range(1,Nom-1):
        j=j+1
        if (j>4):
            j=0
        col.append(cc[j])
    col.append('r')
    
    x_as   = np.linspace(0,it,it)
    my_leg = []
    
    with plot.PyPlot( 'conv_pmsgmres', figsize=(10,10)) as plt:
        for k in range(Nom):
            plt.semilogy(x_as, resvec[:,k]/resvec[0,k],col[k])
            my_leg = my_leg+['f'+str(k)]
        plt.title('Convergence of pmsGMRES')       
        plt.xlabel('Number of matrix-vector multiplications')
        plt.ylabel('Relative residual norm')
        plt.ylim((1e-8,1))
        plt.legend(my_leg)
        plt.grid()       
        
def plot_meconvergence(resvec):
    it = len(resvec)
    x_as = np.linspace(0,it,it)
    
    with plot.PyPlot( 'conv_megmres', figsize=(10,10)) as plt:
        plt.semilogy(x_as, resvec[:]/resvec[0])
        plt.title('Convergence of global GMRES')       
        plt.xlabel('Number of operator applications')
        plt.ylabel('Relative residual norm')
        plt.ylim((1e-8,1))
        plt.grid()        