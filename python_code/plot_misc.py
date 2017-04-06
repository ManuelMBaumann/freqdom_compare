import scipy.sparse as sparse
import matplotlib.pyplot as plt
import scipy.io as io
from math import sqrt, atan, cos, sin, pi, atan2
import numpy as np 
from nutils import *

cc  = list('gbcmy')

#def plot_opt_splitting2(e, w, W):
    
    #df     = 0.1
    #w_span = np.arange(w,W+2.0*pi*df,2.0*pi*df)
    #J_span = np.empty((len(w_span),))
    
    #J_span[0] = J(e, w_span[0], w_span[-1])
    #for i in range(1,len(w_span)-1):
        #J_span[i] = max( J(e, w_span[0], w_span[i]), J(e, w_span[i], w_span[-1]) )
    #J_span[-1] = J(e, w_span[0], w_span[-1])
    
    #with plot.PyPlot( 'opt_splitting', figsize=(10,10)) as plt:
        #plt.title('Splitting frequency range')
        #plt.plot(w_span/(2.0*pi), J_span)
        #plt.xlabel('frequency')
        #plt.ylabel(r'$\mathcal{J} (\tau^\ast$)')
        #plt.plot( np.mean([w,W])/(2.0*pi), max(J(e, w_span[0], np.mean([w,W])), J(e, np.mean([w,W]), w_span[-1])), 'rx' )


#def J(e, w, W, tau=False):
    
    #if not tau:
        #tau = opt_tau_anal(e,w,W)
    
    #r     = 0.5*np.sqrt(1.0 + (tau.real/tau.imag)**2)
    #c1_im = tau.real/(2.0*tau.imag) - ((tau.imag+e*tau.real)*w)/((w-tau.real)**2+(e*w+tau.imag)**2)
    #cN_im = tau.real/(2.0*tau.imag) - ((tau.imag+e*tau.real)*W)/((W-tau.real)**2+(e*W+tau.imag)**2)
    #R     = np.sqrt(tau.real**2+tau.imag**2)*np.sqrt((e**2+1.0))/(2.0*abs(tau.real*e+tau.imag))
    #C_im  = e*(tau.real**2+tau.imag**2)/(2.0*tau.imag*(tau.real*e+tau.imag))
    
    ##c1_re = 0.5 - ((1.0+e**2)*w**2+(e*tau.imag-tau.real)*w)/((w-tau.real)**2+(e*w+tau.imag)**2)
    ##cN_re = 0.5 - ((1.0+e**2)*W**2+(e*tau.imag-tau.real)*W)/((W-tau.real)**2+(e*W+tau.imag)**2)
    ##c1    = c1_re+1j*c1_im
    ##cN    = cN_re+1j*cN_im
    ##print (abs(c1)-abs(cN))
    
    #return np.sqrt(r**2/(R**2-C_im**2+2.0*C_im*c1_im))
    

def opt_tau_anal(e,w,W):    
    r        = sqrt(w*W*(1.0+e**2))
    th       = atan(-sqrt( (e**2*(W+w)**2+(W-w)**2) /(4.0*w*W) ))
    #th       = atan(sqrt( (e**2*(W+w)**2+(W-w)**2) /(4.0*w*W) ))
    tau_anal = r*cos(th) + 1j*(r*sin(th))
    #print('DEBUG -- test tau')
    #print( tau_anal.real )
    #print( 2.0*w*W/(w+W) )
    #print( -sqrt((e**2*(w+W)**2+(W-w)**2)*w*W)/(w+W) )
    #print(tau_anal.imag)
    #print('-----------------')
    return tau_anal


#def plot_tau_im(e, w, W, imgtype='png'):
    #col  = list('bgrcmy')
    #mark = list('ovs*Dh')
    #NOP     = 100
    #e_span  = np.linspace(0.0,1.0,NOP)
    #tau_im  = np.empty((NOP,))
    
    #fmin = np.array([1.0, 3.0, 6.0, 9.0])
    #fmax = (W/(2.0*pi))*np.ones((len(fmin),))
       
    #kk = 0
    #my_leg = []
    
    #with plot.PyPlot( 'opt_tau_im', figsize=(10,10), imgtype=imgtype) as plt:
        #for f,F in zip(fmin, fmax):
            #w = 2.0*pi*f
            #W = 2.0*pi*F
            
            #for i in range(NOP):
                #tau       = opt_tau_anal(e_span[i],w,W)
                #tau_im[i] = tau.imag
            
            ##plt.title('Imaginary part of '+r'$\tau^\ast$')
            #plt.plot(e_span, tau_im/W, color=col[kk], linewidth=2.0, marker=mark[kk], markevery=5, markersize=10)
            ##plt.plot(e, opt_tau_anal(e,w,W).imag/W, 'rx', markersize=10, mew=1 )
            #plt.xlabel('damping parameter '+r'$\epsilon$', fontsize=20)
            ##plt.ylabel('relative imaginary part of '+r'$\tau^\ast$')
            #plt.ylabel('imaginary part of '+r'$\tau^\ast$'+' (scaled)', fontsize=20)
            #my_leg = my_leg+['f = ['+str(round(f,1))+','+str(round(F,1))+']']
            #plt.legend(my_leg)
            #kk = kk + 1


#def plot_obj_fun(e, w, W, imgtype='png'):
    #col  = list('bgrcmy')
    #mark = list('ovs*Dh')
    #NOP     = 100
    #e_span  = np.linspace(0.0,1.0,NOP)
    #J_span  = np.empty((NOP,))
    
    #fmin = np.array([1.0, 3.0, 6.0, 9.0])
    #fmax = (W/(2.0*pi))*np.ones((len(fmin),))
    
    #kk = 0
    #my_leg = []
    
    #with plot.PyPlot( 'obj_fun', figsize=(10,10), imgtype=imgtype) as plt:
        #for f,F in zip(fmin, fmax):
            #w = 2.0*pi*f
            #W = 2.0*pi*F
            
            #for i in range(NOP):
                #tau       = opt_tau_anal(e_span[i],w,W)
                #J_span[i] = J(e_span[i], w, W)
                        
            ##plt.title('Objective function value')
            #plt.plot(e_span, J_span, color=col[kk], linewidth=2.0, marker=mark[kk], markevery=5, markersize=10)
            ##plt.plot(e, J(e, w, W), 'rx', markersize=10, mew=1 )
            #plt.xlabel('damping parameter '+r'$\epsilon$', fontsize=20)
            ##plt.ylabel(r'$\mathcal{J} (\tau^\ast$)')
            #plt.ylabel('GMRES-bound in Corrolary 2.7', fontsize=20)
            #my_leg = my_leg+['f = ['+str(round(f,1))+','+str(round(F,1))+']']
            #plt.legend(my_leg)
            #kk = kk + 1
            
        
#def plot_num_tau(om, tau_anal, imgtype='png'):
    
    #my_eps    = 1e-8
    ##dd        = -om[0].imag/om[0].real
    #alpha1    = om.real
    #beta1     = om.imag
    #alpha_max = np.max(alpha1)
    #step1     = 0.007*alpha_max
    #step2     = 0.007*alpha_max
    #alpha2    = np.arange(my_eps, alpha_max+my_eps, step1)
    #beta2     = np.arange(-my_eps, -alpha_max-my_eps, -step2)
    
    #cf    = np.empty((len(alpha2),len(beta2)))
    #c_fac = np.empty((len(alpha2),len(beta2),len(alpha1)))
    
    #for i in range(len(alpha2)):
        #for j in range(len(beta2)):
            #tau = alpha2[i]+1j*beta2[j]
            #for k in range(len(alpha1)):
                
                #omk    = alpha1[k]+1j*beta1[k]
                #eta    = omk/(omk-tau)
                #c_re   = ((0.0 - np.conj(tau))/(tau - np.conj(tau)) - eta).real
                #c_im   = ((0.0 - np.conj(tau))/(tau - np.conj(tau)) - eta).imag
                #radius = abs((tau - 0)/(tau - np.conj(tau)))
        
                #c_fac[i,j,k] = radius/np.sqrt(c_re**2+c_im**2)
            
            #cf[i,j] = np.max(c_fac[i,j,:])
            
    #with plot.PyPlot( 'opt_tau', figsize=(10,10), imgtype=imgtype) as plt:
        #xx, yy = np.meshgrid(alpha2/alpha_max, beta2/alpha_max, sparse=False, indexing='ij')
        #plt.contourf(xx,yy,cf, 20)
        #plt.axis('equal')
        #plt.xlabel('real part (relative)', fontsize=16)
        #plt.ylabel('imag part (relative)', fontsize=16)
        #plt.xlim((0,1))
        #plt.ylim((0,-1))
        #plt.colorbar()
        
        #ind_tau      = np.argmin(cf)
        #i_ind, j_ind = np.unravel_index(ind_tau, cf.shape)
        #tau_num      = alpha2[i_ind]+1j*beta2[j_ind]
        ##plt.plot(tau_num.real/alpha_max, tau_num.imag/alpha_max, linestyle='None', markersize=15, linewidth=3.0, color='y', marker='x', mew=2)
        
        #plt.plot(tau_anal.real/alpha_max, tau_anal.imag/alpha_max, markersize=15, linewidth=3.0, color='w', marker='x', mew=2)
        #NOP = 1000
        ##th  = np.linspace(0.0, atan2(tau_anal.imag,tau_anal.real), NOP)
        #th  = np.linspace(0.0, -pi/2.0, NOP)
        #x_anal = abs(tau_anal/alpha_max)*np.cos(th)
        #y_anal = abs(tau_anal/alpha_max)*np.sin(th)
        #plt.plot(x_anal, y_anal, 'w--')
        ##plt.plot([tau_anal.real/alpha_max,tau_anal.real/alpha_max], [0.0,-1.0], 'w--' )
        ##plt.plot([0.0,tau_anal.real], [0.0,tau_anal.imag], 'w--' )
       
        #plt.plot(om.real/alpha_max, om.imag/alpha_max, linestyle='None', markersize=10, linewidth=3.0, color='k', marker='x', mew=1)
       
    #return tau_num


#def plot_circles_on_circle(A, B, om, tau, plot_spec=False):
    #NOP = 1000
    #th  = np.linspace(0.0,2.0*pi,NOP)
    #Nom = len(om)
    
    #cc  = list('grcmy')
    #col = list('b')
    #j = -1
    #for k in range(1,Nom-1):
        #j=j+1
        #if (j>4):
            #j=0
        #col.append(cc[j])
    #col.append('b')

    #eta = om/(om-tau)
    #dd  = -om[0].imag/om[0].real
    #C   = 0.0 + 1j*( (dd*abs(tau)**2)/(2.0*tau.imag*(tau.imag+dd*tau.real)) )
    #R   = sqrt( abs(tau)**2*(dd**2+1.0)/(4.0*(tau.imag+dd*tau.real)**2) )
    #X   = R*np.cos(th)+C.real
    #Y   = R*np.sin(th)+C.imag
    
    #with plot.PyPlot( 'circles', figsize=(10,10)) as plt:
        #plt.plot(X, Y, 'k')
        #plt.plot(C.real, C.imag, 'kx', markersize=10)
        
        #for k in range(0,Nom):
            #ck = -np.conj(tau)/(tau-np.conj(tau)) - eta[k]
            #r  = abs(tau/(tau-np.conj(tau)))
            #x = r*np.cos(th)+ck.real
            #y = r*np.sin(th)+ck.imag
            #plt.plot(x, y, col[k]+'--')
            #plt.plot(ck.real, ck.imag, col[k]+'x', markersize=10)
            
            #if plot_spec:
                #n = A.shape[0]
                #I = sparse.identity(n).tocsc()
                #P = (A - tau*B).tocsc()
                #Pinv = sparse.linalg.inv(P)
                #vals, vecs  = sparse.linalg.eigs(A.tocsc()*Pinv.tocsc()-eta[k]*I,k=n-2)
                #plt.plot(vals.real, vals.imag, col[k]+'x', markersize=4)
            
        #plt.axhline(linewidth=0.5, color='k')
        #plt.axvline(linewidth=0.5, color='k')
        #plt.axis('equal')

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
    #dd  = -om[0].imag/om[0].real
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
            
#def plot_convergence(resvec):
    #Nom = resvec.shape[1]
    #it  = resvec.shape[0]

    #cc  = list('grcmy')
    #col = list('b')
    #j = -1
    #for k in range(1,Nom-1):
        #j=j+1
        #if (j>4):
            #j=0
        #col.append(cc[j])
    #col.append('b')
    
    #x_as   = np.linspace(0,it,it)
    #my_leg = []
    
    #with plot.PyPlot( 'conv_pmsgmres', figsize=(10,10)) as plt:
        #for k in range(Nom):
            #plt.semilogy(x_as, resvec[:,k]/resvec[0,k],col[k])
            #my_leg = my_leg+['f'+str(k)]
        #plt.title('Convergence of pmsGMRES')       
        #plt.xlabel('Number of matrix-vector multiplications')
        #plt.ylabel('Relative residual norm')
        #plt.ylim((1e-8,1))
        #plt.legend(my_leg)
        #plt.grid()

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