from __future__ import print_function, division

import matplotlib
matplotlib.use('agg')

from nutils import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.io import mmwrite
from mskrylov import poly_driver, nested_driver
from mekrylov import me_driver
import time
from math import pi
import pyamg
import scipy.sparse.linalg as spla

@log.title
def makeplots( domain, geom, Lx, Lz, value, name, title, ndigits=0, index=None, clim=None, lineOn=False, imgtype=None,):
  points, colors = domain.elem_eval( [ geom, value ], ischeme='bezier3', separate=True )

  with plot.PyPlot( name, ndigits=ndigits, figsize=(10,10), index=index, imgtype=imgtype ) as plt:
    plt.mesh( points, colors, triangulate='bezier', edgecolors='none' )
    if imgtype is not 'eps':
        plt.title(title, fontsize=20)
    plt.xlabel('x [m]', fontsize=20)
    plt.ylabel('z [m]', fontsize=20)
    plt.xticks( [0, 0.5*Lx, Lx], ['0', str(int(0.5*Lx)), str(int(Lx))], fontsize=20 )
    plt.yticks( [-0, -0.5*Lz, -Lz], ['0', str(int(0.5*Lz)), str(int(Lz))], fontsize=20 )

    if clim is not None:
      plt.clim(*clim)
    if imgtype is not 'eps':
        plt.colorbar()

def makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, name):
  Nom = len(freq)
  vtk_geom, vtk_rho, vtk_lam, vtk_mu, vtk_cp, vtk_cs = domain.simplex.elem_eval( [ geom, rho, lam, mu, cp, cs ], ischeme='vtk', separate=True )
  with plot.VTKFile( name ) as vtk:
      vtk.unstructuredgrid( vtk_geom )
      vtk.pointdataarray( 'rho', vtk_rho )
      vtk.pointdataarray( 'lambda', vtk_lam )
      vtk.pointdataarray( 'mu', vtk_mu )
      vtk.pointdataarray( 'cp', vtk_cp )
      vtk.pointdataarray( 'cs', vtk_cs )
      for k in range(0,Nom):
          disp = vec_basis.dot( sol[k,:] ).real
          vtk_disp = domain.simplex.elem_eval( disp, ischeme='vtk', separate=True )
          vtk.pointdataarray( 'disp_f'+str(freq[k]), vtk_disp )

def makespyplot( matrix, name, imgtype=None ):
  if not scipy.sparse.isspmatrix( matrix ):
      matrix = matrix.toscipy()
  with plot.PyPlot( name, ndigits=0, imgtype=imgtype ) as plt:
    plt.spy( matrix, markersize=0.8, color='black')
    plt.title( name+', nnz = '+str(matrix.nnz) )

def point_eval(func, domain, geom, point):
  domain = domain[tuple(slice(0, p) if p > 0 else slice(None) for p in point)]
  for p in point:
      domain = domain.boundary['right' if p > 0 else 'left']
  return numpy.asarray(domain.integrate( func, geometry=geom, ischeme='gauss2' ).toscipy().todense())

def elast_mat(rho, cp, cs, lam, mu, ndims, nx, ny, nz, vec_basis, domain, geom, block):
  # define PDE
  stress = lambda u: lam*u.div(geom)[:,_,_]*function.eye(ndims) + 2.0*mu*u.symgrad(geom)
  elasticity = function.outer( stress(vec_basis), vec_basis.grad(geom) ).sum([2,3])

  w_mass = lambda u: rho*u
  mass = function.outer( w_mass(vec_basis), vec_basis ).sum(-1)

  # define BC
  n = geom.normal()
  t = np.eye(ndims)
  t = t-(t*n[_,:]).sum(1)
  B_bc = cp*n[:,_]*n[_,:]+cs*(t[:,:,_]*t[:,_,:]).sum(0)

  bc_fun = lambda u: rho*(B_bc*u[:,_,:]).sum(-1)
  sommerfeld = function.outer( bc_fun(vec_basis), vec_basis ).sum(-1)

  if ndims == 2:
      sommerfeld_boundary = 'left,right,bottom'
      source_position = nx//2, nz
  else:
      sommerfeld_boundary = 'left,right,bottom,front,back'
      source_position = nx//2, ny//2, nz

  # Build matrices
  K, M  = domain.integrate( [elasticity, mass], geometry=geom, ischeme='gauss2' )
  C     = domain.boundary[sommerfeld_boundary].integrate( sommerfeld, geometry=geom, ischeme='gauss2' )
  
  # Build RHS
  if not block:
      C = 0*C
      source_position = int(np.round(6.0/10.0*nx)), int(np.round(4.0/10.0*nz))
      rhs = point_eval(vec_basis, domain, geom, source_position)[:,-1] #+ point_eval(vec_basis, domain, geom, source_position)[:,0] 
  else:
      rhs = point_eval(vec_basis, domain, geom, source_position)[:,-1] 
             
  return K, C, M, rhs


def test_orig_problem(K, C, M, b, om, x):
    Nom    = len(om)
    normb  = np.linalg.norm(b)
    relerr = np.zeros((Nom,))
    for j in range(Nom):
        xs = x[:,j]
        r  = b - (K*xs + 1j*om[j]*(C*xs) - om[j]**2*(M*xs))
        relerr[j] = np.linalg.norm(r)/normb
    print('Relative residual of original problem:' +str(relerr))


def main( ndims=2,           # problem dimension (2,3) 
          dx=100.0,          # grid size in x-direction 
          dy=100.0,          # grid size in y-direction          
          dz=100.0,          # grid size in z-direction  
          freq=[1.0,9.0],    # frequencies in Hz 
          Nom=5,             # number of equally-spaced freq's
          degree=1,          # degree of FEM splines
          damping=0.7,       # viscous damping param    
          maxit=300,         # is also maxit_o in nested algo
          maxit_i=10,        # maxit of inner method
          tol=1e-8,          # is also tol_o in nested algo
          tol_i=1e-1,        # tol of inner method
          dg_pp=0,           # degree of poly preconditioner
          tau_re=0.7,        # real(seed), if tau.real<0: take 'optimal' tau
          tau_im=-0.3,       # imag(seed)
          iLU=False,         # exact or inexact shif-and-invert
          rot =True,         # spactral rotation for MatrEqn
          fill_factor=1.0,   # fill-in in iLU decomposition
          block=True,        # C=0 if False
          plots=False,       # plots on/off
          plot_resnrm=False, # display residual norm live
          solver_flag=-1):   # -1(python's built-in), 0(MatrEqn), 1(poly_pre), 2(nested)


  # domain size
  Lx = 500.0
  Ly = 500.0
  Lz = 500.0
  
  tau = tau_re+1j*tau_im

  # problem parameters
  freq = np.linspace(freq[0],freq[-1],Nom)
  om   = 2.0*np.pi*freq*(1.0-1j*damping)
  Nom  = len(om)

  # define physical params
  rho0 = 1800.0
  cp0  = 2000.0
  cs0  = 800.0
  
  rho1 = 2100.0
  cp1  = 3000.0
  cs1  = 1600.0

  # define Cartesian grid
  nx = int(np.round(Lx/dx))+1
  nz = int(np.round(Lz/dz))+1
  verts_x = np.linspace( 0, Lx, nx )
  verts_z = np.linspace( -Lz, 0, nz )

  if ndims == 2:
      ny = 1
      dy = 0.
      verts = [verts_x, verts_z]
  elif ndims == 3:
      ny = int(np.round(Ly/dy))+1
      verts_y = np.linspace( 0, Ly, ny )
      verts = [verts_x, verts_y, verts_z]

  domain, geom = mesh.rectilinear(verts)
  vec_basis    = domain.splinefunc( degree=degree ).vector( ndims )
  
  # define block-in-block problem  
  rho = function.select( [function.greater(0.35*Lx,geom[0]), function.greater(-0.85*Lx,-geom[0]), 
                          function.greater(0.35*Lz,-geom[-1]), function.greater(-0.85*Lz,geom[-1])],
                          [rho0, rho0, rho0, rho0], rho1)
  cp  = function.select( [function.greater(0.35*Lx,geom[0]), function.greater(-0.85*Lx,-geom[0]), 
                          function.greater(0.35*Lz,-geom[-1]), function.greater(-0.85*Lz,geom[-1])],
                          [cp0, cp0, cp0, cp0], cp1)
  cs  = function.select( [function.greater(0.35*Lx,geom[0]), function.greater(-0.85*Lx,-geom[0]), 
                          function.greater(0.35*Lz,-geom[-1]), function.greater(-0.85*Lz,geom[-1])],
                          [cs0, cs0, cs0, cs0], cs1)
  mu   = cs**2 * rho
  lam  = rho * (cp**2 - 2.0*cs**2) 
    
  # problem summary
  print( '----     SQUARES PROBLEM     ----' )
  ppw = 20.0  
  print( 'problem size   : ' + str(nx-1+degree)+' x '+str(ny-1+degree)+' x '+str(nz-1+degree) )
  print( '# dofs         : ' + str(len(vec_basis)) )
  print( 'max. frequency : ' + str( min(cs0,cs1,cp0,cp1)/(ppw*max(dx,dy,dz)) ) )
  print( '---------------------------------\n' )

  # Create discretization matrices using nutils
  K, C, M, rhs = elast_mat(rho, cp, cs, lam, mu, ndims, nx, ny, nz, vec_basis, domain, geom, block)

  t0 = time.time()
  if solver_flag==0:           
      print('Use megmres')
      sol, it = me_driver(K.toscipy().tocsc(), C.toscipy().tocsc(), M.toscipy().tocsc(), rhs, freq, tau, damping, tol, maxit, iLU, rot, fill_factor, plot_resnrm)
            
  elif solver_flag==1:
      print('Use poly_msgmres')
      sol, it = poly_driver(K.toscipy().tocsc(), C.toscipy().tocsc(), M.toscipy().tocsc(), rhs, freq, tau, damping, tol, maxit, dg_pp, iLU, fill_factor, plot_resnrm)
      
  elif solver_flag==2:
      maxit_o = maxit
      tol_o   = tol
      print('Use fom-fgmres')
      sol, it = nested_driver(K.toscipy().tocsc(), C.toscipy().tocsc(), M.toscipy().tocsc(), rhs, freq, tau, damping, 
                              tol_i, tol_o, maxit_i, maxit_o, iLU, fill_factor, plot_resnrm)

  else:
      print('Use pythons built-in solver...')
      sol = np.zeros((Nom, len(vec_basis)), dtype=complex)
      it = -1
      for k in range(0,Nom):
          matrix = K + 1j*om[k]*C - om[k]**2*M
          A = matrix.toscipy().tocsc()
          if ndims==2:
              t0_lu = time.time()
              lu = spla.splu(A)
              print('LU decomposition:'+str(time.time()-t0_lu))
              t0_solve = time.time()
              sol[k,:] = lu.solve(rhs)
              print('solve:'+str(time.time()-t0_solve))
          else:             
              print('Use ILU+GMRES')
              class gmres_counter(object):
                  def __init__(self, disp=True):
                      self._disp = disp
                      self.resvec=[]
                      self.niter = 0
                  def __call__(self, rk=None):
                      self.niter += 1
                      self.resvec.append(rk)
                      if self._disp:
                          print('iter %3i\trk = %s' % (self.niter, str(rk)))            
              t0_lu = time.time()
              invA = spla.spilu( A, fill_factor=1.0)
              invA_x = lambda x: invA.solve(x)
              ilu = spla.LinearOperator(A.shape, invA_x)
              print('ilu setup:'+str(time.time()-t0_lu))
              t0_solve = time.time()
              counter = gmres_counter(disp=True)
              sol[k,:], info = spla.gmres(A, rhs, tol=1e-16, restart=200, maxiter=600, M=ilu, callback=counter)
              it = info
              print('GMRES time:'+str(time.time()-t0_solve))
              print('GMRES info:'+str(counter.niter)+' -- '+str(counter.resvec[-1]))
              
          
  te = time.time()
  print('No iterations: '+str(it)+'     CPU time: '+str(te-t0))        
  test_orig_problem(K.toscipy(), C.toscipy(), M.toscipy(), rhs, om, sol.T) 
      
            
  if plots:
      if(ndims ==2):
          makeplots( domain, geom, Lx, Lz, rho, 'rho', 'Density distribution' )
          #makeplots( domain, geom, Lx, Lz, rho, 'rho', 'Density distribution', imgtype='eps' )        
          #makeplots( domain, geom, Lx, Lz, cp, 'cp', 'c_p [m/s]' )
          #makeplots( domain, geom, Lx, Lz, cs, 'cs', 'c_s [m/s]' )
          
          for k in range(0,Nom):
              disp     = vec_basis.dot( sol[k,:] ) 
              disp_x   = disp[0].real
              disp_z   = disp[-1].real
              makeplots( domain, geom, Lx, Lz, disp_x, 'disp_x'+str(k), 'u_x at {} Hz'.format(freq[k]), lineOn=False )
              makeplots( domain, geom, Lx, Lz, disp_z, 'disp_z'+str(k), 'u_z at {} Hz'.format(freq[k]), lineOn=False )
              #makeplots( domain, geom, Lx, Lz, disp_x, 'disp_x'+str(k), 'u_x at {} Hz'.format(freq[k]), lineOn=False, imgtype='eps' )
              #makeplots( domain, geom, Lx, Lz, disp_z, 'disp_z'+str(k), 'u_z at {} Hz'.format(freq[k]), lineOn=False, imgtype='eps' )
          makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, 'wedge2d')

      elif(ndims==3):
          makevtk(domain, geom, rho, lam, mu, cp, cs, sol, freq, vec_basis, 'wedge3d')

util.run( main )
