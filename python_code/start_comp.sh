#! /bin/bash

python3  elast_squares.py --ndims=3 --dx=30 --dy=30 --dz=30 --freq=[1.0,9.0] --Nom=5 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=True --fill_factor=10.0 --block=True \
                            --plots=True --plot_resnrm=True --solver_flag=2 --nprocs=4
