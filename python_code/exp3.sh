#! /bin/bash
. ~/.bashrc
           

# Exp. 3
python3 -u elast_squares.py --ndims=3 --dx=15 --dy=15 --dz=15 --freq=[1.0,3.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp3_iLU0.txt 
                            
python3 -u elast_squares.py --ndims=3 --dx=15 --dy=15 --dz=15 --freq=[1.0,3.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=True --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp3_iLU10.txt 
                            
python3 -u elast_squares.py --ndims=3 --dx=15 --dy=15 --dz=15 --freq=[1.0,3.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=True --rot=True --fill_factor=20.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp3_fiLU20.txt  
                            
python3 -u elast_squares.py --ndims=3 --dx=15 --dy=15 --dz=15 --freq=[1.0,3.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=True --rot=True --fill_factor=30.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp3_iLU30.txt               
             
