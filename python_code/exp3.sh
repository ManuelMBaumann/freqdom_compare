#! /bin/bash
. ~/.bashrc

python3 -u elast_squares.py --ndims=3 --dx=15 --dy=15 --dz=15 --freq=[1.0,3.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=11.1 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1.txt     
python3 -u elast_squares.py --ndims=3 --dx=15 --dy=15 --dz=15 --freq=[1.0,3.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=True --fill_factor=10.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1.txt      
python3 -u elast_squares.py --ndims=3 --dx=15 --dy=15 --dz=15 --freq=[1.0,3.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=True --fill_factor=20.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1.txt     
python3 -u elast_squares.py --ndims=3 --dx=15 --dy=15 --dz=15 --freq=[1.0,3.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=True --fill_factor=30.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1.txt     
                            
# python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,4.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
#                            --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=11.1 --block=True \
#                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1.txt                                
# python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,4.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
#                            --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=True --fill_factor=15.0 --block=True \
#                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1.txt      
# python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,4.0] --Nom=10 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
#                            --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=True --fill_factor=10.0 --block=True \
#                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1.txt    
