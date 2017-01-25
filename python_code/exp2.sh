#! /bin/bash
. ~/.bashrc


python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[6.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp2_flag2_68.txt
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[4.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp2_flag2_48.txt                         
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp2_flag2_18.txt
                            
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[6.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp2_flag1_68.txt
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[4.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp2_flag1_48.txt                            
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp2_flag1_18.txt   
                                                        
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[6.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_flag0_68.txt
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[4.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_flag0_48.txt                           
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                           --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                           --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_flag0_18.txt                               
                