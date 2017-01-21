#! /bin/bash
. ~/.bashrc


# # iLU(fill=1.0)+GMRES(tol=1e-16) faster than LU at 20 iters each: 22sec vs. 41sec
# ##########################################################################################################################################
# python3 elast_squares.py --ndims=3 --dx=30 --dy=30 --dz=30 --freq=[1.0,2.0] --Nom=5 --degree=1 --damping=0.0 --maxit=300 --maxit_i=10 \
#                          --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=True --fill_factor=1.0 --block=True \
#                          --plots=True --plot_resnrm=True --solver_flag=1 --nprocs=8
#  
# python3 elast_squares.py --ndims=3 --dx=30 --dy=30 --dz=30 --freq=[1.0,2.0] --Nom=5 --degree=1 --damping=0.0 --maxit=300 --maxit_i=10 \
#                          --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=1.0 --block=True \
#                          --plots=True --plot_resnrm=True --solver_flag=1 --nprocs=8
# ##########################################################################################################################################



# python3 -u elast_squares.py --ndims=3 --dx=30 --dy=30 --dz=30 --freq=[1.0,9.0] --Nom=5 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=True --fill_factor=10.0 --block=True \
#                             --plots=True --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp1.txt
                            
                            
python3 -u elast_squares.py --ndims=2 --dx=10 --dy=10 --dz=10 --freq=[1.0,4.0] --Nom=5 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=0.7 -tau_im=-0.3 iLU=False --fill_factor=10.0 --block=True \
                            --plots=True --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1.txt