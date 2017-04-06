#! /bin/bash
. ~/.bashrc



# python3 elast_squares.py --ndims=2 --dx=50.0 --dy=50.0 --dz=50.0 --freq=[1.0,9.0] --Nom=5 --degree=1 --damping=0.5 --maxit=300 --maxit_i=20 \
#                          --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=-100 --tau_im=-0.7 --iLU=False --rot=True --fill_factor=10.0 --block=False \
#                          --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8

                         

# Exp. 1.2
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[12.0,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=-100 --tau_im=-0.7 --iLU=False --rot=True --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_0_f1216_Nom5.txt  
                            
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[10.0,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=-100 --tau_im=-0.7 --iLU=False --rot=True --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f1016_Nom5.txt  
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8.0,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=-100 --tau_im=-0.7 --iLU=False --rot=True --fill_factor=10.0 --block=False \
#                             --plots=True --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f816_Nom5.txt            
#                             
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[12.0,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=-100 --tau_im=-0.7 --iLU=False --rot=True --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f1216_Nom15.txt  
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[10.0,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=-100 --tau_im=-0.7 --iLU=False --rot=True --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f1016_Nom15.txt  
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8.0,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=-100 --tau_im=-0.7 --iLU=False --rot=True --fill_factor=10.0 --block=False \
#                             --plots=True --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f816_Nom15.txt   
                            