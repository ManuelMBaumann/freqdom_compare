#! /bin/bash
. ~/.bashrc

                           
                            
# Exp. 2.1 
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[7.0,8.0] --Nom=5 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p1_f78_Nom5.txt   

python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[4.0,8.0] --Nom=5 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p1_f48_Nom5.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,8.0] --Nom=5 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p1_f18_Nom5.txt  


python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[7.0,8.0] --Nom=15 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p1_f78_Nom15.txt  
            
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[4.0,8.0] --Nom=15 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p1_f48_Nom15.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,8.0] --Nom=15 --degree=1 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p1_f18_Nom15.txt                            


# Exp. 2.2
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[7.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p2_f78_Nom15.txt  
            
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[4.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p2_f48_Nom15.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=5 --dy=5 --dz=5 --freq=[1.0,8.0] --Nom=15 --degree=2 --damping=0.0 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 iLU=False --rot=True --fill_factor=10.0 --block=True \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp2_p2_f18_Nom15.txt             
                