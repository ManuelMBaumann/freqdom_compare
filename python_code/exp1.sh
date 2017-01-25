#! /bin/bash
. ~/.bashrc


# # Exp 1.1
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=0 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_dp0.txt
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=10 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_dp10.txt
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=5 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_dp5.txt
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=4 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_dp4.txt
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_dp3.txt
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=2 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_dp2.txt
#                             
# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_dp10.txt   




# python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
#                             --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
#                             --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_dp10.txt   
# 
# 537 - 2.40445508355e-08
# conv_megmres000.png
# No iterations: 537     CPU time: 648.585


python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f816_Nom5.txt   

python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_1_f816_Nom5.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp1_2_f816_Nom5.txt                               
                            
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[10,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f1016_Nom5.txt   

python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[10,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_1_f81016_Nom5.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[10,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp1_2_f1016_Nom5.txt                               
                            
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[12,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f1216_Nom5.txt   

python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[12,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_1_f81216_Nom5.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[12,16.0] --Nom=5 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp1_2_f1216_Nom5.txt                               
                            
                            
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f816_Nom15.txt   

python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_1_f816_Nom15.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[8,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp1_2_f816_Nom15.txt                               
                            
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[10,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f1016_Nom15.txt   

python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[10,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_1_f81016_Nom15.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[10,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp1_2_f1016_Nom15.txt                               
                            
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[12,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=0 --nprocs=8  | tee experm/exp1_0_f1216_Nom15.txt   

python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[12,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=1 --nprocs=8  | tee experm/exp1_1_f81216_Nom15.txt  
                            
python3 -u elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[12,16.0] --Nom=15 --degree=1 --damping=0.05 --maxit=300 --maxit_i=20 \
                            --tol=1e-8 --tol_i=1e-1 --dg_pp=1 --tau_re=-100 -tau_im=-0.7 iLU=False --fill_factor=10.0 --block=False \
                            --plots=False --plot_resnrm=True --solver_flag=2 --nprocs=8  | tee experm/exp1_2_f1216_Nom15.txt                              
                            
                            
                            
                            

           