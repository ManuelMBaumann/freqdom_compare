[![DOI](https://zenodo.org/badge/79431988.svg)](https://zenodo.org/badge/latestdoi/79431988)

Python code for the numerical experiments presented in [BvG17](http://www.sciencedirect.com/science/article/pii/S1877050917306294).

Usage
------
The code is stored in the subfolder `/python_code`. All experiments can be reproduced with a single run from the command line. For example:

`python3 elast_squares.py --ndims=2 --dx=2.5 --dy=2.5 --dz=2.5 --freq=[1.0,9.0] --Nom=5 --degree=1 --damping=0.005 --maxit=300 --maxit_i=20
                          --tol=1e-8 --tol_i=1e-1 --dg_pp=3 --tau_re=0.7 --tau_im=-0.3 --iLU=True --fill_factor=10.0 --rot=True --block=True
                          --plots=True --plot_resnrm=True --solver_flag=0 --nprocs=8`

Dependencies
-------------
* [nutils](http://www.nutils.org/):  `pip install git+https://github.com/joostvanzwieten/nutils@955bc67d219496e26b037f47855709a222850d7c`
* [PyAMG](http://pyamg.org/): `pip install --upgrade pyamg`
* NumPy [v 1.8.2], SciPy [v 0.14.0], matplotlib [v 1.4.2]

Declaration
-----------
The [author](http://www.manuelbaumann.de) is a PhD student in Numerical Analysis at TU Delft. My research is focused on linear solvers and preconditioning for the elastic wave equation. Feel free to [get in touch](mailto:m.m.baumann@tudelft.nl).

References
----------
* Manuel Baumann and Martin B. van Gijzen (2017). [Efficient iterative methods for multi-frequency wave propagation problems: A comparison study.](http://www.sciencedirect.com/science/article/pii/S1877050917306294) Procedia Computer Science, Volume 108, Pages 645-654.
