Python code for the numerical experiments presented in [BvG17].

Usage
------
The code is stored in the subfolder `/python_code`. All experiments can be reproduced with a single run from the command line. For example:

`python3 elast_squares.py --dx=30 --dz=30 --dy=30 --damping=0.0 --dg_pp=0 --freq=[1.0,2.0] --Nom=5 --ndims=3 --block=True --tol=1e-8 --solver_flag=1`

Dependencies
-------------
* [nutils](http://www.nutils.org/):  `pip install git+https://github.com/joostvanzwieten/nutils@955bc67d219496e26b037f47855709a222850d7c`
* [PyAMG](http://pyamg.org/): `pip install --upgrade pyamg`
* NumPy [v 1.8.2], SciPy [v 0.14.0], matplotlib [v 1.4.2]

Declaration
-----------
The [author](http://www.manuelbaumann.de) is a PhD student in Numerical Analysis at TU Delft. My research is focused on linear solvers and preconditioning. I highly encourage experts in geophysics to comment on the numerical results and to [get in touch](mailto:m.m.baumann@tudelft.nl).

References
----------
* [Manuel Baumann and Martin B. Van Gijzen. Efficient iterative methods for multi-frequency wave propagation problems: A comparison study. To appear: Proceedings of INTERNATIONAL CONFERENCE ON COMPUTATIONAL SCIENCE 2017](/literature/iccs17_report.pdf)

[BvG17]: /literature/iccs17_report.pdf
