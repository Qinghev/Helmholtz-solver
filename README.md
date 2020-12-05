# Helmholtz-solver
Poisson and Helmholtz equation solver.

The sovler is based on fishpack90, you can find it through https://www2.cisl.ucar.edu/resources/legacy/fishpack90/documentation#hwscrt.txt

Because it is single precision, so I changed all the original variables to real(kind=8) 'double precision', and 'mytest_hwscrt.f90' is the file tests the equation
$$u_{xx}+u_{yy}+\lambda u = g,$$
the exact solution is $u(x,y) = \exp(x)\sin(\pi y)$.
