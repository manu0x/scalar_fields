# Mock folder for some test cases
For testing any developing code on 1-d simplified equations

# Relevant files
classic_imexRK_1*.* : any file with such name tests ImEx codes on 1-d simplified SP-like equation
classic_imexRK_1.cpp: implements Fourier based solution of linear system.
classic_imexRK_1_trapezoid.cpp: implements

# Compiling 
See given example compile script. One needs to provide path to hdf5 parallel libs and fftw3 libs.
# Running
Just run generated objective file(e.g. a.out) using 
```
mpirun -np n path/to/a.out p.txt 

```
where n i no. of processes and "p.txt" is parameter file. See example folder for examples of compile scripts and parameter file.
# Generating initial conditions
Initial conditions are needed in simulation folder as file "dc.hdf5". The code yet does not have a well tested initial conditions generator. We are currently using NGenIC (https://www.h-its.org/2014/11/05/ngenic-code/) through Gadget4's (https://gitlab.mpcdf.mpg.de/vrs/gadget4) initial condition generator. Gadget produces snapshot files. Using utilities in folder [py_ut/initial_cond_utitlities/](https://github.com/manu0x/scalar_fields/tree/main/py_ut/initial_cond_utitlities) one can generate initial conditions file needed by this code. For details of using initial condition utilities, read the docs in that folder.  
