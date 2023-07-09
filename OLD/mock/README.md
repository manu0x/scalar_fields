# Mock folder for some test cases
For testing any developing code on 1-d simplified equations

# Relevant files
classic_imexRK_1*.* : any file with such name tests ImEx codes on 1-d simplified SP-like equation

[classic_imexRK_1.cpp](https://github.com/manu0x/scalar_fields/blob/main/mock/classic_imexRK_1.cpp): implements Fourier based solution of linear system.

[classic_imexRK_1_trapezoid.cpp](https://github.com/manu0x/scalar_fields/blob/main/mock/classic_imexRK_1_trapezoid.cpp): implements trapezoid/implicit Eueler using FD

[classic_imexRK_1_predator_lin_eqs.cpp](https://github.com/manu0x/scalar_fields/blob/main/mock/classic_imexRK_1_lin_eqs.cpp) is FD based implementation which currently uses one time matrix inversion as discussed.

# Compiling
Use 
```
g++ filename.cpp -lm  -lblas -lgsl -lfftw3
```

# Structure of [classic_imexRK_1_predator_lin_eqs.cpp](https://github.com/manu0x/scalar_fields/blob/main/mock/classic_imexRK_1_lin_eqs.cpp)
main function calls run()

run() takes in dt,dx etc., calls initialise() to initialize the system, then calls cr_invert_mat() which creates necessary matrices and does one-time matrix inversion and subsequetly performs ImEx integration.



