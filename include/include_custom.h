using namespace std;

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string.h> 
#include <omp.h>

#include "../other_source/mt19937ar.c"
#include <fftw3.h>
#include "../other_source/spline.c"


enum code1 {give_f,give_f_t,give_f_x,give_f_y,give_f_z,give_f_lap};
#include "scalar_f_global_constants.h"
#include "scalar_f_back.h"

#include "../source/my_classes.cpp"

#include "utilities.h"
#include "scalar_f_initial.h"
#include "scalar_f_evolve.h"
