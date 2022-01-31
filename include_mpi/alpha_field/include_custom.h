using namespace std;

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string.h> 
#include <string>
#include <cctype>
#include <algorithm> 

#include <omp.h>
#include <mpi.h>
#include <time.h>

#include "../../other_source/mt19937ar.c"
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <hdf5.h>
#include "../../other_source/spline.c"


enum code1 {give_f,give_f_t,give_f_x,give_f_y,give_f_z,give_f_lap};
#include "scalar_f_global_constants.h"

#include "../../source_mpi/alpha_field/potn_approx.cpp"
#include "../../source_mpi/alpha_field/field_approx.cpp"

#include "../../source_mpi/alpha_field/classes_alpha_field.cpp"
#include "../../source_mpi/alpha_field/param_parser.cpp"


#include "scalar_f_back.h"


bool hdf5_format;


//#include "utilities.h"
#include "scalar_f_initial.h"
#include "scalar_f_evolve.h"
