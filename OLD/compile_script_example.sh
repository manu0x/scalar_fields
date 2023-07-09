#!/bin/bash
# Need to speicify hdf5 and fftw3 parallel libraries.
mpic++ /address/to_folder_containig/source_mpi/fdm3/scalar_f.cpp -I/usr/include/hdf5/openmpi -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.a -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lfftw3_mpi -lfftw3
