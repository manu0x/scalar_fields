



import sys

import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
import importlib

pie = np.pi




n_high = int(sys.argv[1])
n_low = int(sys.argv[2])

f_in = str(sys.argv[3])
#f_out = str(sys.argv[4])

use_zeldo = int(sys.argv[4])

finlist = f_in.split('_')
n_name = finlist[1]

foutlist = []
i=0
for i in range(len(finlist)):
    if i==1:
        foutlist.append(str(n_low))
    else:
        foutlist.append(finlist[i])

f_out = "_".join(foutlist)

print(f_in,f_out)
fp_in = h5py.File(f_in+".hdf5", "r")

psi_in_zeldo = fp_in["psi_zeldo"]

psi_in = fp_in["psi"]

psi_in_zeldo_t  = np.reshape(psi_in_zeldo,(n_high,n_high,n_high,2))
psi_in_t  = np.reshape(psi_in,(n_high,n_high,n_high,2))

skp = int(n_high/n_low)

psi_out_zeldo  = psi_in_zeldo_t[::skp,::skp,::skp,:]
psi_out  = psi_in_t[::skp,::skp,::skp,:]

psi_out_zeldo = np.reshape(psi_out_zeldo,(n_low*n_low*n_low,2))
psi_out = np.reshape(psi_out,(n_low*n_low*n_low,2))

with h5py.File(f_out+".hdf5", "w") as ff:
    ff.create_dataset("psi_zeldo",psi_out_zeldo.shape,data=psi_out_zeldo)
    ff.create_dataset("psi",psi_out.shape,data=psi_out)

print(psi_in.shape,psi_out.shape,psi_in_zeldo.shape,psi_out_zeldo.shape,skp,n_name)






