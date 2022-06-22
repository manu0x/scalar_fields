


import tensorflow as tf
from tensorflow import keras
import matplotlib.pyplot as plt
import numpy as np
import h5py
pie = np.pi

from fdm_evolver import make_k_grid 
from fdm_evolver import make_x_grid 
from fdm_evolver import fdm_evolver_a 
from fdm_evolver import cal_dc_from_psi 
from fdm_evolver import a_t_calc


################# SIMULATION PARAMETERS  #################################################

n=128
L=1.0



hbar_box = 6.582
c_box = 2.99
pc_box = 3.0857
h = 0.7
alpha = 0.175
hbar_by_m_val = (hbar_box*h*c_box*c_box/(alpha*pc_box))*0.001
H0_val = 100.0
da_val = 0.00001
a0_val = 1.0
ai_val = 0.0078125
a_print = (int)((a0_val-ai_val)/(10.0*da_val))
print(a_print,hbar_by_m_val)
hbar_by_m_val*pie*128


Pie = tf.constant(np.pi,dtype=tf.float64)


a = tf.constant(0.0,dtype=tf.float64)
b = tf.constant(1.0,dtype=tf.float64)
comp_i = tf.complex(a,b)


a0 =  tf.constant(a0_val,dtype=tf.float64)
ai =  tf.constant(ai_val,dtype=tf.float64)
m=   tf.constant(1.0,dtype=tf.float64)
omega_dm0 =  tf.constant(0.3,dtype=tf.float64)

H0 =  tf.constant(100.000,dtype=tf.float64)



hbar_by_m = tf.constant(hbar_by_m_val,dtype=tf.float64)
hbar_by_m_c = tf.complex(hbar_by_m,tf.zeros_like(hbar_by_m))


##############################################################################################################################################



########################### CALCULATING da #####################################


att = a_t_calc(ai,H0,omega_dm0)
da_lim = (4.0/(3.0*pie))*ai*ai*att*(L/n)*(L/n)/hbar_by_m
da = tf.constant(da_lim,dtype=tf.float64)
print(ai.numpy(),da_lim.numpy(),att.numpy(),hbar_by_m.numpy())

##############################################################################################################################################

f = h5py.File("dc.hdf5", "r")
dc = tf.constant(f["/dc"],dtype=tf.float64)
#dc = tf.zeros_like(dc)
v = tf.constant(f["/v"],dtype=tf.float64)
v = tf.transpose(v,perm=[3,0,1,2])
print(dc.shape,v.shape)
h5py.File.close(f)





def loaf_ini_file(file="initial1.txt",n=128):
    pl_r = np.genfromtxt(file,delimiter="\t",usecols=4)
    pl_i = np.genfromtxt(file,delimiter="\t",usecols=5)
    potn = np.genfromtxt(file,delimiter="\t",usecols=6)
    #plt.plot(pl_r)
    
    pl_r = np.reshape(pl_r,(n,n,n))
    pl_i = np.reshape(pl_i,(n,n,n))
    potn = np.reshape(potn,(n,n,n))
    
    pc = tf.complex(pl_r,pl_i)
    potn_c = tf.complex(potn,tf.zeros_like(potn))
    
    return pc,potn





#f,p = loaf_ini_file()





dx = L/n
nh = (int)(n/2+1)
print(nh)

k_grid,kk = make_k_grid(L,n)

k_grid_den = 1.0*k_grid
k_grid_den[:,0,0,0] = 1.0

k_grid_half = k_grid[:,:,:,:nh]
k_grid_half_den = k_grid_den[:,:,:,:nh]



k_grid_sqr = np.sum(np.square(k_grid),axis=0)
k_grid_sqr_c = tf.complex(k_grid_sqr,tf.zeros_like(k_grid_sqr))
k_grid_sqr_den = np.sum(np.square(k_grid_den),axis=0)
k_grid_sqr_half_den = np.sum(np.square(k_grid_half_den),axis=0)

k_grid_sqr_den[0,0,0] = 1.0;
k_grid_sqr_half_den[0,0,0] = 1.0;

k_grid_sqr_half_den_c = tf.complex(k_grid_sqr_half_den,tf.zeros_like(k_grid_sqr_half_den))
k_grid_sqr_den_c = tf.complex(k_grid_sqr_den,tf.zeros_like(k_grid_sqr_den))


x_grid = make_x_grid(L,n)

k_grid_sqr_half = tf.reduce_sum(tf.square(k_grid_half),axis=0)
k_grid_sqr_half_c = tf.complex(k_grid_sqr_half,tf.zeros_like(k_grid_sqr_half))


h_filt = np.ones(shape=(n,n,nh))
h_filt_z = np.zeros_like(h_filt)
h_filt[0,0,0]=0.0
h_filt_c = tf.complex(h_filt,h_filt_z)

filt = np.ones(shape=(n,n,n))
filt_z = np.zeros_like(filt)
filt[0,0,0]=0.0
filt_c = tf.complex(filt,filt_z)






def write_fields(a,psi,V):
    name = "data_"+".h5"
    #dc = tf.abs(psi)
    
    #with h5py.File(name, "w") as f:
    dset_psi = f.create_dataset("psi", psi.shape,data=psi, dtype='lf')
     #   dset_V = f.create_dataset("V", V.shape,data = V, dtype='lf')
      #  dset_dc = f.create_dataset("dc", dc.shape, data = dc, dtype='lf')



psi_ini,V_ini = loaf_ini_file()


V_ini = ai*V_ini

my_fdm_a = fdm_evolver_a(hbar_by_m,da,k_grid,k_grid_sqr_c,k_grid_sqr_half_den_c,h_filt_c)



psi_f,V_f = my_fdm_a.evolve(ai,0.008,psi_ini,V_ini,H0,omega_dm0)
            




dc_f  = np.array(cal_dc_from_psi(psi_f,1.00,H0,omega_dm0))

with h5py.File("f_final.hdf5", "w") as ff:
    dset = ff.create_dataset("psi_f", psi_f.shape,data=psi_f)
    dsetdc = ff.create_dataset("dc_f", dc_f.shape,data=dc_f)
    dset2 = ff.create_dataset("V_f", V_f.shape,data=V_f)
    dset3 = ff.create_dataset("V_ini", V_ini.shape,data=V_ini)





fl = h5py.File("f_final.hdf5","r")




dc_f=np.array(fl['dc_f'])


print("###########","   ","RAN","  ","#############")


#ap = 100
#plt.imshow(dc_f[:,ap,:])




