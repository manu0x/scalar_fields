


import tensorflow as tf
from tensorflow import keras
import matplotlib.pyplot as plt
import numpy as np
import h5py
pie = np.pi


a = tf.constant(0.0,dtype=tf.float64)
b = tf.constant(1.0,dtype=tf.float64)
comp_i = tf.complex(a,b)


############################################


def a_t_calc(a,H0,omega_dm0,a0=1.0):
    a_t = a*H0*tf.sqrt(omega_dm0*tf.pow(a0/a,3.0)+(1.0-omega_dm0))
    return a_t







def make_k_grid(L,n):
    dk = 1.0/(L)
    k1 = np.arange(0,n/2+1)
    k2 = np.arange(-n/2+1,0,1)
    k = 2.0*np.pi*dk*np.concatenate((k1,k2))
    k_grid = np.meshgrid(k,k,k,indexing="ij")
    return np.array(k_grid),k

def make_x_grid(L,n):
    dx = 1.0/(L)
    x = np.arange(0,n)
    x = dx*x
    x_grid = np.meshgrid(x,x,x,indexing="ij")
    return np.array(x_grid)














class fdm_evolver_a():
    def __init__(self,hbar_by_m,da,k_grid,k_grid_sqr_c,k_grid_sqr_half_den_c,h_filt_c):
        self.hbar_by_m = hbar_by_m
        self.da = da
        
        
        self.hbar_by_m_c = tf.complex(hbar_by_m,tf.zeros_like(hbar_by_m))
        self.da_c = tf.complex(da,tf.zeros_like(da))
        self.k_grid_sqr_c = k_grid_sqr_c
        self.k_grid_sqr_half_den_c = k_grid_sqr_half_den_c
        self.h_filt_c = h_filt_c
        
        
    #@tf.function    
    def evolve(self,a_ini,a_fin,psi_ini,pot_ini,H0=tf.constant(100.0,dtype=tf.float64),omega_dm0=tf.constant(0.3,dtype=tf.float64),a0=tf.constant(1.0,dtype=tf.float64)):
        a = a_ini
        psi = psi_ini
        pot = pot_ini
        
        #psi_st = tf.expand_dims(psi,axis=0)
        #pot_st = tf.expand_dims(pot,axis=0)
        
        a_cnt = 0
        
        a_t = a_t_calc(a,H0,omega_dm0,a0)
        
        a_t_c = tf.complex(a_t,tf.zeros_like(a_t))
        
        
        epow1 = (0.5*self.da/(a*a_t*self.hbar_by_m))*pot
        efac1 = tf.exp(-comp_i*tf.complex(epow1,tf.zeros_like(epow1)))
        psi = efac1*psi
        
        print("################ STARTING EVOLUTION ######################")
        
        while(a<=(a_fin-self.da)):
            if(a_cnt%10==0):
                tf.print("time is ",a)
                #psi_st_new = tf.expand_dims(psi,axis=0)
                #psi_st = tf.concat([psi_st,psi_st_new],axis=0)
            
            a_t = a_t_calc(a,H0,omega_dm0,a0)
            a_t_c = tf.complex(a_t,tf.zeros_like(a_t))
            a_c = tf.complex(a,tf.zeros_like(a))
            
            

            
            
            efac_diff = tf.exp(-comp_i*0.5*self.da_c*self.k_grid_sqr_c*self.hbar_by_m_c/(a_c*a_c*a_t_c))
            fft_psi = efac_diff*tf.signal.fft3d(psi)
            ifft_psi = tf.signal.ifft3d(fft_psi)
            abs2_psi = tf.square(tf.abs(ifft_psi))
            pois_arg = tf.signal.rfft3d(0.5*abs2_psi-1.5*H0*H0*omega_dm0*a0*a0*a0)
            pot = tf.signal.irfft3d(self.h_filt_c*(-pois_arg/self.k_grid_sqr_half_den_c))
            epow1 = (self.da/(a*a_t*self.hbar_by_m))*pot
            efac1 = tf.exp(-comp_i*tf.complex(epow1,tf.zeros_like(epow1)))
            psi = efac1*ifft_psi
            
            a=a+self.da
            
            a_cnt = a_cnt+1
            
            

            
        efac_diff = tf.exp(-comp_i*0.5*self.da_c*self.k_grid_sqr_c*self.hbar_by_m_c/(a_c*a_c*a_t_c))
        fft_psi = efac_diff*tf.signal.fft3d(psi)
        ifft_psi = tf.signal.ifft3d(fft_psi)
        abs2_psi = tf.square(tf.abs(ifft_psi))
        pois_arg =tf.signal.rfft3d(0.5*abs2_psi-1.5*H0*H0*omega_dm0*a0*a0*a0)
        pot = tf.signal.irfft3d(self.h_filt_c*(-pois_arg/self.k_grid_sqr_half_den_c))
        epow1 = (0.5*self.da/(a*a_t*self.hbar_by_m))*pot
        efac1 = tf.exp(-comp_i*tf.complex(epow1,tf.zeros_like(epow1)))
        psi = efac1*ifft_psi

        print("################ ENIDING EVOLUTION ######################")
        
        return psi,pot




def cal_dc_from_psi(psi,a,H0,omega_dm0):
    abs2_psi = tf.square(tf.abs(psi))
    dc = abs2_psi/(3.0*H0*H0*omega_dm0) - 1.0
    return dc



