#!/usr/bin/env python
# coding: utf-8



import sys
import tensorflow as tf
from tensorflow import keras
import matplotlib.pyplot as plt
import numpy as np
import h5py
import os
import importlib

pie = np.pi




n = int(sys.argv[1])
L = float(sys.argv[2])
snap_name = str(sys.argv[3])
param_file = str(sys.argv[4])

par_path = os.getcwd()
sys.path.append(par_path)


param = importlib.import_module(param_file, package=None)
print(param.c_unit)
print("sys_args ",n,L,snap_name)



hbar_box = param.hbar_unit
c_box = param.c_unit
pc_box = param.pc_unit
h = param.h
alpha = param.m_alpha
hbar_by_m_val = (hbar_box*h*c_box*c_box/(alpha*pc_box))*0.001
H0_val = 100.0
#da_val = 0.00001
a0_val = 1.0
ai_val = param.ai

print("\n#############PARAM loaded############")
print('hbar',hbar_box)
print('c_box',c_box)
print('pc_box',pc_box)
print('h',h)
print('alpha',alpha)
print('hbar_by_m',hbar_by_m_val)
print("ai",ai_val)
print("#####################################\n")




a0 =  tf.constant(a0_val,dtype=tf.float64)
ai =  tf.constant(ai_val,dtype=tf.float64)

omega_dm0 =  tf.constant(0.3,dtype=tf.float64)

H0 =  tf.constant(100.000,dtype=tf.float64)




def pcube(p):
    pivot = np.floor(p/dx).astype(int)
    a001 = np.array([0,0,1])
    a010 = np.array([0,1,0])
    a100 = np.array([1,0,0])
    
    a011 = np.array([0,1,1])
    a110 = np.array([1,1,0])
    a101 = np.array([1,0,1])
    
    a111 = np.array([1,1,1])
    
    p001 = pivot+a001
    p010 = pivot+a010
    p100 = pivot+a100
    p011 = pivot+a011
    p110 = pivot+a110
    p101 = pivot+a101
    
    p111 = pivot+a111
    
    part_cube =np.stack([pivot,p100,p010,p001,p110,p011,p101,p111],axis=0)

    return part_cube


def pic_dc_construct(p,v,n,dx):
    
    
    n_part = np.max(p.shape)
    print("n_part",n_part)
    
    dc=np.zeros((n,n,n))
    vel=np.zeros((n,n,n,3))
    
    for i in range(n_part):
        if i%(n*n)==0:
            print("i",i/(n*n))

        p_cube_ind = pcube(p[i,:])
        p_cube = p_cube_ind*dx
        
        d = np.prod((1.0-np.absolute(p[i,:]-p_cube)/dx),axis=1)
        
        for j in range(8):
            m = (n+(p_cube_ind[j,0]))%n
            q = (n+(p_cube_ind[j,1]))%n
            k = (n+(p_cube_ind[j,2]))%n
            dc[m,q,k] = dc[m,q,k] + d[j]
            vel[m,q,k] = vel[m,q,k]+d[j]*v[i]
            
    dc = dc-1.0
    print("dc shape",dc.shape,"v shape",vel.shape)
            
    return dc,vel



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







dx = L/(float(n))

pdata = h5py.File(snap_name)
par_p = np.array(pdata['PartType1']['Coordinates'])
par_v = np.array(pdata['PartType1']['Velocities'])

print("Loaded\n",par_v.shape,par_p.shape)

dc,v = pic_dc_construct(par_p,par_v,n,dx)




outname = "dc_"+str(n)




def cal_psi_from_dc_v_2(dc,v,ai,n,L,hbar_by_m,H0,f_ini,omega_dm0=0.3,a0=1.0):
    rho_b = 3.0*H0*H0*omega_dm0*np.power(a0/ai,3.0) 
    print("dc type",type(dc),"vtype",type(v),type(H0),type(omega_dm0),type(ai),type(rho_b))
    psi_amp = np.sqrt(rho_b*(1.0+dc))

    H_ini = H0*np.sqrt( omega_dm0*np.power(a0/ai,3.0) + (1.0-omega_dm0)   )
    v = ai*v
   
    v_ft_x = np.fft.fftn(v[:,:,:,0]+0*1j,v[:,:,:,0].shape)
    print(v[:,:,:,0].shape)
    v_ft_y = np.fft.fftn(v[:,:,:,1]+0*1j,v[:,:,:,1].shape)
    v_ft_z = np.fft.fftn(v[:,:,:,2]+0*1j,v[:,:,:,2].shape)

    xi = np.fft.fftfreq(n)*n*2*np.pi/L
    xix,xiy,xiz = np.meshgrid(xi,xi,xi)

    theta_zel_rhs = -H_ini*f_ini*ai*ai*dc/hbar_by_m
    theta_zel_ft = -np.fft.fftn(theta_zel_rhs,theta_zel_rhs.shape)/(xix*xix+xiy*xiy+xiz*xiz)
    theta_zel_ft[0,0,0] = 0.0+0.0*1j

    theta_zel = np.fft.ifftn(theta_zel_ft,theta_zel_ft.shape)
    
    
    
    grad_v = 1j*(xix*v_ft_x +xiy*v_ft_y+xiz*v_ft_z)
    #print("GRAD",grad_v.dtype,hbar_by_m_c)
    alpha_rhs = -grad_v/hbar_by_m
    

    #print("ALPHA",alpha_rhs.shape,h_filt_c.shape,alpha_rhs[0,0,0])
    
    alpha = np.fft.ifftn(alpha_rhs,alpha_rhs.shape)
    
    alpha = ai*alpha
    
    
    psi =  psi_amp*np.cos(np.real(alpha))+1j*psi_amp*np.sin(np.real(alpha))
    psi_zeldo =  psi_amp*np.cos(np.real(theta_zel))+1j*psi_amp*np.sin(np.real(theta_zel))

    psi = np.power(ai_val,1.5)*psi/H0
    psi_zeldo = np.power(ai_val,1.5)*psi_zeldo/H0

    psi_zeldo_r = psi_zeldo.real
    psi_zeldo_i = psi_zeldo.imag

    psi_zeldo_r = np.reshape(psi_zeldo_r,(n*n*n))
    psi_zeldo_i = np.reshape(psi_zeldo_i,(n*n*n))
    psi_zeldo = np.stack([psi_zeldo_r,psi_zeldo_i],axis=-1)


    psi_r = psi.real
    psi_i = psi.imag

    psi_r = np.reshape(psi_r,(n*n*n))
    psi_i = np.reshape(psi_i,(n*n*n))
    psi = np.stack([psi_r,psi_i],axis=-1)


    alpha_r = alpha.real
    alpha = np.reshape(alpha_r,(n*n*n,1))

    theta_zel_r = theta_zel.real
    theta_zel = np.reshape(theta_zel_r,(n*n*n,1))
    

    print("PSI shapes",psi.shape,psi_zeldo.shape,alpha.shape,theta_zel.shape)
    
    
    
    return psi,alpha,psi_zeldo,theta_zel











psi_ini2,alpha_ini2,psi_zeldo,theta_zel = cal_psi_from_dc_v_2(dc,v,ai_val,n,L,hbar_by_m_val,H0_val,1.0,omega_dm0.numpy() )









dc = np.array(dc).flatten()

dc = np.reshape(dc,(dc.shape[0],1))





with h5py.File(outname+"_dc_theta_psi_zeldo.hdf5", "w") as ff:
    dsetdc = ff.create_dataset("psi",psi_ini2.shape,data=psi_ini2)
    dsetdc = ff.create_dataset("psi_zeldo",psi_zeldo.shape,data=psi_zeldo)
    dsetdc1 = ff.create_dataset("dc", dc.shape,data=dc)
    dset2 = ff.create_dataset("theta", alpha_ini2.shape,data=alpha_ini2)
    dset3 = ff.create_dataset("theta_zel", theta_zel.shape,data=theta_zel)



    











