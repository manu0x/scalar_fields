import numpy as np

def pwr_spectrum(f,n,L):
    dk = 2.0*np.pi/L
    xi  = np.fft.fftfreq(n)*n*2*np.pi/L
    xix,xiy,xiz = np.meshgrid(xi,xi,xi)
    kcount = np.zeros((n))
    f_ft = np.square(np.abs(np.fft.fftn(f,f.shape)))
    
    print(f_ft.dtype,f_ft.shape)
    
    k2 = xix*xix + xiy*xiy + xiz*xiz
    pwr_spec = np.zeros((n))
    k_vec = np.zeros((n))
    
    
    
    max_ind = -1
    for i in range(n):
        for j in range(n):
            for l in range(n):
                ind = int(np.sqrt(k2[i,j,l])/dk)
                #print(i,j,l,ind,k2[i,j,l],dk)
                kcount[ind] = kcount[ind]+1
                pwr_spec[ind] = pwr_spec[ind]+f_ft[i,j,l]
                k_vec[ind] = np.sqrt(k2[i,j,l])
                if ind>max_ind:
                    max_ind = ind
    print("max_ind",max_ind)
    
    return pwr_spec[1:max_ind]/kcount[1:max_ind],k_vec[1:max_ind],kcount[1:max_ind]