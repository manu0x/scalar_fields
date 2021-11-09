#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include "mt19937ar.c"


#include <fftw3.h>

using namespace std;

//////////////////////////////Constants/////////////////////////
#define twopie  2.0*M_PI
#define mass 1.0
#define G 1.0

///////////////////////////////////////////Todo List///////////////////////////////////////////
//////////////   1) Add error check for successful memory allocation in class scalar_field_3d    /////////////////////////////////////////
//////////////   2) Check cosmo ini conditions as per fdm					 ////////////////////////////////////////

double hbar_by_m;

enum code1 {give_f,give_f_t,give_f_x,give_f_y,give_f_z,give_f_lap};
#include "my_classes.cpp"


void ini_rand_field(int * ,double *,double *,double * ,double ,double ,double );
void initialise(int * ,fdm_psi,double ,double ,double ,double);
double ini_power_spec(double );
double dlogD_dloga(double );








int main()
{




   int a[3][3][3];
    int n[3]={2,2,2};

	

 int t[3] {0,0,0};
   scalar_field_3d b(n);
	b.get_field(t,give_f_t);



}


double ini_power_spec(double ksqr)
{
	return (1e-5);
}


double dlogD_dloga(double a)
{
	return (1.0);

}

void initialise(int * ind,fdm_psi &psi,double a0,double ai,double Hi,double omega_dm_ini)
{
      

    
      int ci,i,j,k;
      int xcntr[3]={-1,-1,-1};
      
      double ktmp,maxkmagsqr = 0.0,minkmagsqr = 1e10;
     
      double a = ai;
      double a_t = a*Hi;
      double a_ti = ai*Hi;
     
 
      double L[3],dx[3];
      double dk; 
      int kbins;     
      double lenfac = 1.0;
      int tN = ind[0]*ind[1]*ind[2];
      int kbin_count[tN],kbin_grid[tN];
      double k_grid[tN][3],kmag_grid[tN];
      int n[3]{ind[0],ind[1],ind[2]};
      int loc_ind[3],err_hold;

      double ini_dc[tN],ini_theta[tN];
      double psi_amp,psi_r_val,psi_i_val;
       
	     
	 
	

	double kf = twopie*lenfac/(64.0);

	dx[0] = 1.0; dx[1] =1.0; dx[2] = 1.0;
        L[0] = dx[0]*((double) ind[0]);  L[1] = dx[1]*((double) (ind[1]));  L[2] = dx[2]*((double) (ind[2]));
	dk = 0.01/dx[0]; kbins = 0;
	
	//ini_rand_field();
	//  read_ini_rand_field();
        
	for(ci = 0;ci <tN; ++ci)
	{
		kbin_count[ci]=0;
	}
	
	for(ci = 0;ci <tN; ++ci)
	{
		
		
		if((ci%(n[2]*n[1]))==0)
		 ++xcntr[0];
		if((ci%(n[2]))==0)
		 ++xcntr[1];
		 ++xcntr[2];
		ktmp=0.0;
		
		for(j=0;j<3;++j)
		{	
			

			if((xcntr[j]%n[j])<=(n[j]/2))
				{
					
				 k_grid[ci][j] = ((double)(xcntr[j]%n[j]))/L[j];

				  ktmp+= k_grid[ci][j]*k_grid[ci][j];

					
				}
			else
				{ 
				 k_grid[ci][j] = ((double)((xcntr[j]%n[j])-n[j]))/L[j];

				 ktmp+= k_grid[ci][j]*k_grid[ci][j];

				
				  
				}
		
			 
			
			
		
		}
		
		
		


		
	
			
		if(ktmp>maxkmagsqr)
		maxkmagsqr = (ktmp);
		if((ktmp>0.0)&&(minkmagsqr>ktmp))
		minkmagsqr = ktmp;
		

		kbin_grid[ci] = (int)(sqrt(ktmp)/(dk));
		kmag_grid[ci] = sqrt(ktmp);
		 //printf("yo  %d  %lf\n",kmag_grid[ci],sqrt(ktmp));
		++kbin_count[kbin_grid[ci]];

		if(kbin_grid[ci]>kbins)
		kbins=kbin_grid[ci];
		
			

		
		
		
      	}

	
	ini_rand_field(n,kmag_grid, ini_dc,ini_theta,ai,a0,a_ti);

	for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			loc_ind[0] = i;  loc_ind[1] = j;  loc_ind[2] = k;
			psi_amp = sqrt(omega_dm_ini*pow(a0/ai,3.0)*(1.0+ini_dc[ci]));
			psi_r_val = psi_amp*cos(ini_theta[ci]);
			psi_i_val = psi_amp*sin(ini_theta[ci]);

			err_hold =  psi.update(loc_ind, psi_r_val, psi_i_val);
			


		    }
		  }

		}

	
	

	printf("Initialization Complete.\n");
	printf("\nK details:\n	dk is %lf  per MPc",dk/lenfac);
	/*printf("\n Nyquist Wavenumber is %lf",M_PI/dx[0]);
	printf("\n	Min k_mag is %lf per MPc:: corr lmbda is %.16lf MPc",1.0/(dx[0]*lenfac*((double) n)),dx[0]*lenfac*((double) n));
	printf("\n	Max k_mag is %lf per MPc:: corr lmbda is %.16lf Mpc",sqrt(maxkmagsqr)/lenfac,lenfac/sqrt(maxkmagsqr));
	printf("\n	kbins is %d\n",kbins);

	printf("\nLengthscales:");
	printf("\n	Grid Length is %.6lf MPc",dx[0]*lenfac*((double) n));
	printf("\n	dx is %.16lf MPc\n",dx[0]*lenfac);

	*/
	
	
	  

}


void ini_rand_field(int * ind,double *kmag_grid,double * ini_dc,double * ini_theta_return,double a,double a0,double a_t)
{	init_genrand(time(0));
	int i,cnt,tN; 
	double ksqr,muk,sigk;
	double a1,a2,b1,b2,a_rand,b_rand;
	double f_ini;
	


	FILE *fpinirand = fopen("initial_rand_field.txt","w");

	
	
	tN = ind[0]*ind[1]*ind[2];

	fftw_plan ini_del_plan;
	fftw_plan ini_theta_plan;

	fftw_complex *F_ini_del;
	fftw_complex *ini_del;
	fftw_complex *F_ini_theta;
	fftw_complex *ini_theta;

	F_ini_del = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	ini_del = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	F_ini_theta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);
	ini_theta = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) *tN);

	init_genrand(time(0));

	f_ini  = dlogD_dloga(a);
	
	for(cnt=0;cnt<tN;++cnt)
	{	
		 	    ksqr = kmag_grid[cnt];
			    sigk  = sqrt(ini_power_spec(sqrt(ksqr)));
			    muk = sigk/sqrt(2.0);
		 	    a1 = genrand_res53();
 			    a2 = genrand_res53(); 
			   // b1 = genrand_res53();
 			  //  b2 = genrand_res53();
			    a_rand = (muk*(sqrt(-2.0*log(a1))*cos(2.0*M_PI*a2)));
			    b_rand = (muk*(sqrt(-2.0*log(a1))*cos(2.0*M_PI*a2)));
				
			    F_ini_del[cnt][0] = a_rand;	F_ini_del[cnt][1] = b_rand;
			    if(ksqr>0.0)
			    { F_ini_theta[cnt][0] =  mass*(a_t/a)*f_ini*(a/a0)*(a/a0)*F_ini_del[cnt][0];
			      F_ini_theta[cnt][1] =  mass*(a_t/a)*f_ini*(a/a0)*(a/a0)*F_ini_del[cnt][1];
			    }	



	}




	ini_del_plan = fftw_plan_dft_3d(ind[0],ind[1],ind[2], F_ini_del, ini_del, FFTW_BACKWARD, FFTW_ESTIMATE);
	ini_theta_plan = fftw_plan_dft_3d(ind[0],ind[1],ind[2], F_ini_theta, ini_theta, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	

	fftw_execute(ini_del_plan);
	fftw_execute(ini_theta_plan);
	

	
	for(cnt=0;cnt<tN;++cnt)
	{
		
		ini_del[cnt][0] = ini_del[cnt][0]/sqrt(tN); ini_del[cnt][1] = ini_del[cnt][1]/sqrt(tN); 
		ini_theta[cnt][0] = ini_theta[cnt][0]/sqrt(tN); ini_theta[cnt][1] = ini_theta[cnt][1]/sqrt(tN); 
		ini_dc[cnt] = ini_del[cnt][0];
		ini_theta_return[cnt] = ini_theta[cnt][0];
		

		fprintf(fpinirand,"%d\t%.16lf\t%.16lf\t%.16lf\n",
					cnt,kmag_grid[cnt],ini_del[cnt][0],ini_theta[cnt][0]);


		

	}
    

	 fftw_free(F_ini_del);
	 fftw_free(ini_del);
	 fftw_destroy_plan(ini_del_plan);
	
	 fftw_free(F_ini_theta);
	 fftw_free(ini_theta);
	 fftw_destroy_plan(ini_theta_plan);
	

	
}



