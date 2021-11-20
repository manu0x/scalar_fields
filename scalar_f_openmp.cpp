#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <string.h> 
#include "mt19937ar.c"
#include <omp.h>

#include <fftw3.h>

using namespace std;

//////////////////////////////Constants/////////////////////////
#define twopie  2.0*M_PI
#define mass 1.0
#define G 1.0

///////////////////////////////////////////Todo List///////////////////////////////////////////
//////////////   1) Add error check for successful memory allocation in class scalar_field_3d    /////////////////////////////////////////
//////////////   2) Check cosmo ini conditions as per fdm					 /////////////////////////////////////////
//////////////   3) Check k grid related things							 /////////////////////////////////////////


double hbar_by_m,h,H0,lenfac;
double c_box,pc_box,hbar_box;

enum code1 {give_f,give_f_t,give_f_x,give_f_y,give_f_z,give_f_lap};
#include "my_classes.cpp"


void ini_rand_field(int * ,double *,double *,double * ,double ,double ,double,ini_power_generator );
void initialise(int * ,fdm_psi &,metric_potential &,double [][3],double ,double ,double ,double,double *,ini_power_generator);
double ini_power_spec(double );
double dlogD_dloga(double );
void set_back_cosmo(double &,double &,double &,double &);
int evolve_kdk(int *,fdm_psi &,metric_potential &,double [][3],double ,double,double,double,double *,double dt=1e-4);







int main()
{
	int c1=1,c2=1,fail=0;
	int ind[3]{128,128,128};	
	int tN = ind[0]*ind[1]*ind[2];
	fdm_psi psi(ind,true);
	metric_potential phi(ind,true);

	const char name[] = "ax_test_matterpower.dat";
	ini_power_generator gen(name);
	gen.check(3.1);

	

	double a0,ai,Hi,omega_dm_ini;
	double k_grid[tN][3],dx[3];
	
	set_back_cosmo(a0,ai,Hi,omega_dm_ini);
	printf("Hi %lf\nOmega_dm_ini %lf\nai %lf\n",Hi,omega_dm_ini,ai);
	initialise(ind,psi,phi,k_grid,a0,ai,Hi,omega_dm_ini,dx,gen);
	fail = evolve_kdk(ind,psi,phi,k_grid,a0,ai,a0,omega_dm_ini,dx,0.8e-4);
	printf("fail is %d\n",fail);
	
	
}


void set_back_cosmo(double &a0,double &ai,double &Hi,double &omega_dm_ini)
{
	double z = 99.0;
	double alpha = 1.0;// Mass in 10^(-22) eV;	

	c_box = 2.99;
	pc_box = 3.0857;
	hbar_box = 6.582119569;//eVs 
	lenfac = 1.0;
	omega_dm_ini = 0.29;
	h = 0.7;
	hbar_by_m = hbar_box*c_box*(1e-8)/(alpha*pc_box);	
		
	
	
	 H0 = lenfac*(h/c_box)*0.001;
	printf("H0 %lf\n",H0);
	
	
	a0 = 1.0;
	ai = a0/(1.0+z);
	Hi =   H0*sqrt(omega_dm_ini*pow(a0/ai,3.0)+ (1.0-omega_dm_ini));


}





int evolve_kdk(int *n,fdm_psi &psi,metric_potential &phi,double k_grid[][3],
					double a_final,double a_ini,double a0,double omega_dm_ini,double *dx,double dt)
{	printf("Yo\n");
	double a,a_t,t,ak,a3a03omega,dti=dt;
	double a_print;
	double psi_vel[2],psi_val[n[0]*n[1]*n[2]][2],psi_upd[2],psi_k[2],potn,poisson_rhs,psi_amp;
	double psi_retrive[2];
	int i,j,k,ci,ind[3];
	int c1,c2,fail=0;
	int step_cnt;	
	int num_prints = 14;
	double z_print_list[num_prints]{99.0,40.0,10.0,5.0,4.0,5.0,3.0,2.0,1.0,0.5,0.25,0.1,0.000001,0.0};
	double z_cur;
	
	int prn=0;	


	FILE *fp_psi = fopen("psi.txt","w");

	for(a=a_ini,a_print=a_ini,step_cnt=0;(a<a_final)&&(!fail);t+=dt,++step_cnt)
	{
	   //dt=dti*sqrt(a/a_ini);
	 if(a>=a_print)
	  {
		
		a_print+=1e-3;
		a3a03omega = pow(a/a0,3.0)/omega_dm_ini;
		z_cur = ((a0/a) -1.0);
		
		printf("a %lf %lf\n",a,a3a03omega);
		
		if(z_cur<=z_print_list[prn])
		{ 
			psi.write_psi(fp_psi,dx,a3a03omega,a,true, true);
			printf("\nWriting at z = %lf\n",z_cur);
			++prn;

		}
		

	  }
	a_t = a*sqrt(omega_dm_ini*pow(a0/a,3.0)+ (1.0-omega_dm_ini));
	  
	#pragma openmp parallel for private(j,k,ci,ind,c1,psi_k,psi_amp,poisson_rhs)
	 for(i=0;i<n[0];++i)
	 {
		  for(j=0;j<n[1];++j)
		  {
		   
		     for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;
			potn = phi.get_potential(ci);
			c1  = psi.get_psi(ind,psi_val[ci]);///Check if this as intented
			
			c1 = psi.calc_vel(ind,psi_vel, potn, a, a_t,dx);
			//printf("i %d j %d k %d psi %lf\n",i,j,k,psi_val[ci][0]);
			psi_k[0] = psi_val[ci][0]+ psi_vel[0]*dt;
			psi_k[1] = psi_val[ci][1]+ psi_vel[1]*dt;
			psi_amp = sqrt(psi_k[0]*psi_k[0] + psi_k[1]*psi_k[1]);
			
			psi.update(ind,psi_k[0],psi_k[1]);

			poisson_rhs = 1.5*H0*H0*a*a*(psi_amp*psi_amp - omega_dm_ini*pow(a0/a,3.0));
			phi.update_4pieGpsi(ci,poisson_rhs);
			
			//fail =1;

		    }

		   }
	}
	phi.solve_poisson(psi,k_grid);

	ak = a+a_t*dt;
	a_t = ak*sqrt(omega_dm_ini*pow(a0/ak,3.0)+ (1.0-omega_dm_ini));

	#pragma openmp parallel for private(j,k,ci,ind,c1,psi_k,psi_amp,poisson_rhs)
	for(i=0;(i<n[0])&&(!fail);++i)
	 {
		  for(j=0;(j<n[1])&&(!fail);++j)
		  {
		    for(k=0;(k<n[2])&&(!fail);++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			ind[0] = i;ind[1] = j;ind[2] = k;
			
			c1  = psi.get_psi(ind,psi_retrive);			
			c1 = psi.calc_vel(ind,psi_vel, potn, a, a_t,dx);
			psi_k[0] = 0.5*(psi_retrive[0]+psi_val[ci][0]+ psi_vel[0]*dt);
			psi_k[1] = 0.5*(psi_retrive[1]+psi_val[ci][1]+ psi_vel[1]*dt);
			psi.update(ind,psi_k[0],psi_k[1]);

			psi_amp = sqrt(psi_k[0]*psi_k[0] + psi_k[1]*psi_k[1]);
			if(isnan(psi_amp))
			{	fail=1;
				printf("FAIL %d %d %d\n",i,j,k);break;
			}

			poisson_rhs = 1.5*H0*H0*a*a*(psi_amp*psi_amp - omega_dm_ini*pow(a0/ak,3.0));
			phi.update_4pieGpsi(ci,poisson_rhs);
			


		    }

		   }
	}

	phi.solve_poisson(psi,k_grid);


	a = 0.5*(ak+a+a_t*dt);
	 if(isnan(a))
	  {printf("FAILED  \n"); fail = 1;	break;
	  }

	}

	a3a03omega = pow(a/a0,3.0)/omega_dm_ini;
	psi.write_psi(fp_psi,dx,a3a03omega,a,true, true);


	fclose(fp_psi);

	return(fail);

}


double ini_power_spec(double ksqr)
{
	return (1e-5);
}


double dlogD_dloga(double a)
{
	return (1.0);

}

void initialise(int * ind,fdm_psi &psi,metric_potential &phi,double k_grid[][3],double a0,double ai,double Hi,double omega_dm_ini,double *dx,ini_power_generator gen)
{
      

    
      int ci,i,j,k;
      int xcntr[3]={-1,-1,-1};
      
      double ktmp,maxkmagsqr = 0.0,minkmagsqr = 1e10;
     
      double a = ai;
      double a_t = a*Hi;
      double a_ti = ai*Hi;
     

      double L[3];
      double dk; 
      int kbins;     
      int tN = ind[0]*ind[1]*ind[2]; 
      int kbin_count[tN],kbin_grid[tN];
      double kmag_grid[tN];
      int n[3]{ind[0],ind[1],ind[2]};
      int loc_ind[3],err_hold;

      double ini_dc[tN],ini_theta[tN];
      double psi_amp,psi_r_val,psi_i_val,poisson_rhs;
       
	     
	 
	

	double kf = twopie*lenfac/(64.0);

	dx[0] = 0.02; dx[1] =0.02; dx[2] = 0.02;
        L[0] = dx[0]*((double) ind[0]);  L[1] = dx[1]*((double) (ind[1]));  L[2] = dx[2]*((double) (ind[2]));
	dk = 0.01/dx[0]; kbins = 0; printf("dk %lf\n",dk);
	
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
	printf("kmin %lf kmax %lf\n",sqrt(minkmagsqr),sqrt(maxkmagsqr));
	
	ini_rand_field(n,kmag_grid, ini_dc,ini_theta,ai,a0,a_ti,gen);

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

			poisson_rhs = 1.5*H0*H0*a*a*(psi_amp*psi_amp - omega_dm_ini*pow(a0/ai,3.0));

			//printf("P rhs %.10lf %.10lf %.10lf %.10lf %.10lf\n",poisson_rhs,psi_amp*psi_amp,psi_r_val,psi_i_val,ini_theta[ci]);
			phi.update_4pieGpsi(ci,poisson_rhs);
			


		    }
		  }

		}

	phi.solve_poisson(psi,k_grid);

/*	for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;
			loc_ind[0] = i;  loc_ind[1] = j;  loc_ind[2] = k;
			double potn = phi.get_potential(ci);

			printf("Checking ini potn %.13lf\n",potn);

	            }
		  }

		}
*/

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


void ini_rand_field(int * ind,double *kmag_grid,double * ini_dc,double * ini_theta_return,double a,double a0,double a_t,ini_power_generator gen)
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
			    //sigk  = sqrt(ini_power_spec(sqrt(ksqr)));
			if(ksqr>0.0)	
			    sigk  = gen.get_ini_spectrum(h*sqrt(ksqr));
			    muk = sigk/sqrt(2.0);
		 	    a1 = genrand_res53();
 			    a2 = genrand_res53(); 
			   // b1 = genrand_res53();
 			  //  b2 = genrand_res53();
			    a_rand = (muk*(sqrt(-2.0*log(a1))*cos(2.0*M_PI*a2)));
			    b_rand = (muk*(sqrt(-2.0*log(a1))*cos(2.0*M_PI*a2)));
				
			    F_ini_del[cnt][0] = a_rand;	F_ini_del[cnt][1] = b_rand;
			    if(ksqr>0.0)
			    { F_ini_theta[cnt][0] =  mass*(a_t/a)*f_ini*(a/a0)*(a/a0)*F_ini_del[cnt][0]/ksqr;
			      F_ini_theta[cnt][1] =  mass*(a_t/a)*f_ini*(a/a0)*(a/a0)*F_ini_del[cnt][1]/ksqr;
			    }	
			    else
			    { F_ini_theta[cnt][0] =  0.0;
			      F_ini_theta[cnt][1]  = 0.0;
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



