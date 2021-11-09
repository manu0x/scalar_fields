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

class scalar_field_3d
{

	protected:
	 double ***f,***f_t,****f_x,***f_lap;
         int n[3];


	public:

	scalar_field_3d(int *n_arr)
	{
	   int i,j;
	   n[0] = n_arr[0];	n[1] = n_arr[1];	n[2] = n_arr[2];
	   f = new double** [n[0]] ;
	   f_t = new double** [n[0]] ;
	   f_lap = new double** [n[0]] ;
	   f_x = new double *** [3];
	   f_x[0] = new double** [n[0]];
	   f_x[1] = new double** [n[0]];
	   f_x[2] = new double** [n[0]];  

	   for(i=0;i<n[0];++i)
	   {
		f[i] = new double* [n[1]];
		f_t[i] = new double* [n[1]];
		f_lap[i] = new double* [n[1]];

		f_x[0][i] = new double* [n[1]];
		f_x[1][i] = new double* [n[1]];
		f_x[2][i] = new double* [n[1]];

		for(int j=0;j<n[1];++j)
	     	{
		  f[i][j] = new  double[n[2]] ;
		  f_t[i][j] = new  double[n[2]] ;
		  f_lap[i][j] = new  double[n[2]] ;
		  f_x[0][i][j] = new  double[n[2]] ;
		  f_x[1][i][j] = new  double[n[2]] ;
		  f_x[2][i][j] = new  double[n[2]] ;
	

		 }

	    }
	cout<<"Field allocated\n";
	}

	int cal_lap(int *ind,double *dx,bool spt_grad=true,bool laplacian=true)
	{
	   int i,j;
	   double m[3],lapsum=0.0;
	   int ind_l1[3]{ind[0],ind[1],ind[2]},ind_l2[3]{ind[0],ind[1],ind[2]},ind_r1[3]{ind[0],ind[1],ind[2]},ind_r2[3]{ind[0],ind[1],ind[2]};
	   for(i=0;i<3;++i)
	   {
	     
	     
	     ind_l1[i] = (n[i]+ind[i]-1)%n[i];
	     ind_l2[i] = (n[i]+ind[i]-2)%n[i];

	     ind_r1[i] = (ind[i]+1)%n[i];
	     ind_r2[i] = (ind[i]+2)%n[i];

	    
	     if(spt_grad==true)
	     {
	       m[0] = (-8.0*f[ind_l1[0]][ind_l1[1]][ind_l1[2]]+8.0*f[ind_r1[0]][ind_r1[1]][ind_r1[2]]); 
	       m[1] = (f[ind_l2[0]][ind_l2[1]][ind_l2[2]]-f[ind_r2[0]][ind_r2[1]][ind_r2[2]]);
	       m[2] = m[0] + m[1] ;
		
		
	       f_x[i][ind[0]][ind[1]][ind[2]] = m[2]/(dx[i]);
	      }

	     if(laplacian==true)
	      {	m[0] = (16.0*f[ind_l1[0]][ind_l1[1]][ind_l1[2]]+16.0*f[ind_r1[0]][ind_r1[1]][ind_r1[2]]); 
	      	m[1] = (-f[ind_l1[0]][ind_l1[1]][ind_l1[2]]-f[ind_r2[0]][ind_r2[1]][ind_r2[2]]);
	      	m[2] = m[0] + m[1] -30.0*f[ind[0]][ind[1]][ind[2]];
	     	lapsum+= (m[2]/(dx[i]));

	     	
	      }
		ind_l1[i] = ind[i];
	    	ind_l2[i] = ind[i];

	     	ind_r1[i] = ind[i];
	     	ind_r2[i] = ind[i];

	   } 

	if(laplacian==true)
	  f_lap[ind[0]][ind[1]][ind[2]]=lapsum;

	  if(isnan(lapsum))
		return(0);
	  else
		return (1);
	
	}

	int update_field(int * ind,double fu,double f_tu=0.0)
	{
		f[ind[0]][ind[1]][ind[2]] = fu;
		f_t[ind[0]][ind[1]][ind[2]] = f_tu;

		if(isnan(fu+f_tu))
			return(0);
		else
			return(1);

	}

	double get_field(int *ind,code1 c)
	{
		
		switch(c){
				case give_f:
					return(f[ind[0]][ind[1]][ind[2]]);
				case give_f_t:
					cout<<"Reddd\n";
					return(f_t[ind[0]][ind[1]][ind[2]]);
				case give_f_x:
					cout<<"Green\n";
					return(f_x[0][ind[0]][ind[1]][ind[2]]);
				case give_f_y:
					return(f_x[1][ind[0]][ind[1]][ind[2]]);
				case give_f_z:
					return(f_x[2][ind[0]][ind[1]][ind[2]]);
				case give_f_lap:
					return(f_lap[ind[0]][ind[1]][ind[2]]);
					

			}
		
	}

};

class fdm_psi
{
	private:
	scalar_field_3d psi_r;
	scalar_field_3d psi_i;
	
	public:
	//int n[3];
	fdm_psi(int *ind):psi_r(ind),psi_i(ind)
	{
	}
	
	

	int calc_vel(int * ind,double *v,double potn,double a,double a_t)
	{
		double psi_r_lap,psi_i_lap,psi_r_val,psi_i_val;
		psi_r_val = psi_r.get_field(ind,give_f);
		psi_i_val = psi_i.get_field(ind,give_f);
		psi_r_lap = psi_r.get_field(ind,give_f_lap);
		psi_i_lap = psi_i.get_field(ind,give_f_lap);
		v[0] = -1.5*(a_t/a)*psi_i_val;
		v[0]+= (-0.5*hbar_by_m*psi_i_lap/(a*a) + potn/hbar_by_m);
		v[1] = -1.5*(a_t/a)*psi_r_val;
		v[1]+= (0.5*hbar_by_m*psi_r_lap/(a*a) - potn/hbar_by_m);

		if(isnan(v[0]+v[1]))
			return (-1);
		else
			return(1);
	}

	double cal_mod(int * ind)
	{
		double psi_r_val,psi_i_val,mag;		
		
		psi_r_val = psi_r.get_field(ind,give_f);
		psi_i_val = psi_i.get_field(ind,give_f);

		
		mag = sqrt(psi_r_val*psi_r_val+psi_r_val*psi_r_val);

		return(mag);

	}

	int update(int * ind,double psir,double psii)
	{
		int c1,c2;
		c1 = psi_r.update_field(ind,psir);
		c2 = psi_i.update_field(ind,psii);
		
		return (c1*c2);
	}



};


class metric_potential
{

	private:
	int n[3];
	scalar_field_3d phi;

	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	
	public:
	
	metric_potential(int *ind):phi(ind)
	{
		int l = ind[0]*ind[1]*ind[2];
		n[0]=ind[0];n[1]=ind[1];n[2]=ind[2];
		
		fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
		fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *l);


		plan_pois_f = fftw_plan_dft_3d(ind[0],ind[1],ind[2], fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
		plan_pois_b = fftw_plan_dft_3d(ind[0],ind[1],ind[2],fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);

	}


	void solve_poisson(fdm_psi psi,double **k_grid)
	{
		int i,j,k,ci,ind[3]{0,0,0},r;
		double k2fac;
		fftw_execute(plan_pois_f);

		for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;			
			k2fac = twopie*twopie*(k_grid[ci][0]*k_grid[ci][0]+k_grid[ci][1]*k_grid[ci][1]+k_grid[ci][2]*k_grid[ci][2]);
			fpGpsi_ft[ci][0] = 2.0*twopie*G*mass*fpGpsi_ft[ci][0]/k2fac;
			fpGpsi_ft[ci][1] = 2.0*twopie*G*mass*fpGpsi_ft[ci][1]/k2fac;

		    }

		  }

		}
		
		fftw_execute(plan_pois_f);
		

	}



};



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

void initialise(int * ind,fdm_psi * psi)
{
      int l1,l2,r1,r2;
      double Vvlb;

    
      int px,py,pz,ci,pgi,j;
      int xcntr[3]={-1,-1,-1},anchor[3];
      double gamma, v, gradmagf;
      double ktmp,maxkmagsqr = 0.0,minkmagsqr = 1e10;
      double a0 = 1.00;
      double ai = 0.001;
      double a = ai;
      double Hi = 1.0;//H0;
      double a_t = a*Hi;
      double a_ti = ai*Hi;
      double cpmc = 0.3;
      double omega_dm_ini= (cpmc)*pow((a0/ai),3.0)/(cpmc*a0*a0*a0/(ai*ai*ai) + (1.0-cpmc));
 
      double L[3],dx[3];
      double dk; 
      int kbins;     
      double lenfac = 1.0;
      int tN = ind[0]*ind[1]*ind[2];
      int kbin_count[tN],kmag_grid[tN];
      double k_grid[tN][3];
      int n[3]{ind[0],ind[1],ind[2]};
      int loc_ind[3];

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
		

		kmag_grid[ci] = (int)(sqrt(ktmp)/(dk));
		 //printf("yo  %d  %lf\n",kmag_grid[ci],sqrt(ktmp));
		++kbin_count[kmag_grid[ci]];

		if(kmag_grid[ci]>kbins)
		kbins=kmag_grid[ci];
		
			

		
		
		
      	}


	ini_rand_field(n,kmag_grid, ini_dc,ini_theta,ai,a0,a_ti)

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
		

		fprintf(fpinirand,"%d\t%.16lf\t%.16lf\n",
					cnt,ini_del[cnt][0],ini_theta[cnt][0]);


		

	}
    

	 fftw_free(F_ini_del);
	 fftw_free(ini_del);
	 fftw_destroy_plan(ini_del_plan);
	
	 fftw_free(F_ini_theta);
	 fftw_free(ini_theta);
	 fftw_destroy_plan(ini_theta_plan);
	

	
}



