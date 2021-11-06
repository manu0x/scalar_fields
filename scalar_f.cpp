#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

#include <fftw3.h>

using namespace std;

//////////////////////////////Constants/////////////////////////
#define twopie  2.0*M_PI
#define mass 1.0
#define G 1.0

///////////////////////////////////////////Todo List///////////////////////////////////////////
//////////////   1) Add error check for successful memory allocation in class scalar_field_3d    /////////////////////////////////////////

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
		return(-1);
	  else
		return (1);
	
	}

	int update_field(int * ind,double fu,double f_tu)
	{
		f[ind[0]][ind[1]][ind[2]] = fu;
		f_t[ind[0]][ind[1]][ind[2]] = f_tu;

		if(isnan(fu+f_tu))
			return(-1);
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
