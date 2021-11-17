#include <fftw3.h>
#include "spline.c"


class scalar_field_3d
{

	protected:
	 double ***f,***f_t,****f_x,***f_lap;
         int n[3];


	public:

	scalar_field_3d(int *n_arr,bool need_lap=false,bool need_space_grads=false)
	{
	   int i,j;
	   n[0] = n_arr[0];	n[1] = n_arr[1];	n[2] = n_arr[2];
	   f = new double** [n[0]] ;
	   f_t = new double** [n[0]] ;
	   if(need_lap)
	   f_lap = new double** [n[0]] ;
	   if(need_space_grads)
	   { f_x = new double *** [3];
	     f_x[0] = new double** [n[0]];
	     f_x[1] = new double** [n[0]];
	     f_x[2] = new double** [n[0]];  
	   }

	   for(i=0;i<n[0];++i)
	   {
		f[i] = new double* [n[1]];
		f_t[i] = new double* [n[1]];
		if(need_lap)
		f_lap[i] = new double* [n[1]];
		if(need_space_grads)
		{ f_x[0][i] = new double* [n[1]];
		  f_x[1][i] = new double* [n[1]];
		 f_x[2][i] = new double* [n[1]];
		}
		for(int j=0;j<n[1];++j)
	     	{
		  f[i][j] = new  double[n[2]] ;
		  f_t[i][j] = new  double[n[2]] ;
		  if(need_lap)
		   f_lap[i][j] = new  double[n[2]] ;
		  if(need_space_grads)
		  { f_x[0][i][j] = new  double[n[2]] ;
		     f_x[1][i][j] = new  double[n[2]] ;
		     f_x[2][i][j] = new  double[n[2]] ;
		  }
	

		 }

	    }
	if((need_lap)||(need_space_grads))
	  {	
		if((need_lap)&&(need_space_grads))
 		 cout<<"Field allocated with arrays for space der and laplacian\n";
		else
		  if(need_lap)
			cout<<"Field allocated with array for laplacian\n";	
		  else
			cout<<"Field allocated with arrays for space der\n";

	 }
	else
	   cout<<"Field allocated withOUT arrays for space der and laplacian\n";
		

	}

	int cal_spt_grads(int *ind,double *dx,bool laplacian=false,bool spt_grad=false)
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
				default:
					return(f[ind[0]][ind[1]][ind[2]]);

			}
		
	}

};

class fdm_psi
{
	private:
	scalar_field_3d psi_r;
	scalar_field_3d psi_i;
	int n[3];
	
	public:  
	//int n[3];
	fdm_psi(int *ind,bool lb=false,bool sgb=false):psi_r(ind,lb,sgb),psi_i(ind,lb,sgb)
	{
		n[0] = ind[0];  n[1] = ind[1];  n[2] = ind[2];
	}
	
	

	int calc_vel(int * ind,double *v,double potn,double a,double a_t,double *dx)
	{
		int c1;		
		double psi_r_lap,psi_i_lap,psi_r_val,psi_i_val;
		psi_r_val = psi_r.get_field(ind,give_f);
		psi_i_val = psi_i.get_field(ind,give_f);
		c1 = psi_r.cal_spt_grads(ind,dx,true);
		c1 = psi_i.cal_spt_grads(ind,dx,true);
		psi_r_lap = psi_r.get_field(ind,give_f_lap);
		psi_i_lap = psi_i.get_field(ind,give_f_lap);
		v[0] = -1.5*(a_t/a)*psi_i_val;
		v[0]+= (-0.5*hbar_by_m*psi_i_lap/(a*a) + potn*psi_i_val/hbar_by_m);
		v[1] = -1.5*(a_t/a)*psi_r_val;
		v[1]+= (0.5*hbar_by_m*psi_r_lap/(a*a) - potn*psi_r_val/hbar_by_m);

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

	int get_psi(int *ind,double *psi_ret,code1 c = give_f)
	{
		int c1=1;		
		
		psi_ret[0]= psi_r.get_field(ind,c);
		psi_ret[1] = psi_i.get_field(ind,c);

		return c1;




	}


	void write_psi(FILE *fp_psi,double *dx,double a3a03omega,double a,bool get_dc=false, bool get_psi=true)
	{
		int i,j,k,locind[3];
		double psi_r_val,psi_i_val,psi_amp2,dc;
		for(i=0;i<n[0];++i)
		{
			for(j=0;j<n[1];++j)
			{
				for(k=0;k<n[2];++k)
				{
					locind[0]=i; locind[1]=j; locind[2]=k;
					psi_r_val=psi_r.get_field(locind,give_f);
					psi_i_val=psi_i.get_field(locind,give_f);	
							

					if(get_dc)
					{
					  psi_amp2 = psi_r_val*psi_r_val + psi_i_val*psi_i_val;		
					  dc= a3a03omega*psi_amp2 - 1.0;
					  fprintf(fp_psi,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",a,dx[0]*i,dx[1]*j,dx[2]*k,psi_r_val,psi_i_val,dc);
					}

					else
					{
					  fprintf(fp_psi,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",a,dx[0]*i,dx[1]*j,dx[2]*k,psi_r_val,psi_i_val);


					}


				}
	
			}

		}
		
		fprintf(fp_psi,"\n\n\n\n");

	}



};


class metric_potential
{

	private:
	int n[3];
	//scalar_field_3d phi;

	fftw_complex *fpGpsi;
	fftw_complex *fpGpsi_ft;

	fftw_plan plan_pois_f;
	fftw_plan plan_pois_b;
	
	public:
	
	metric_potential(int *ind,bool lb=false,bool sgb=false)//:phi(ind,lb,sgb)
	{
		int l = ind[0]*ind[1]*ind[2];
		n[0]=ind[0];n[1]=ind[1];n[2]=ind[2];
		
		fpGpsi = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * l);
		fpGpsi_ft =(fftw_complex*) fftw_malloc(sizeof(fftw_complex) *l);


		plan_pois_f = fftw_plan_dft_3d(ind[0],ind[1],ind[2], fpGpsi, fpGpsi_ft, FFTW_FORWARD, FFTW_ESTIMATE);
		plan_pois_b = fftw_plan_dft_3d(ind[0],ind[1],ind[2],fpGpsi_ft, fpGpsi, FFTW_BACKWARD, FFTW_ESTIMATE);

	}


	void solve_poisson(fdm_psi psi,double k_grid[][3])
	{
		int i,j,k,ci,ind[3]{0,0,0},r;
		double k2fac;
		fftw_execute(plan_pois_f);
		double sqrt_tN = sqrt((double)(n[0]*n[1]*n[2])); 

		for(i=0;i<n[0];++i)
		{
		  for(j=0;j<n[1];++j)
		  {
		    for(k=0;k<n[2];++k)
		    {
			ci = (n[2]*n[1])*i + n[2]*j + k;			
			k2fac = twopie*twopie*(k_grid[ci][0]*k_grid[ci][0]+k_grid[ci][1]*k_grid[ci][1]+k_grid[ci][2]*k_grid[ci][2]);
			
			if(k2fac>0.0)
			{fpGpsi_ft[ci][0] = fpGpsi_ft[ci][0]/(k2fac*sqrt_tN*sqrt_tN);
			 fpGpsi_ft[ci][1] = fpGpsi_ft[ci][1]/(k2fac*sqrt_tN*sqrt_tN);
			}	
			else
			{fpGpsi_ft[ci][0] = 0.0;
			 fpGpsi_ft[ci][1] = 0.0;
			}

			

		    }

		  }

		}
		
		fftw_execute(plan_pois_b);
	

	}


	void update_4pieGpsi(int ci,double val)
	{
		
		fpGpsi[ci][0] = val;
		fpGpsi[ci][1] = 0.0;
	
	}


	double get_potential(int ci)
	{

		return (fpGpsi[ci][0]);	

	}


/*	int update(int * ind,double phi_val)
	{
		int c1;
		c1 = phi.update_field(ind,phi_val);
		
		return (c1);
	}

*/


};


class ini_power_generator
{

	private:
	double *p,*k;
	const char *fname;
	double **c,max_dx=0.0,x_min;
	int point_cnt,checkspl;
	
	
	public:

	ini_power_generator(const char *name)
	{
		double a,b,pre_x=0.0;
		int i = 0,fre=1;
		fname =  name;		

		FILE *fp = fopen(fname,"r");
		while(fre!=EOF)
		{
			fre=fscanf(fp,"%lf\t%lf\n",&a,&b);

			if(i==0)
			{
			
				max_dx = a;	
				x_min = a;
				pre_x = a;	


			}

			else
			if(fre!=EOF)
			{
			   if(fabs(a-pre_x)>max_dx)
				max_dx = fabs(a-pre_x);
			   if(a<x_min)
				x_min = a;
			   pre_x = a;

			}
			
			if(fre!=EOF)
			{
			// printf("Reading %d %d %lf %.8lf\n",i,fre,a,b);
			 ++i;
			}

			
		
			
		}

		fclose(fp);
		point_cnt = i;
		p = new double[point_cnt];	
		k = new double[point_cnt];
		c = new double* [point_cnt];	
		for(i=0;i<point_cnt;++i)
			c[i] =  new double[3];	
		fp = fopen(fname,"r");
		i=0;
		fre=1;
		double cload[point_cnt][3];
		while(fre!=EOF)
		{
			fre=fscanf(fp,"%lf\t%lf\n",&a,&b);
			if(fre!=EOF)
			{*(k+i)=a;
			 *(p+i)=b;
			 //printf("ff %d %lf %.10lf\n",i,*(k+i),*(p+i));
			}
			++i;
		}
		
		fclose(fp);
		
		checkspl = spline(k,p,point_cnt,cload); 
	
		for(i=0;i<point_cnt;++i)
		{
			c[i][0] = cload[i][0];
			c[i][1] = cload[i][1];
			c[i][2] = cload[i][2];
			 

		}

		
		if(checkspl)
		printf("\nALERT !!!! Spline fit failed....%.10lf\n",max_dx);
	
		

	}

	double get_ini_spectrum(double x)
	{

	  double Val;
	  int cur_ind = (int)((x-x_min)/max_dx);
	  int ind_match=0;//printf("Yo %d %lf %lf\n",cur_ind,x,max_dx);
	while(cur_ind<point_cnt)
	 { if((x>=k[cur_ind])&&(x<=k[cur_ind+1]))
    	  {//printf("cur_ind is %d\n",cur_ind);
            ind_match = 1;
            Val = (p[cur_ind]+c[cur_ind][0]*(x-k[cur_ind])+c[cur_ind][1]*(x-k[cur_ind])*(x-k[cur_ind])
                                                  +c[cur_ind][2]*(x-k[cur_ind])*(x-k[cur_ind])*(x-k[cur_ind]));
          	break;
          }
	  ++cur_ind;
	  
	 }

	if(!ind_match)
	printf("\nALERT index overrun for spline...for x %lf\n",x);

	return Val;
		

	}

	void check(double x)
	{
		int i;
		double v;
		FILE *fpcheck = fopen("spline_check.txt","w");
		for(i=0;i<point_cnt;++i)
		{
			v = get_ini_spectrum(*(k+i));
			fprintf(fpcheck,"%.10lf\t%.10lf\t%.10lf\n",*(k+i),v,*(p+i));	
			
			//printf("Check result is v %lf @ x %lf\n",v,x);
		}
		
		
		

	}



};

