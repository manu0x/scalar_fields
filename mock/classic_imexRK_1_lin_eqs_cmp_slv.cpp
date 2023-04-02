using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <limits>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>

#define pie M_PI 



//////////GLobal constants/////////////

double m,n,T;

///////////////////////////////////////


//////////// ImEx RK Butcher Tableau /////

const int imex_s = 3;
const int im_s = imex_s;
const int ex_s = imex_s;

double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  205.0/462.0,2.0/11.0,0.0,  2033.0/4620.0,21.0/110.0,2.0/11.0};
double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  5.0/6.0,0.0,0.0,  11.0/24.0,11.0/24.0,0.0};

double im_c[im_s] = {2.0/11.0,289.0/462.0,751.0/924.0};
double ex_c[ex_s] = {0.0,5.0/6.0,11.0/12.0};

double im_b[im_s] = {24.0/55.0,1.0/5.0,4.0/11.0};
double ex_b[ex_s] = {24.0/55.0,1.0/5.0,4.0/11.0};


//	/	/	/	/	/	/	/	/	/	/	
/*
double im_a[im_s][im_s] = {2.0/11.0,0.0,0.0,  41.0/154.0,2.0/11.0,0.0,  289.0/847.0,42.0/121.0,2.0/11.0};
double ex_a[ex_s][ex_s] = {0.0,0.0,0.0,  0.5,0.0,0.0,  0.5,0.5,0.0};

double im_c[im_s] = {2.0/11.0,69.0/154.0,67.0/77.0};
double ex_c[ex_s] = {0.0,0.5,1.0};

double im_b[im_s] = {1.0/3.0,1.0/3.0,1.0/3.0};
double ex_b[ex_s] = {1.0/3.0,1.0/3.0,1.0/3.0};
*/

////////////////////////////////////////////



void d_xx(double *p,double *der,double dx,int N)
{

	int i,j,l1,r1;
	double vl1,vr1,vc;
   for(i=0;i<(N);++i)
   {	l1 = ((N)+ (i-1))%(N);
	

	r1 = (i+1)%(N);
	


	vl1 = p[l1];
	

	vr1 = p[r1];
	

	vc = p[i];
	
	der[i] =  (vr1 +vl1 - 2.0*vc )/(dx*dx);

	

	

    }		
	



}

////////////////////  Function for calculating potential //////////////////////////////////////////////////////////////

double V(double psi_amp)
{
	//return(0.0);
	return( 1.0 - 2.0*pie*pie*n*n/(m*m));


}

///////////////////////////  Explicit term for veolcity //////////////////////////////////////////////////////////////////
void ex_vel(double v[2],double psi[2],double Vval)
{

	v[0] = m*Vval*psi[1];
	v[1] = -m*Vval*psi[0];;



}

/////////////////////////// Function to create and invert matrix (1+w^2P^2) stored in Mp and matrix Inv((1+w^2P^2))P ////////

void cr_invert_mat(gsl_matrix_complex *Mp,gsl_permutation * p,int N, double dt, double dx, double diff,double gamma)
{	///////////////////// This assumed Diagonal with same element gamma
	




}






//////////////////////////////// Function to set initial conditions /////////////////////////////////////////////////

void initialise(double *psi,double *x,double *k,double dx,int N)
{

	int i,j; double di,dk;
	dk = 2.0*pie/(((double)(N))*dx);
	for(i=0;i<N;++i)
	{
	   di = (double)i;	
	   *(psi+i*2) = sin(2.0*pie*n*di*dx);
	   *(psi+i*2+1) = 0.0;

	    if(i<=(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk; 

		x[i] = 0.0+di*dx;

	//	printf("k %lf  %lf  dx %lf\n",dk,*(k+i),dx );

	}


}

////////////////////////////////////  Function to evolve the system given dt and dx. Returns energy errors //////////////////////////// 
			///// stb_any = 1 does the evoln. only for short duration ///////
double run(double dt,double dx,double *abs_err,double *stb_avg,int stb_any,int printfp,int prt)
{

	int N,t_steps;
	double box_len,t_end,t_start,diff,stb_ini;


///////////////////////////////Files for data/////////////////////////////////////////////////////

	FILE *fp = fopen("data_linq.txt","w");
	FILE *fp2 = fopen("data2_linq.txt","w");
	FILE *fptemp = fopen("temp_linq.txt","w");
	
	FILE *fpmass = fopen("mass_linq.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	m  = 1.0;
	n  = 1.0;
	T = 2.0*pie/m;
	diff = 0.0;
	N = 100;





/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 2.0;
	//dx = box_len/(double(N));
	N = ((int)(box_len/dx)) + 1;

////////////////////////////// Time & dt settings ////////////////////////////////////////////////
	
	t_start = 0.0;
	t_end = 2.0*T;
	if(stb_any)
	t_end = 0.01*t_end;
	t_steps = (int)((t_end-t_start)/dt);
	//dt  = (t_end-t_start)/((double)t_steps);
	printf("dt %lf N %d\n",dt,N);
	
/////////////////////////////////////////RK things/////////////////////////////////////////


	double im_K_P[ex_s][N][2],ex_K_P[im_s][N][2];

////////////////////////////// Psi variables  /////////////////////////////////////////////////////

	///////////// Psi stores calculated solution, Psik stores intermediate solutions at substage and Psib stores rhs of linear system /////
	////i.e.////////////     Psik = M(I+iwP)Psib	/// where M defined as Inverse(I+w^2^P^2)

	double Psi[N][2],lap_val[2*N],lambda;

	
	double x[N],k_grid[N];
	
	double Psib[2],im_K[im_s][N][2],ex_K[ex_s][N][2],dP_xx[N];
	double vel_val[2],c_psi[2],c_psi_amp,Vval;

	int i,j;


///////////////	 Initialise the system	///////////////////////////////////////////////////////////

	initialise(&Psi[0][0],x,k_grid,dx,N);

	
	
	 

///// Get inverted matrix ///////////////////////////


  size_t tn  = N;

  gsl_matrix_complex * MLU = gsl_matrix_complex_alloc(tn,tn);
  gsl_permutation * prm = gsl_permutation_alloc (tn);

  gsl_matrix_complex_set_zero(MLU);


	
	double mu = dt/(dx*dx);
	
	double omega = dt*im_a[0][0]/(dx*dx*2.0*m);
	
	gsl_complex cmp;
	
	
	
	
	int sig, out;
	
	double P[N][N];
	
	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {
	    
	     P[i][j] = 0.0;
		
	   }
	}


	for(i=1;i<N-1;++i)
	{
	  for(j=0;j<N;++j)
	  {
		//M[i][j] = 0.0;
		if(i==j)
		{
		  P[i][j] = -2.0;	
		  P[i][j-1] = 1.0;
		  P[i][j+1] = 1.0;

		
		}
		

	  }



	}

	P[0][N-1] = 1.0;
	P[N-1][0] = 1.0;
	P[0][1] = 1.0;
	P[N-1][N-2] = 1.0;



	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {
		//M[i][j] = 0.0;
		if(i==j)
		{ 
		  GSL_REAL(cmp) = 1.0;
		  
		}
		else
		 GSL_REAL(cmp) = 0.0;

	
		GSL_IMAG(cmp) = -omega*P[i][j]; 

		  gsl_matrix_complex_set(MLU, i, j,cmp);
		  

		

	  }

	

	}




////////////////////// LU decomp and inversion ////////////////////////////////////////////////////

	out = gsl_linalg_complex_LU_decomp(MLU, prm, &sig);
	











///////////////////////  Evolution ///////////////////////////////////////////////////////////////
	srand(time(0));

	int s_cntr,tcntr,printcntr,fail=0;
	int l1,r1;

	printcntr = (int)(((double) t_steps)/100.0);	//printf("%d\n",printcntr);
	if(t_steps<=100)
		printcntr = 1.0;
	double t,vl1,vr1,vc,amp,avg_amp,sol[2],dbi;
	double drc = 2.0*pie*n*2.0*pie*n;
	double fdt,amp_ini;

	
	*stb_avg = 0.0;

	gsl_vector_complex * Psib_cmp = gsl_vector_complex_alloc(tn);
	gsl_vector_complex * Psik = gsl_vector_complex_alloc(tn);
	
	for(i=0;i<N;++i)
	{
		GSL_REAL(cmp) = Psi[i][0];
		GSL_IMAG(cmp) = Psi[i][1];
	
		gsl_vector_complex_set (Psib_cmp, i, cmp);

	}
	
	gsl_linalg_complex_LU_solve(MLU,prm, Psib_cmp,Psik);


	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(tcntr<1);t+=dt,++tcntr)
	{	
		avg_amp = 0.0;
	
	 	for(i=0;i<N;++i)
		{	
			
			dbi = (double)(i);
			//////////////////////////////	Printing data at some regular steps   ///////////////////////////
			if((tcntr%printcntr)==0||(stb_any))  
			{
			  if(i<N)
			  {amp  = sqrt(Psi[i][0]*Psi[i][0] + Psi[i][1]*Psi[i][1]);				
			   avg_amp+=amp;
			  }
			  sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			  sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			 if(printfp)
			  fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,Psi[i][0],Psi[i][1],sol[0],sol[1],amp);
			  if((i==20)&&(printfp))
			   fprintf(fptemp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,Psi[i][0],Psi[i][1],sol[0],sol[1],amp);
			}
			
		
			

			
		}

		stb_ini = avg_amp;
		
		
		if(t==t_start)
		{amp_ini = avg_amp;
		 //printf("av ini %lf g %lf\n",amp_ini,avg_amp);
		 if(prt)
		  printf("mu %lf dt %lf\tt %lf\tamp_ini %lf\n\n",dt/(dx*dx),dt,t,amp_ini);
		}

		for(s_cntr=1;s_cntr<imex_s;++s_cntr)
		{

			 for(i=0;i<N;++i)
			 {
			    for(j=0;j<s_cntr;++j)
			    {
			  	
				if(j==0)
				{
					// Calc. explicit contribution from Psik of last substage //////////////////////
					if(s_cntr==1)
					  {c_psi[0] = Psi[i][0];
					   c_psi[1] = Psi[i][1];
					   c_psi_amp = sqrt(c_psi[0]*c_psi[0] + c_psi[1]*c_psi[1]);

					  }
					else
					  {cmp = gsl_vector_complex_get(Psik,i);
					   c_psi_amp=gsl_complex_abs(cmp);	
		
					  }
		    			
		    			Vval = V(c_psi_amp);
		    			ex_vel(vel_val,c_psi,Vval);
					

		    			ex_K_P[s_cntr-1][i][0] = vel_val[0];
					ex_K_P[s_cntr-1][i][1] = vel_val[1]; 

					//////////////////////////////////////////////////////////////////////////////////
					/////////	Calc spatial derivatives using last Psik	/////////////////

					l1 = i-1;
					r1 = (i+1);
					if(i==0)
					 l1 = N-2;
					if(i==(N-1))	
			 		 r1 = 1;
	
					vl1 = GSL_REAL(gsl_vector_complex_get(Psik,l1));
					vr1 = GSL_REAL(gsl_vector_complex_get(Psik,r1));	
					vc = GSL_REAL(gsl_vector_complex_get(Psik,i));
	
					im_K_P[s_cntr-1][i][1]  = (vr1 +vl1 - 2.0*vc )/(2.0*m*dx*dx);

					vl1 = GSL_IMAG(gsl_vector_complex_get(Psik,l1));
					vr1 = GSL_IMAG(gsl_vector_complex_get(Psik,r1));	
					vc = GSL_IMAG(gsl_vector_complex_get(Psik,i));
	
					im_K_P[s_cntr-1][i][0]  =  -(vr1 +vl1 - 2.0*vc )/(2.0*m*dx*dx);
					//printf("%d %d %lf\n",s_cntr,i,im_K_P[s_cntr-1][i] );
					/////////////////////////////////////////////////////////////////////////////////////////
					/////////////// Start addition to rhs from previous substeps ///////////////////////////

					Psib[0] = Psi[i][0] + dt*ex_a[s_cntr][j]*ex_K_P[j][i][0]+ dt*im_a[s_cntr][j]*im_K_P[j][i][0];
					Psib[1] = Psi[i][1] + dt*ex_a[s_cntr][j]*ex_K_P[j][i][1]+ dt*im_a[s_cntr][j]*im_K_P[j][i][1];
					

				}
			
				
				else
				{ Psib[0]+=  dt*ex_a[s_cntr][j]*ex_K_P[j][i][0]+ dt*im_a[s_cntr][j]*im_K_P[j][i][0];
				  Psib[1]+=  dt*ex_a[s_cntr][j]*ex_K_P[j][i][1]+ dt*im_a[s_cntr][j]*im_K_P[j][i][1];
				}
				
	
	
			   }

			   GSL_REAL(cmp) = Psib[0];
			   GSL_IMAG(cmp) = Psib[1];
			   gsl_vector_complex_set(Psib_cmp, i, cmp);
	
			}

			
			gsl_linalg_complex_LU_solve(MLU,prm, Psib_cmp,Psik);
				
			


		}//RK stages concluded
	
		avg_amp=0.0;

	
		*abs_err = 0.0;

		for(i=0;i<N;++i)
		{	
	
			cmp = gsl_vector_complex_get(Psik,i);
			c_psi_amp=gsl_complex_abs(cmp);	 

    			Vval = V(c_psi_amp);
    			ex_vel(vel_val,c_psi,Vval);
					

    			ex_K_P[imex_s-1][i][0] = vel_val[0];
			ex_K_P[imex_s-1][i][1] = vel_val[1]; 

			l1 = i-1;
			r1 = (i+1);
			if(i==0)
			 l1 = N-2;
			if(i==(N-1))	
			 r1 = 1;
	
			vl1 = GSL_REAL(gsl_vector_complex_get(Psik,l1));
			vr1 = GSL_REAL(gsl_vector_complex_get(Psik,r1));	
			vc = GSL_REAL(gsl_vector_complex_get(Psik,i));
	
			im_K_P[imex_s-1][i][1]  = (vr1 +vl1 - 2.0*vc )/(2.0*m*dx*dx);

			vl1 = GSL_IMAG(gsl_vector_complex_get(Psik,l1));
			vr1 = GSL_IMAG(gsl_vector_complex_get(Psik,r1));	
			vc  = GSL_IMAG(gsl_vector_complex_get(Psik,i));
	
			im_K_P[imex_s-1][i][0]  =  -(vr1 +vl1 - 2.0*vc )/(2.0*m*dx*dx);

				
			   for(j=0;j<imex_s;++j)
			   {	Psi[i][0]+=  dt*ex_b[j]*ex_K_P[j][i][0]+ dt*im_b[j]*im_K_P[j][i][0];
				
				Psi[i][1]+=  dt*ex_b[j]*ex_K_P[j][i][1]+ dt*im_b[j]*im_K_P[j][i][1];
			
				//if(i==51)
				//printf("!!!+tivity broken!!!j  %d\ti %d\t%lf\t%lf\n",j,i,P[i],im_K_P[j][i]);
			    }

			////////////// Check to stop if something blows up ////////////
				//if(isnan(P[i]))				
			if(isnan(Psi[i][0])||isnan(Psi[i][1]))
			 {

				if(isnan(Psi[i][0]))
				printf("!!!GONE Nan: Real!!!!  %d\tt %lf\t%lf\n",tcntr,t,Psi[i][0]);
				if(isnan(Psi[i][1]))
				printf("!!!GONE Nan: Imag!!!!  %d\tt %lf\t%lf\n",tcntr,t,Psi[i][1]);
					
				fail = 1;

				break;
			}
		


			if(fail)
			break;

			GSL_REAL(cmp) = Psi[i][0];
			GSL_IMAG(cmp) = Psi[i][1];
	
			gsl_vector_complex_set (Psib_cmp, i, cmp);
			


			amp  = sqrt(Psi[i][0]*Psi[i][0] + Psi[i][1]*Psi[i][1]);
			avg_amp+=amp;
			//printf("avg_amp %lf\n",avg_amp);
	
			sol[0] = cos(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi);
			sol[1] = -sin(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi);
	
			(*abs_err)+=(fabs(sol[0]-Psi[i][0])+fabs(sol[1]-Psi[i][1]));

			if((tcntr%printcntr)==0)  
			{ 
			  if(printfp)
			  fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,Psi[i][0],Psi[i][1],sol[0],sol[1],amp);
			}
			

			

		   }

		*abs_err = (*abs_err)/((double)N);
		(*stb_avg)+=(avg_amp/stb_ini);

		if((stb_any)&&(fail))
		{

			*abs_err = 1e5;
			*stb_avg = std::numeric_limits<double>::infinity();		
			return(1e5);
			
		}		  

		else//if(fail)
		if(((100.0*fabs(avg_amp-amp_ini)/amp_ini)>=1e3)||(fail)||((*abs_err)>=1e3))
		{		if(prt)
				printf("Overload %lf %lf %lf\n",((100.0*fabs(avg_amp-amp_ini)/amp_ini)),avg_amp,*abs_err);
			
				*abs_err = -1e3;
					
				return(-1e3);
			//// If energy or absolute errs are too huge, stop evolution.
	
		
		}
	




		if((tcntr%printcntr)==0)  
		{ 
		 if(printfp)
		  {fprintf(fp,"\n\n\n");
		   fprintf(fp2,"\n\n\n");
		   fprintf(fpmass,"%lf\n",t/t_end);	
		  }
		if(prt)
		  printf("mu %lf dt %lf\tt %lf\tamp %lf\n\n",dt/(dx*dx),dt,t,avg_amp);
		}
		


		  

	}///// ENd f Time Evolution /////////////////////////
	
	*stb_avg = (*stb_avg)/((double) tcntr);
	if(prt)
	printf("Error in conserv. %lf\n mu is %lf", 100.0*fabs(avg_amp-amp_ini)/amp_ini,dt/(dx*dx));

	fclose(fp);
	fclose(fp2);
	fclose(fptemp);
	fclose(fpmass);

	return(100.0*fabs(avg_amp-amp_ini)/amp_ini);


}



int main()
{
	double dt = 3e-4;
	double dx = 4e-3;
	double abs_err,en_loss,stb_avg;

	//double dx_l=4e-3,dx_u = 4e-2;
	//double dt_l= 1e-5,dt_u = 1e-3;

	double dx_l=4e-3,dx_u = 4e-2;
	double dt_l= 1e-5,dt_u = 1e-2;

	double ddx = (dx_u-dx_l)/(10.0);	
	double ddt = (dt_u-dt_l)/(10.0);

	FILE *fp = fopen("imex_linq_cmplx.txt","w");

	//for(dt=dt_l;dt<=dt_u;dt+=ddt)
	{
		dt = 1e-4;
		dx = 2e-2;
		
	//	for(dx = dx_l;dx<=dx_u;dx+=ddx)
		{
			
			en_loss = run(dt,dx,&abs_err,&stb_avg,0,1,1);

			printf("%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err,stb_avg);
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%.10lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err,stb_avg);
		}

	}


}
