using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <omp.h>

#include <gsl/gsl_linalg.h>

#define pie M_PI 



//////////GLobal constants/////////////

double m,n,T;
//double theta;	// How much implicit;

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


double V(double psi_amp)
{
	//return(0.0);
	return( 1.0 - 2.0*pie*pie*n*n/(m*m));


}


void ex_vel(double v[2],double psi[2],double Vval)
{

	v[0] = m*Vval*psi[1];
	v[1] = -m*Vval*psi[0];;



}

void cr_invert_mat(double *Mp,double *MpP,int N, double dt, double dx, double diff,double Vval,double theta)
{	///////////////////// This assumed Diagonal with same element gamma
	int i,j,k;
	double mu = dt/(dx*dx);
	//double omega = mu*gamma/(2.0*m);
	size_t tn  = N;
	gsl_matrix *gM = gsl_matrix_alloc(tn, tn);
	gsl_matrix *gMinv = gsl_matrix_alloc(tn, tn);
	gsl_matrix_set_zero(gM);
	int sig, out;
	double M[N][N]={0.0};
	double P[N][N]={0.0};
	double P2[N][N]={0.0};
	double res;
	FILE *fptest = fopen("test_inv1.txt","w");

///////////////create P^2  ///////////////////////////////////////////////////////

	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {
	     M[i][j] = 0.0;
	     P[i][j] = 0.0;
	     if(i==j)
		M[i][j] = dt*theta*(m*Vval + 2.0/(2.0*m*dx*dx));
		
	   }
	}


	for(i=1;i<N-1;++i)
	{
	  for(j=0;j<N;++j)
	  {
		//M[i][j] = 0.0;
		if(i==j)
		{
		  
		  M[i][j-1] = -dt*theta/(2.0*m*dx*dx);
		  M[i][j+1] = -dt*theta/(2.0*m*dx*dx);

		
		}
		

	  }



	}

	M[0][N-1] = -theta*dt/(2.0*m*dx*dx);
	M[N-1][0] = -theta*dt/(2.0*m*dx*dx);
	M[0][1] = -theta*dt/(2.0*m*dx*dx);
	M[N-1][N-2] = -theta*dt/(4.0*m*dx*dx);

	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {

		P2[i][j]=0.0;
		P[i][j] = M[i][j];
	  	for(k=0;k<N;++k)
		{
		  P2[i][j]+=M[i][k]*M[k][j];

		}

		
	   }	
	  

	}

	



///////////////Create Matrix////////////////////////////////////////////////////////
	
	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {
		//M[i][j] = 0.0;
		if(i==j)
		 M[i][j] = 1.0+P2[i][j];	
		else
		  M[i][j] = P2[i][j];

		  gsl_matrix_set(gM, i, j,   M[i][j]);
		  
		fprintf(fptest,"%lf ",M[i][j]);
		

	  }

	 fprintf(fptest,"\n");

	}

	fprintf(fptest,"\n\n");


////////////////////// LU decomp and inversion ////////////////////////////////////////////////////

	gsl_permutation *p = gsl_permutation_alloc(tn);
	out = gsl_linalg_LU_decomp(gM, p, &sig);
	out = gsl_linalg_LU_invert(gM, p,gMinv);

////////////////////////////////////////////////////////////////////////////////////////////////////	
	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {
		Mp[i*N+j] = gsl_matrix_get(gMinv,i,j);
		MpP[i*N+j] = 0.0;
		for(k=0;k<N;++k)
		{
			MpP[i*N+j]+=(gsl_matrix_get(gMinv,i,k)*P[k][j]);
	

		}
	  }

	}


	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {	res=0.0;
		for(k=0;k<N;++k)
			res+= Mp[i*N+k]*M[k][j];

		//if(i==j)
			

		fprintf(fptest,"%lf ",res);
			 
	  }

		fprintf(fptest,"\n");

	}

	fprintf(fptest,"\n\n\n");
	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {	

		fprintf(fptest,"%.10lf ",Mp[i*N+j]);
			 
	  }

		fprintf(fptest,"\n");

	}

	fprintf(fptest,"\n\n\n");
	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {	

		fprintf(fptest,"%.10lf ",MpP[i*N+j]);
			 
	  }

		fprintf(fptest,"\n");

	}




fclose(fptest);


}








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

double run(double dt,double dx,double *abs_err,int printfp,int prt,double theta=1.0)
{
	int N,t_steps;
	double box_len,t_end,t_start,diff;

///////////////////////////////File for data/////////////////////////////////////////////////////

	FILE *fp = fopen("data_trp.txt","w");
	FILE *fp2 = fopen("data2_trp.txt","w");
	FILE *fptemp = fopen("temp_trp.txt","w");
	
	FILE *fpmass = fopen("mass_trp.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	m  = 1.0;
	n  = 1.0;
	T = 2.0*pie*m;
	

	theta = 1.0;





/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 2.0;
//	dx = box_len/(double(N));
	N = ((int)(box_len/dx)) + 1;

////////////////////////////// Time & dt settings ////////////////////////////////////////////////
	
	t_start = 0.0;
	t_end = 2.0*T;
	t_steps = (int)((t_end-t_start)/dt);
	//dt  = (t_end-t_start)/((double)t_steps);
	printf("dt %lf N %d\n",dt,N);	

/////////////////////////////////////////RK things/////////////////////////////////////////


	double im_K_P[ex_s][N][2],ex_K_P[im_s][N][2];


////////////////////////////// P variables  /////////////////////////////////////////////////////

	double Psi[N][2],lap_val[2*N],lambda;

	
	double x[N],k_grid[N];
	
	double Psib[N][2],Psik[N][2],im_K[im_s][N][2],ex_K[ex_s][N][2],dP_xx[N];
	double vel_val[2],c_psi[2],c_psi_amp,Vval;

	initialise(&Psi[0][0],x,k_grid,dx,N);

	int i,j;
	
	 

///// Get inverted matrix ///////////////////////////



  double *Mp = new double [N*N];
  double *MpP = new double [N*N];
  double mat[N][N];
  double matP[N][N];

  Vval = V(0.0);
  cr_invert_mat(Mp,MpP, N,  dt,  dx, diff,Vval,theta);

  

  for(i=0;i<N;++i)
  { for(j=0;j<N;++j)
    {
	mat[i][j] = Mp[i*N+j];
	matP[i][j] = MpP[i*N+j];	
    }

	Psib[i][0] = Psi[i][0];
	Psib[i][1] = Psi[i][1];
  }

  delete [] Mp;
  delete [] MpP;




///////////////////////  Evolution ///////////////////////////////////////////////////////////////
	srand(time(0));

	int s_cntr,tcntr,printcntr,fail=0;
	int l1,r1;

	printcntr = (int)(((double) t_steps)/100.0);	printf("%d\n",printcntr);
	if(t_steps<=100)
		printcntr = 1.0;
	double t,vl1,vr1,vc,amp,avg_amp,sol[2],dbi;
	double drc = 2.0*pie*n*2.0*pie*n;
	double fdt,amp_ini;

	double omega = dt*im_a[0][0]/(dx*dx*2.0*m);

	

		
	

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(1);t+=dt,++tcntr)
	{	
	
		if((tcntr%printcntr)==0) 
		avg_amp = 0.0;


		#pragma omp parallel for private(l1,r1,vl1,vr1,vc)
		for(i=0;i<N;++i)
		{	
			
			l1 = i-1;
			r1 = (i+1);
			if(i==0)
			 l1 = N-2;
			if(i==(N-1))	
			  r1 = 1;
	
			vl1 = Psi[l1][1];
			vr1 = Psi[r1][1];	
			vc = Psi[i][1];

			Psib[i][0] = Psi[i][0]+ (1.0-theta)*dt*m*Vval*Psi[i][1] - (1.0-theta)*dt*(vr1 +vl1 - 2.0*vc )/(2.0*m*dx*dx);

			
	
			vl1 = Psi[l1][0];
			vr1 = Psi[r1][0];	
			vc = Psi[i][0];

			Psib[i][1] = Psi[i][1] - (1.0-theta)*dt*m*Vval*Psi[i][0] + (1.0-theta)*dt*(vr1 +vl1 - 2.0*vc )/(2.0*m*dx*dx); 

		}
		(*abs_err) = 0.0;
	 	for(i=0;i<N;++i)
		{	
			
			dbi = (double)(i);
			
			
			
			
			Psi[i][0] = 0.0;
			Psi[i][1] = 0.0;

			for(j=0;j<N;++j)	
			{	
				Psi[i][0]+=( mat[i][j]*Psib[j][0] + matP[i][j]*Psib[j][1] );

			        Psi[i][1]+=(  -matP[i][j]*Psib[j][0] + mat[i][j]*Psib[j][1] );


			}

			if(isnan((Psi[i][0])+Psi[i][1]))
			{
				fail =1;
				if(prt)
				printf("FAILED %d tcntr %d %lf %lf\n",i,tcntr,Psi[i][0],Psi[i][0]);
			
				break;

			}

			amp  = sqrt(Psi[i][0]*Psi[i][0] + Psi[i][1]*Psi[i][1]);				
			avg_amp+=amp;
			  
			sol[0] = cos(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi);
			sol[1] = -sin(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi);


			if((tcntr%printcntr)==0)  
			{
			  
			  fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,Psi[i][0],Psi[i][1],sol[0],sol[1],amp);
			  if(i==20)
			   fprintf(fptemp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,Psi[i][0],Psi[i][1],sol[0],sol[1],amp);
			}


			
			(*abs_err)+=(fabs(sol[0]-Psi[i][0])+fabs(sol[1]-Psi[i][1]));
			
		}
	
		*abs_err = (*abs_err)/((double)N);

		
		
		
		if(t==t_start)
		{amp_ini = avg_amp;
		 
			if((tcntr%printcntr)==0)  
			{ 
			  sol[0] = cos(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi);
			  sol[1] = -sin(2.0*pie*(t+dt)/T)*sin(2.0*pie*n*dx*dbi);	
			  if(printfp)
			  fprintf(fp2,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,Psi[i][0],Psi[i][1],sol[0],sol[1],amp);
			}
			

			

		  }

		if(((100.0*fabs(avg_amp-amp_ini)/amp_ini)>=1e3)||(fail)||((*abs_err)>=1e3))
		{

			
				*abs_err = 1e3;
					
				return(1e3);
	
		
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
	

	avg_amp = 0.0;
	for(i=0;i<N;++i)
		{	
			
			dbi = (double)(i);
			
			amp  = sqrt(Psi[i][0]*Psi[i][0] + Psi[i][1]*Psi[i][1]);				
			avg_amp+=amp;
			  
			  sol[0] = cos(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			  sol[1] = -sin(2.0*pie*t/T)*sin(2.0*pie*n*dx*dbi);
			 if(printfp)
			  fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",dx*dbi,Psi[i][0],Psi[i][1],sol[0],sol[1],avg_amp);
			 if((i==20)&&(printfp))
			   fprintf(fptemp,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,Psi[i][0],Psi[i][1],sol[0],sol[1],amp);
			
		}
	
	if(prt)
	printf("Error in conserv. %lf\n mu is %lf", 100.0*fabs(avg_amp-amp_ini)/amp_ini,dt/(dx*dx));



	return(100.0*fabs(avg_amp-amp_ini)/amp_ini);


}




int main()
{
	double dt = 3e-4;
	double dx = 4e-3;
	double abs_err,en_loss;

	double dx_l=4e-3,dx_u = 4e-2;
	double dt_l= 1e-5,dt_u = 1e-2;

	double ddx = (dx_u-dx_l)/(20.0);	
	double ddt = (dt_u-dt_l)/(20.0);

	FILE *fp = fopen("trapezoid_fd.txt","w");

	for(dt=dt_l;dt<=dt_u;dt+=ddt)
	{

		
		for(dx = dx_l;dx<=dx_u;dx+=ddx)
		{
			en_loss = run(dt,dx,&abs_err,0,0);

			printf("%lf\t%lf\t%lf\t%lf\t%lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err);
			fprintf(fp,"%lf\t%lf\t%lf\t%lf\t%lf\n",dx,dt,dt/(dx*dx),en_loss,abs_err);
		}

	}


}
