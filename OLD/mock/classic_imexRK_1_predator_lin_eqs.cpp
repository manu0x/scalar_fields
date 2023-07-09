using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

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


double b(double x,double Pv)
{
	double eps = 0.005;	

	if(x<=0.5)
	return(eps/(eps+Pv));
	else
	return(100.0*eps/(eps+Pv));


}

double f(double t,double x)
{
	double r;
	
	if(t==0.0)
	{
		
		r = ((double)rand())/((double) RAND_MAX);
		r = 0.8 + r*(1.2-0.8);
		//printf("%lf %lf\n",x,r);
	}

	

	else
	r=0.0;
	return(r);


}


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

double ex_vel(double t,double x,double Pv)
{
	double fval,bval,res;
	double rd = 1.0;
	
	fval = f(t,x);
	bval = b(x,Pv);

	return(fval+bval-rd*Pv);



}

void cr_invert_mat(double *Mp,int N, double dt, double dx, double diff,double gamma)
{	///////////////////// This assumed Diagonal with same element gamma
	int i,j,k;
	double mu = dt/(dx*dx);
	size_t tn  = N;
	gsl_matrix *gM = gsl_matrix_alloc(tn, tn);
	gsl_matrix *gMinv = gsl_matrix_alloc(tn, tn);
	gsl_matrix_set_zero(gM);
	int sig, out;
	double M[N][N]={0.0};
	double res;
	FILE *fptest = fopen("test_inv1.txt","w");

///////////////Create Matrix////////////////////////////////////////////////////////
	
	for(i=1;i<N-1;++i)
	{
	  for(j=0;j<N;++j)
	  {
		//M[i][j] = 0.0;
		if(i==j)
		{
		  M[i][j] = 1.0 + 2.0*diff*gamma*mu;	
		  M[i][j-1] =  - diff*gamma*mu;
		  M[i][j+1] =  - diff*gamma*mu;

		  gsl_matrix_set(gM, i, j,   M[i][j]);
		  gsl_matrix_set(gM, i, j-1,   M[i][j-1]);
		  gsl_matrix_set(gM, i, j+1,   M[i][j+1]);
		}
		

	  }



	}

	for(j=0;j<N;++j)
	{ M[0][j] = 0.0;
	  M[N-1][j] = 0.0;

	  
	  M[0][0] = 1.0 + 2.0*diff*gamma*mu;
	  M[N-1][N-1] = 1.0 + 2.0*diff*gamma*mu;	

	 gsl_matrix_set(gM, 0, j,   M[0][j]);
	 gsl_matrix_set(gM, N-1, j,   M[N-1][j]);
	}

	M[0][N-1] =  - diff*gamma*mu;
	M[0][1] =  - diff*gamma*mu;
	M[N-1][0] =  - diff*gamma*mu;
	M[N-1][N-2] =  - diff*gamma*mu;

	gsl_matrix_set( gM, 0, N-1,   M[0][N-1]);
	gsl_matrix_set( gM, 0, 1,   M[0][1]);
	gsl_matrix_set( gM, N-1, 0,   M[N-1][0]);
	gsl_matrix_set( gM, N-1, N-2,   M[N-1][N-2]);


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

		fprintf(fptest,"%.10lf ",M[i][k]);
			 
	  }

		fprintf(fptest,"\n");

	}




fclose(fptest);


}


void test_invert_mat(double *Mp,int N, double dt, double dx, double diff,double gamma)
{	///////////////////// This assumed Diagonal with same element gamma
	int i,j,k;
	double mu = dt/(dx*dx);
	double res;
	size_t tn  = N;
	
	int sig, out;
	double Mn[N][N];//={0.0};

	FILE *fptest = fopen("test_inv.txt","w");
///////////////Create Matrix////////////////////////////////////////////////////////
	printf("M[0] %lf %lf %lf\n",Mn[0][0],Mn[0][1],Mn[0][2]);

	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {
		Mn[i][j] = 0.0;
		

	  }



	}
	printf("M[0] %lf %lf %lf\n",Mn[0][0],Mn[0][1],Mn[0][2]);
	for(i=1;i<N-1;++i)
	{
	  for(j=0;j<N;++j)
	  {
		
		if(i==j)
		{
		  Mn[i][j] = 1.0 + 2.0*diff*gamma*mu;	
		  Mn[i][j-1] = 1.0 - diff*gamma*mu;
		  Mn[i][j+1] = 1.0 - diff*gamma*mu;
		}
		

	  }



	}
	
	
	  
	  Mn[0][0] = 1.0 + 2.0*diff*gamma*mu;
	  Mn[N-1][N-1] = 1.0 + 2.0*diff*gamma*mu;

	
	

	Mn[0][N-1] = 1.0 - diff*gamma*mu;
	Mn[0][1] = 1.0 - diff*gamma*mu;
	Mn[N-1][0] = 1.0 - diff*gamma*mu;
	Mn[N-1][N-2] = 1.0 - diff*gamma*mu;


////////////////////////////////////////////////////////////////////////////////////////////////////	
	for(i=0;i<N;++i)
	{
	  for(j=0;j<N;++j)
	  {	res=0.0;
		for(k=0;k<N;++k)
			res+= Mp[i*N+k]*Mn[k][j];

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



fclose(fptest);

}






void initialise(double *P,double *x,double *k,double dx,int N)
{

	int i,j; double di,dk;
	dk = 2.0*pie/(((double)(N))*dx);
	for(i=0;i<N;++i)
	{
	   di = (double)i;	
	   P[i] = 0.0;

	    if(i<=(N/2))
		*(k+i) = di*dk;
	    else
		*(k+i) = (di-((double) N))*dk; 

		x[i] = 0.0+di*dx;

	//	printf("k %lf  %lf  dx %lf\n",dk,*(k+i),dx );

	}


}






int main()
{

	int N,t_steps;
	double box_len,dx,t_end,t_start,dt,diff;


///////////////////////////////File for data/////////////////////////////////////////////////////

	FILE *fp = fopen("data_pred.txt","w");
	FILE *fp2 = fopen("data2_pred.txt","w");
	FILE *fptemp = fopen("temp.txt","w");
	
	FILE *fpmass = fopen("mass_pred.txt","w");

/////////////////////////////// Parameter setting ////////////////////////////////////////////////
	m  = 1.0;
	n  = 1.0;
	T = 2.0*pie*m;
	diff = 0.02;
	N = 10000;



/////////////////////////////////////////RK things/////////////////////////////////////////


	double im_K_P[ex_s][N],ex_K_P[im_s][N];


/////////////////////////////// Box & res. setting ///////////////////////////////////////////////

	box_len = 1.0;
	dx = box_len/(double(N));

////////////////////////////// Time & dt settings ////////////////////////////////////////////////
	
	t_start = 0.0;
	t_end = 10.0;
	t_steps = 40;
	dt  = (t_end-t_start)/((double)t_steps);
	printf("dt %lf\n",dt);
	

////////////////////////////// P variables  /////////////////////////////////////////////////////

	double P[N],lap_val[2*N],lambda;

	
	double x[N],k_grid[N];
	
	double Pb[N],Pk[N],im_K[im_s][N],ex_K[ex_s][N],dP_xx[N];

	initialise(P,x,k_grid,dx,N);

	int i,j;
	
	 

///// Get inverted matrix ///////////////////////////



  double *Mp = new double [N*N];
  double mat[N][N];


  cr_invert_mat(Mp, N,  dt,  dx, diff,im_a[0][0]);

  test_invert_mat(Mp,N, dt, dx, diff,im_a[0][0]);

  for(i=0;i<N;++i)
  { for(j=0;j<N;++j)
    {
	mat[i][j] = Mp[i*N+j];	
    }

	Pb[i] = P[i];
  }

  delete [] Mp;




///////////////////////  Evolution ///////////////////////////////////////////////////////////////
	srand(time(0));

	int s_cntr,tcntr,printcntr,fail=0;
	int l1,r1;

	printcntr = (int)(((double) t_steps)/100.0);	printf("%d\n",printcntr);
	if(t_steps<=100)
		printcntr = 1.0;
	double t,vel_val,vl1,vr1,vc;



	

		
	

	for(t=t_start,tcntr=0;(t<=t_end)&&(!fail)&&(1);t+=dt,++tcntr)
	{	
		
	
	 	for(i=0;i<N;++i)
		{	
			
			if((tcntr%printcntr)==0)  
			{
			  
			  
			  fprintf(fp,"%lf\t%lf\n",x[i],P[i]);
				if(i==5)
			  fprintf(fptemp,"%lf\t%lf\n",t,P[i]);

			   
			}
			
			
			Pk[i] = 0.0;

			for(j=0;j<N;++j)	
			{	
				Pk[i]+=mat[i][j]*Pb[j];

			


			}
			

			
		}
		
		


		for(s_cntr=1;s_cntr<imex_s;++s_cntr)
		{

			for(j=0;j<s_cntr;++j)
			{
			   for(i=0;i<N;++i)
			   {
				
				if(j==0)
				{

					if(s_cntr==1)
		    			vel_val = ex_vel(t+ex_c[s_cntr-1]*dt,x[i],P[i]);
					else
					vel_val = ex_vel(t+ex_c[s_cntr-1]*dt,x[i],Pk[i]);

		    			ex_K_P[s_cntr-1][i] = vel_val; 

					l1 = ((N)+ (i-1))%(N);
					r1 = (i+1)%(N);
	
					vl1 = Pk[l1];
					vr1 = Pk[r1];	
					vc = Pk[i];
	
					im_K_P[s_cntr-1][i]  =  diff*(vr1 +vl1 - 2.0*vc )/(dx*dx);
					//printf("%d %d %lf\n",s_cntr,i,im_K_P[s_cntr-1][i] );

					Pb[i] = P[i] + dt*ex_a[s_cntr][j]*ex_K_P[j][i]+ dt*im_a[s_cntr][j]*im_K_P[j][i];
					

				}
			
				
				else
				Pb[i]+=  dt*ex_a[s_cntr][j]*ex_K_P[j][i]+ dt*im_a[s_cntr][j]*im_K_P[j][i];
				 
				
				
	
	
			   }

	
			}
	

			  for(i=0;i<N;++i)
			  {	
				Pk[i] = 0.0;

				for(j=0;j<N;++j)	
				{
					Pk[i]+=mat[i][j]*Pb[j];


				}
							

			     if(Pk[i]<0.0)
				{

						printf("+ity broken interim %d %d %d %lf\n",s_cntr,j,i,Pk[i]);
						fail = 1;				
						break;

				}
		
				
			
			   }
				
			


		}//RK stages concluded

		if((tcntr%printcntr)==0)  
		{
			  
			   printf("%lf\n",t/t_end);
		}
			
		

		for(i=0;i<N;++i)
		{	
	
			vel_val = ex_vel(t+ex_c[imex_s-1]*dt,x[i],Pk[i]);

		    	ex_K_P[imex_s-1][i] = vel_val; 

			l1 = ((N)+ (i-1))%(N);
			r1 = (i+1)%(N);
			vl1 = Pk[l1];
			vr1 = Pk[r1];
			vc = Pk[i];
	
			im_K_P[imex_s-1][i]  =  diff*(vr1 +vl1 - 2.0*vc )/(dx*dx);


			

	
							
			   for(j=0;j<imex_s;++j)
			   {	P[i]+=  dt*ex_b[j]*ex_K_P[j][i]+ dt*im_b[j]*im_K_P[j][i];
				

			
				//if(i==51)
				//printf("!!!+tivity broken!!!j  %d\ti %d\t%lf\t%lf\n",j,i,P[i],im_K_P[j][i]);
			    }


				//if(isnan(P[i]))				
			if(isnan(P[i])||(P[i]<0.0))
			 {

				if(isnan(P[i])&&(prt))
				printf("!!!GONE NAN!!!!  %d\tt %lf\t%lf\n",tcntr,t,P[i]);
				if(P[i]<0.0&&(prt))
				printf("!!!+tivity broken!!!j  %d\ti %d\t%lf\t%lf\n",j,i,P[i],im_K_P[j][i]);
					
				fail = 1;

				break;
			}
		


			if(fail)
			break;

			Pb[i] = P[i];


			


			

			if((tcntr%printcntr)==0)  
			fprintf(fp2,"%lf\t%lf\t%lf\n",x[i],P[i],Pb[i]);
		

			if(((100.0*fabs(avg_amp-amp_ini)/amp_ini)>=1e3)||(fail)||((*abs_err)>=1e3))
			{

			
				*abs_err = 1e3;
					
				return(1e3);
	
		
			}
			

		   }



		if((tcntr%printcntr)==0)  
		{ fprintf(fp,"\n\n\n");
		  fprintf(fp2,"\n\n\n");
		  fprintf(fpmass,"%lf\n",t/t_end);
		  printf("dt %lf\tt %lf\n\n",dt,t);
		}
		


		  

	}///// ENd f Time Evolution /////////////////////////
	
	





	return(0);

}
