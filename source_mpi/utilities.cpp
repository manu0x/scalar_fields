int calculate_vel_from_psi(int *n,double *dx,fdm_psi_mpi psi,double v[][3],double &vmax,double a)
{	
	int i,j,k,ci,tN,ind[3],cnt,ret =1;
	double psi_ret[2];
	double psi_amp,dtN,ksqr,theta_x[3],m[3];
	int ind_l1[3]{0,0,0},ind_l2[3]{0,0,0},ind_r1[3]{0,0,0},ind_r2[3]{0,0,0};

	 double ***theta;
	

	tN = n[0]*n[1]*n[2];
	dtN = (double) tN;
	
	theta = new double** [n[0]+4] ;
	vmax = 0.0;   
	
	   for(i=0;(i<n[0]+4);++i)
	   {
		theta[i] = new double* [n[1]];
		
		for( j=0;j<n[1];++j)
	     	{
		  theta[i][j] = new  double[n[2]] ;
		  
		 }

	    }



	for(i=0;i<n[0];++i)
	{
		for(j=0;j<n[1];++j)
		{
		  for(k=0;k<n[2];++k)
		  {
			
			ind[0] = i; ind[1] = j; ind[2] = k;
			psi.get_psi(ind,psi_ret);
			psi_amp = sqrt(psi_ret[0]*psi_ret[0] + psi_ret[1]*psi_ret[1]);
			theta[i+2][j][k] = 0.5*(asin(psi_ret[1]/psi_amp)+ acos(psi_ret[0]/psi_amp));			
			


		   }

		}

	}


//////////////////////////  MPI_TO_DO : Write SEND_RECEIVE here...................///////////////////////

	
	
	for(i=0;i<n[0];++i)
	{	ind_l1[0] = i+2;
		ind_l2[0] = i+2;
		ind_r1[0] = i+2;
		ind_r2[0] = i+2;
		ind[0] = i+2; 

		for(j=0;j<n[1];++j)
		{ ind_l1[1] = j;
		  ind_l2[1] = j;
		  ind_r1[1] = j;
		  ind_r2[1] = j;
		  ind[1] = j;
		  for(k=0;k<n[2];++k)
		  {
			ind_l1[2] = k;
			ind_l2[2] = k;
			ind_r1[2] = k;
			ind_r2[2] = k;
			ind[2] = k;
			
			cnt = (n[2]*n[1])*i + n[2]*j + k;
						
			
			for(ci=0;ci<3;++ci)
	   		{
	     
	     			if(ci>0)
	     			{ind_l1[ci] = (n[ci]+ind[ci]-1)%n[ci];
	     			 ind_l2[ci] = (n[ci]+ind[ci]-2)%n[ci];

	     			 ind_r1[ci] = (ind[ci]+1)%n[ci];
	     			 ind_r2[ci] = (ind[ci]+2)%n[ci];
				}
				else
				{ind_l1[ci] = (ind[ci]-1);
	     			 ind_l2[ci] = (ind[ci]-2);

	     			 ind_r1[ci] = (ind[ci]+1);
	     			 ind_r2[ci] = (ind[ci]+2);
				}

	    
	     
	       			m[0] = (-8.0*theta[ind_l1[0]][ind_l1[1]][ind_l1[2]]+8.0*theta[ind_r1[0]][ind_r1[1]][ind_r1[2]]); 
	       			m[1] = (theta[ind_l2[0]][ind_l2[1]][ind_l2[2]]-theta[ind_r2[0]][ind_r2[1]][ind_r2[2]]);
	       			m[2] = m[0] + m[1] ;
		
		
	       			theta_x[ci] = m[2]/(dx[ci]);
	      

	    
				ind_l1[ci] = ind[ci];
	    			ind_l2[ci] = ind[ci];

	     			ind_r1[ci] = ind[ci];
	     			ind_r2[ci] = ind[ci];

				v[cnt][ci] = hbar_by_m*theta_x[ci];
				if(vmax<fabs(v[cnt][ci]))
					vmax = fabs(v[cnt][ci]);
				if(isnan(v[cnt][ci]))
					ret = 0;

	   		}			
			


		   }

		}

	}

/////////////////// MPI_TO_DO : Write a collective MAX call....................////////////////////////


	//printf("Hey vmax  %lf\n",vmax);
	return(ret);	  

}
