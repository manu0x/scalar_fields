int calculate_vel_from_psi(int *n,double *dx,fdm_psi_mpi psi,double v[][3],double &vmax,double a)
{	
	int i,j,k,ci,tN,ind[3],cnt,ret =1;
	double psi_ret[2];
	double psi_amp,dtN,ksqr,theta_x[3],m[3];
	int ind_l1[3]{0,0,0},ind_l2[3]{0,0,0},ind_r1[3]{0,0,0},ind_r2[3]{0,0,0};

	 double ***theta;
	

	tN = n[0]*n[1]*n[2];
	dtN = (double) tN;

	printf("velcal ...  %d %d %d\n",n[0],n[1],n[2]);
	
	theta = new double** [n[0]+4] ;
	double *theta_pool;
	vmax = 0.0;   
	
	   for(i=0;(i<n[0]+4);++i)
	   {
		theta[i] = new double* [n[1]];
		theta_pool =  new double [n[1]*n[2]];
		
		for( j=0;j<n[1];++j)
	     	{
		  theta[i][j] = theta_pool;
		  //f_t[i][j] = f_t_pool;

		  theta_pool+=n[2];
		  
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
			  //theta[i][j][k] = 0.0;
					
			
			

		   }

		}

	}


//////////////////////////  MPI_TO_DO : Write SEND_RECEIVE here...................///////////////////////


		int mpi_check = 0;	
		MPI_Request send_req,recv_req;
	
	
		int left_rank,right_rank,my_cart_rank;
		MPI_Status status;
		mpi_check=MPI_Cart_shift(cart_comm,0,+1,&my_cart_rank,&right_rank);
		mpi_check=MPI_Cart_shift(cart_comm,0,-1,&my_cart_rank,&left_rank);

		int sendtag = 1,recvtag = 1;


		mpi_check = MPI_Isend(&theta[2][0][0],1, c_x_plain, left_rank,sendtag, cart_comm,&send_req);
		mpi_check = MPI_Irecv(&theta[n[0]+2][0][0],1, c_x_plain,right_rank , recvtag,cart_comm, &recv_req);

		mpi_check = MPI_Isend(&theta[3][0][0],1, c_x_plain, left_rank,sendtag, cart_comm,&send_req);
		mpi_check = MPI_Irecv(&theta[n[0]+3][0][0],1, c_x_plain,right_rank , recvtag,cart_comm, &recv_req);

		

	


		mpi_check = MPI_Isend(&theta[n[0]][0][0],1, c_x_plain, right_rank,sendtag, cart_comm,&send_req);
		mpi_check = MPI_Irecv(&theta[0][0][0], 1, c_x_plain,left_rank , recvtag,cart_comm, &recv_req);

		mpi_check = MPI_Isend(&theta[n[0]+1][0][0],1, c_x_plain, right_rank,sendtag, cart_comm,&send_req);
		mpi_check = MPI_Irecv(&theta[1][0][0], 1, c_x_plain,left_rank , recvtag,cart_comm, &recv_req);

		MPI_Barrier(cart_comm);
		






/////////////////////////////////////////////////////////////////////////////////////////////////////////////


	
	
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
	     			if(ci==0)
	     			{ ind_l1[ci] = (ind[ci]-1);
	     			  ind_l2[ci] = (ind[ci]-2);

	     			  ind_r1[ci] = (ind[ci]+1);
	     			  ind_r2[ci] = (ind[ci]+2);
				}

				else
				{ ind_l1[ci] = (n[ci]+ind[ci]-1)%n[ci];
	     			  ind_l2[ci] = (n[ci]+ind[ci]-2)%n[ci];

	     			  ind_r1[ci] = (ind[ci]+1)%n[ci];
	     			  ind_r2[ci] = (ind[ci]+2)%n[ci];
				}
				

	    
	     
	       			m[0] = (-8.0*theta[ind_l1[0]][ind_l1[1]][ind_l1[2]]+8.0*theta[ind_r1[0]][ind_r1[1]][ind_r1[2]]); 
	       			m[1] = (theta[ind_l2[0]][ind_l2[1]][ind_l2[2]]);
	       			m[2] = m[0] + m[1] ;
		
				if(m[2]>1e8)
				printf("RRR %d %d %d %e %e\n",ind_l2[0],ind_l2[1],ind_l2[2],theta[ind_l2[0]][ind_l2[1]][ind_l2[2]],theta[ind_r2[0]][ind_r2[1]][ind_r2[2]]);

		
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


	printf("Hey vmax  %lf %lf %lf %lf\n",vmax,dx[0],dx[1],dx[2]);
	return(ret);	  

}
