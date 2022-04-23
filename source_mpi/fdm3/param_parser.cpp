int parser(char * param_file_name, param_fdm &p)
{
	FILE *fp_param = fopen(param_file_name,"r");
	
	//p.load_default();

	int max_len = 80;
	char read_c_line[max_len];
	string read_str, read_key, read_val;

	double num_val,h=-1.0;

	if(fp_param!=NULL)
	{
	  while(fgets(read_c_line,max_len,fp_param)!=NULL)
		{
		  read_str = read_c_line;

		 
		  read_key = read_str.substr(0,read_str.find("="));
		  remove_if(read_key.begin(),read_key.end(),::isspace);
		  read_key.erase(remove_if(read_key.begin(),read_key.end(),::isspace),read_key.end());
		 

		if(read_key.length()>0)
		{  read_val = read_str.substr(read_str.find("=")+1);
		   

		   if(read_key!="fini_dc")
			num_val = stod(read_val);
		
		   if(read_key=="h")
			p.loc_h = num_val;
		   else
		   if(read_key=="omega_dm_0")
			p.omega_dm_0 = num_val;
		   else
		   if(read_key=="alpha")
			p.loc_alpha_m22 = num_val;
		   else
		   if(read_key=="z_ini")
			p.z_ini = num_val;
		
		   
		
		   else
			if(read_key=="box_n")
			p.box_n = num_val;

		   else
			if(read_key=="box_length")
			p.box_length = num_val;
		   else
			if(read_key=="fini_dc")
			p.fini_dc = read_val.c_str();
		
			

		 } 

		  
			


		}

		
		return(1);


	}
	
	else	
	{

		printf("\n\nERROR from parser(param_parser.cpp)...No parameter file found\n\n");
		return(0);



	}





}
