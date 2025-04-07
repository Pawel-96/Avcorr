#include "Parsfnames.h"


string Fin(string model) //input
{
    return "Data/"+model+EXT;
}




string Mout_onerad(string model, int nnsize) //moment output for one radius
{
    string fin=Fin(model);
    string mout=Replace_string(fin,EXT,"_moments_nnR_"+conv(nnsize)+".txt");
    return Replace_string(mout,"Data/","Results/");
}




string Merrout(string model) //output for one model and all moments
{
    string f1=Replace_string(Fin(model),"Data/","Results/");
    return Replace_string(f1,EXT,"_merrout.txt");
}




void Writelog(string text, string option) //log
{
	ofstream lfile;
    lfile.open(logfile.c_str(),fstream::app);
	
	if(option=="write")
	{
		time_t now = time(0); // get current date and time  
		tm* ltm = localtime(&now);  
		string h=conv(ltm->tm_hour)+":"+conv(ltm->tm_min)+":"+conv(ltm->tm_sec);
		
		lfile<<h<<" "<<text<<endl;
	}
	
    if(option=="intro")
    {
		auto T= chrono::system_clock::now();
		time_t T_time = std::chrono::system_clock::to_time_t(T);
		string t=ctime(&T_time);
        
		lfile<<"---------------START: "<<VERSION<<" --------------- "<<t<<endl;
        lfile<<"[Parameters]:"<<endl;
		
		lfile<<"kappa: "<<kappa<<endl;
        lfile<<"Cmin: "<<Cmin<<endl;
        lfile<<"Cmax: "<<Cmax<<endl;
        lfile<<"nreals: "<<nreals<<endl;
		lfile<<"Random_provided: "<<Random_provided<<endl;
		if(Random_provided==1){lfile<<"Random_file: "<<Random_file<<endl;}
        lfile<<"Datafiles ["<<nmodels<<"]:"<<endl;
        for(int i=0;i<nmodels;++i)
        {
            lfile<<"            "<<Model[i]<<endl;
        }
		lfile<<"Recalc: "<<Recalc<<endl;
		lfile<<"Cleaning: "<<Clean<<endl;
		lfile<<"Poisson errors: "<<ErrP<<endl;
		lfile<<"Combine_reals: "<<combine_reals<<endl;
		if(combine_reals==1){lfile<<"Real_template: "<<real_template<<endl;}
		
		if(VERSION=="angular")
		{
			lfile<<"RA: ["<<ramin<<":"<<ramax<<"] deg"<<endl;
			lfile<<"DEC: ["<<decmin<<":"<<decmax<<"] deg"<<endl;
			lfile<<"Areaf: "<<setprecision(15)<<Areaf<<" deg2"<<endl;
		}
		
		if(VERSION=="angular" or VERSION=="BOX")
		{
			lfile<<"Scales: "<<nR<<" in range: ["<<Rmin<<":"<<Rmax<<"] deg"<<endl;
		}
		
		if(VERSION=="BOX" or VERSION=="BOX_ellipses")
		{
			lfile<<"Boxsize: "<<boxsize<<endl;
		}
		
		if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses")
		{
			lfile<<"LOS parallel axis scales: "<<naxa<<" in range: ["<<axamin<<":"<<axamax<<"] Mpc"<<endl;
			lfile<<"LOS perpendicular axis scales: "<<naxb<<" in range: ["<<axbmin<<":"<<axbmax<<"]"<<endl;
		}

    }
	
    if(option=="outro")
    {
		auto T= chrono::system_clock::now();
		time_t T_time = std::chrono::system_clock::to_time_t(T);
		string t=ctime(&T_time);
		
		lfile<<"----------------END---------------- "<<t<<" have a good day :)"<<endl;
    }
	
	lfile.close();
    return;
}






//counting overall progress
double Progress(int iimax)
{
	string fname;
	double res=0.,val;
	regex pattern("(.*)(.fnd)(.*)");
	
	for (const auto &entry : fs::directory_iterator("Results/"))  //counting files with *.fnd*
	{
		fname=entry.path().filename().string();
        if(regex_match(fname, pattern))
		{
			res+=1;
		}
    }
	
	val=Round_cut(1.*res/(1.*iimax),4);
	if(val<1e-4 and res>0){return Round_cut(1.*res/(1.*iimax),8);}
	else {return val;}
	return 0;
}







void Error(int &err, string msg) //raising error with message
{
	err=1;
	cerr<<msg<<endl;
	Writelog(msg);
	return;
}






//pointing errors in paramfile
int Error_param()
{
	int err=0;
	string msg="";
	
	
	if(err_common!=0) //one or more of common parameters are not specified
	{
		Error(err,"[Error]: check "+paramfile+", "+conv(err_common)+" common parameter(s) not specified. How dare you!");
		return err;
	}
	
	
	if(err_ver!=0) //one or more of common parameters are not specified
	{
		Error(err,"[Error]: check "+paramfile+", "+conv(err_ver)+" VERSION-dependent parameter(s) not specified.\
		Some may not be needed for current VERSION, but better keep them in "+paramfile);
		return err;
	}
	
	if(err_conv!=0) //numeric parameters are empty or strings
	{
		Error(err,"[Error]: check "+paramfile+", "+conv(err_conv)+" numeric parameter(s) not specified/are strings");
		return err;
	}
	
	
	
	if(VERSION!="angular" and VERSION!="BOX" and VERSION!="BOX_ellipses" and VERSION!="LC_ellipses")
	{
		Error(err,"[Error]: check "+paramfile+", VERSION should be angular/BOX/BOX_ellipses/LC_ellipses.");
		return err;
	}
	
	
	if(USE_HDF5!=0 and USE_HDF5!=1)
	{
		Error(err,"[Error]: check "+paramfile+", USE_HDF5 must be either 0 or 1.");
	}
	
	if(nreals<3)
	{
		Error(err,"[Error]: check "+paramfile+", nreals should be >=3 for reliable results");
	}
	
	if(kappa>.5 or kappa<=0)
	{
		Error(err,"[Error]: check "+paramfile+", kappa should be in range (0,0.5] for reliable results");
	}
	
	
	if(Cmin<=0)
	{
		Error(err,"[Error]: check "+paramfile+", Cmin should be positive");
	}
	
	if(Cmin>=Cmax)
	{
		Error(err,"[Error]: check "+paramfile+", Cmin should be < Cmax");
	}
	
	
	if(Cmax>=CMAX)
	{
		Error(err,"[Error]: check "+paramfile+", Cmax should be < "+conv(CMAX));
	}
	
	
	if(Random_provided==1)
	{
		vector<string> frandoms;
		Get_fnames(frandoms,"Randoms");
		int nfrandoms=frandoms.size();
		
		if((Random_file=="*" and nfrandoms<nmodels) or nfrandoms==0)
		{
			Error(err,"[Error]: check "+paramfile+" or Randoms/ directory: some randoms may be missing");
		}
		
		if(Random_file!="*" and Fexist(Random_file)==0)
		{
			Error(err,"[Error]: check "+paramfile+" or Randoms/ directory: Randoms/"+Random_file+" does not exist");
		}
		
		if(Random_file=="*")
		{
			string fr;
			int ncols;
			for(int i=0;i<nmodels;++i)
			{
				fr="Randoms/Random_"+Model[i]+EXT;
				if(Fexist(fr)==0)
				{
					Error(err,"[Error]: check "+paramfile+", file: "+fr+" does not exist.");
				}
				else
				{
					ncols=Fncols(fr);
					if((VERSION=="angular" and ncols!=2) or (VERSION!="angular" and ncols!=3))
					{
						Error(err,"[Error]: file: "+fr+" has wrong number of columns, ");
					}
				}
			}
		}
	}
	
	if(Random_provided!=0 and Random_provided!=1)
	{
		Error(err,"[Error]: check "+paramfile+", Random_provided should be either 0 or 1.");
	}
	
	
	if(Recalc!=0 and Recalc!=1)
	{
		Error(err,"[Error]: check "+paramfile+", Recalc should be either 0 or 1");
	}
	
	
	if(Clean!=0 and Clean!=1)
	{
		Error(err,"[Error]: check "+paramfile+", Clean should be either 0 or 1");
	}
	
	
	if(ErrP!=0 and ErrP!=1)
	{
		Error(err,"[Error]: check "+paramfile+", ErrPoisson should be either 0 or 1");
	}
	
	
	
	if(combine_reals!=0 and combine_reals!=1)
	{
		Error(err,"[Error]: check "+paramfile+", Combine_reals should be either 0 or 1");
	}
	
	
	
	
	
	if(USE_HDF5==0) //only for ASCII: columns beyond data file(s)
	{
		int ncols=Fncols(Fin(Model[0]));
		
		if(stoi(cols_pos[0])<0 or stoi(cols_pos[0])>=ncols
		or stoi(cols_pos[1])<0 or stoi(cols_pos[2])>=ncols
		or stoi(cols_pos[1])<0 or stoi(cols_pos[2])>=ncols)
		{
			err=1;
			if(ncols!=-1) //if ncols==-1, first file doesnt exist, another error shows up
			{
				Error(err,"[Error]: check "+paramfile+", cols_pos exceeding the file - should be in range [0,"+conv(ncols-1)+"].");
			}
		}	
	}
	
	
	if(nmodels==0) //Data/ directory empty
	{
		Error(err,"[Error]: check "+paramfile+", no Datafiles specified or Data/ empty.");
		return err;
	}
	
	
	string modelsymbol=Get_parameter(paramfile,"Datafiles",err)[0];
	
	if(modelsymbol!="*") //file names written directly in paramfile and some doesnt exist
	{
		for(int i=0;i<nmodels;++i)
		{
			if(Fexist(Fin(Model[i]))==0)
			{
				Error(err,"[Error]: check "+paramfile+", file: Data/"+Model[i]+EXT+" does not exist.");
			}
		}
	}
	
	
	
	
	if(combine_reals==1)
	{
		if(real_template.size()==0)
		{
			Error(err,"[Error]: check "+paramfile+", Combine_reals=1 is used, but Real_template is empty");
		}
		else
		{
			for(int i=0;i<nmodels;++i)
			{
				if(Model[i].find(real_template) == string::npos)
				{
					Error(err,"[Error]: check "+paramfile+", Real_template not found in datafile names");
				}
			}
		}
	
	}
	
	
	
	
	
	if(VERSION=="BOX" or VERSION=="BOX_ellipses")
	{
		if(boxsize<=0)
		{
			Error(err,"[Error]: check "+paramfile+", Boxsize should be >0.");
		}
	}
	
	
	if(VERSION=="angular" or VERSION=="BOX")
	{
		if(Rmin<=0 or Rmax<=0 or nR<=0)
		{
			Error(err,"[Error]: check "+paramfile+", Rmin, Rmax and nR should be >0.");
		}
		
		
		if(Rmin>=Rmax)
		{
			Error(err,"[Error]: check "+paramfile+", Rmax should be >= Rmin.");
		}
		
		if(Rmax>=boxsize)
		{
			Error(err,"[Error]: check "+paramfile+", Rmax should be < boxsize");
		}
		
		if(nR>NRMAX)
		{
			Error(err,"[Error]: check "+paramfile+", too big nR; should be <="+conv(NRMAX)+".");
		}
	}
	
	
	if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses")
	{
		if(axamin<=0 or axamax<=0 or axbmin<=0 or axbmax<=0 or naxa<=0 or naxb<=0)
		{
			Error(err,"[Error]: check "+paramfile+", axamin/max, axbmin/max, naxa/b should be >0.");
		}
		
		if(axamin>=axamax or axbmin>=axbmax)
		{
			Error(err,"[Error]: check "+paramfile+", axamax should be >= axamin, same with axb.");
		}
		
		if(axamax>=boxsize or axbmax>=boxsize)
		{
			Error(err,"[Error]: check "+paramfile+", axamax and axbmax should be < boxsize");
		}
		
		if(naxa*naxb>NRMAX)
		{
			Error(err,"[Error]: check "+paramfile+", too big naxa*naxb; should be <="+conv(NRMAX)+".");
		}
	}
	
	
	if(VERSION=="LC_ellipses")
	{
		if(DCMIN<=0 or DCMAX<=0)
		{
			Error(err,"[Error]: check "+paramfile+", DCMIN and DCMAX should be >0.");
		}
		
		if(DCMIN>=DCMAX)
		{
			Error(err,"[Error]: check "+paramfile+", DCMAX should be >= DCMIN.");
		}
	}
	
	
	if(VERSION=="angular")
	{
		if(Areaf<0 and Areaf!=-1)
		{
			Error(err,"[Error]: check "+paramfile+", Areaf should be either >0 or -1.");
		}
		
		if(Areaf>FULLSKYDEG2)
		{
			Error(err,"[Error]: check "+paramfile+", Areaf can not be larger than entire sky: "+conv(FULLSKYDEG2)+" deg2.");
		}
		
		if(ramin<-360. or ramin>360. or ramax<-360. or ramax>360.)
		{
			Error(err,"[Error]: check "+paramfile+", ramin,ramax should not exceed [-360,360] deg");
		}
		
		if(decmin<-90. or decmin>90. or decmax<-90. or decmax>90.)
		{
			Error(err,"[Error]: check "+paramfile+", decmin,decmax should not exceed [-90,90] deg");
		}
	}
	
	
	

	
	
	
	return err;
}