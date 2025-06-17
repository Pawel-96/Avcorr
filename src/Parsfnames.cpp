#include "Parsfnames.h"


int err_common=0; //error in common params
int err_ver=0; //error in params from one of versions
int err_conv=0; //error in conversion from string


const string VERSION=Get_parameter(paramfile,"VERSION",err_common)[0];	    //code version: angular/BOX/BOX_ellipses/LC_ellipses

const int USE_HDF5=conv_check(Get_parameter(paramfile,"USE_HDF5",err_common)[0],err_conv); 	//HDF5 [1] or ASCII [0] input data
const string POS_DSET=Get_parameter(paramfile,"POS_DSET",err_common)[0];		//dataset with positions (if USE_HDF5=1)
const vector<string> cols_pos=Get_parameter(paramfile,"cols_pos",err_common);  //columns with positions (if USE_HDF5=0)

//********************common parameters:********************
const int moment_min=2;                                 			//must be 2
const int moment_max=9;                                			    //must be 9
const double kappa=conv_check(Get_parameter(paramfile,"kappa",err_common)[0],err_conv);	    //multiplicity factor for number of circles
const int Cmin=conv_check(Get_parameter(paramfile,"Cmin",err_common)[0],err_conv);       //minimum number of circles drawn
const int Cmax=conv_check(Get_parameter(paramfile,"Cmax",err_common)[0],err_conv);       //maximum number of circles drawn
const int nreals=conv_check(Get_parameter(paramfile,"nreals",err_common)[0],err_conv);	    //number of independent randoms sets (randoms must be at least C*nreals)
const int NPIX_min=50; 									            //minimum pixelization [NEW]
const int Mvar_Poisson=1000;										//probing for Poisson error
const int Recalc=conv_check(Get_parameter(paramfile,"Recalc",err_common)[0],err_conv);        //if =1, calculating moments again
const int Clean=conv_check(Get_parameter(paramfile,"Clean",err_common)[0],err_conv);        //if =1, cleaning in the end
const int ErrP=conv_check(Get_parameter(paramfile,"ErrPoisson",err_common)[0],err_conv);        //compute Poisson errors? [0/1]
const int combine_reals=conv_check(Get_parameter(paramfile,"Combine_reals",err_common)[0],err_conv); //[0/1] are Data/* files models with different reals?
const string real_template=Get_parameter(paramfile,"Real_template",err_common)[0]; //how are realisations marked in names, e.g. _BOX*_

string EXT="";
string PATH="";
const vector<string> Model=Conditional_modelreading("Data/",paramfile,EXT,"Datafiles",PATH);
const int nmodels=Model.size();

const int Random_provided=conv_check(Get_parameter(paramfile,"Random_provided",err_common)[0],err_conv);  //random file provided? [0/1]
const string Random_file=Get_parameter(paramfile,"Random_file",err_common)[0]; //random file name

//********************parameters - angular:********************
const double ramin=conv_check(Get_parameter(paramfile,"ramin",err_ver)[0],err_conv); 		//[deg]catalog rightascension lower range
const double ramax=conv_check(Get_parameter(paramfile,"ramax",err_ver)[0],err_conv); 		//[deg]catalog rightascension lower range
const double decmin=conv_check(Get_parameter(paramfile,"decmin",err_ver)[0],err_conv);		//[deg]catalog declination lower range
const double decmax=conv_check(Get_parameter(paramfile,"decmax",err_ver)[0],err_conv); 	//[deg] catalog declination upper range

const double Areaf_read=conv_check(Get_parameter(paramfile,"Areaf",err_ver)[0],err_conv);
const double Areaf=Areaf_read==-1 ? Area_spherical(ramin,ramax,decmin,decmax)*str2deg2 : Areaf_read; //catalog sky area [deg2]

//********************parameters - angular and BOX:********************
const double Rmin=conv_check(Get_parameter(paramfile,"Rmin",err_ver)[0],err_conv); 		//[deg] smallest angular scale considered
const double Rmax=conv_check(Get_parameter(paramfile,"Rmax",err_ver)[0],err_conv); 		//[deg] biggest angular scale considered
const int nR=conv_check(Get_parameter(paramfile,"nR",err_ver)[0],err_conv); 				//number of angular scales considered
const double qR=pow(10.,1.*log10(Rmax/Rmin)/nR);        			//R multiplicity factor

//********************parameters - BOX and BOX_ellipses:********************
const double boxsize=conv_check(Get_parameter(paramfile,"Boxsize",err_ver)[0],err_conv); 	//box size in the same units as random spheres


//********************parameters - BOX_ellipses and LC_ellipses:********************
const double axamin=conv_check(Get_parameter(paramfile,"axamin",err_ver)[0],err_conv); 	//smallest semi-major axis
const double axamax=conv_check(Get_parameter(paramfile,"axamax",err_ver)[0],err_conv); 	//biggest semi-major axis
const double axbmin=conv_check(Get_parameter(paramfile,"axbmin",err_ver)[0],err_conv); 	//smallest semi-minor axis
const double axbmax=conv_check(Get_parameter(paramfile,"axbmax",err_ver)[0],err_conv); 	//biggest semi-minor axis
const int naxa=conv_check(Get_parameter(paramfile,"naxa",err_ver)[0],err_conv); 			//number of semi-major axes
const int naxb=conv_check(Get_parameter(paramfile,"naxb",err_ver)[0],err_conv); 			//number of semi-minor axes
const double qaxa=pow(10.,1.*log10(axamax/axamin)/naxa);        	//semi-major axis multiplicity factor
const double qaxb=pow(10.,1.*log10(axbmax/axbmin)/naxb);        //semi-minor axis multiplicity factor
const int n_allax=naxa*naxb;										//number of grid points

//********************parameters - LC_ellipses:********************
const double DCMIN=conv_check(Get_parameter(paramfile,"DCMIN",err_ver)[0],err_conv);		//minimum distance analyzed
const double DCMAX=conv_check(Get_parameter(paramfile,"DCMAX",err_ver)[0],err_conv);		//maximum distance analyzed

/*-------------------------------------------------------------------------------------------------------------------------*/
const int ncols_outfile=ErrP==0?17:25; //number of columns (check Args2grid)








string Fin(string model) //input
{
    return PATH+model+EXT;
}




string Mout_onerad(string model, int nnsize) //moment output for one radius
{
	return "Results/"+model+"_moments_nnR_"+conv(nnsize)+".txt";
}




string Merrout(string model) //output for one model and all moments
{
	return "Results/"+model+"_merrout.txt";
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
        lfile<<"Data ["<<nmodels<<"]:"<<endl;
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
			lfile<<"Scales: "<<nR<<" in range: ["<<Rmin<<":"<<Rmax<<"]"<<endl;
		}
		
		if(VERSION=="BOX" or VERSION=="BOX_ellipses")
		{
			lfile<<"Boxsize: "<<boxsize<<endl;
		}
		
		if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses")
		{
			lfile<<"LOS parallel axis scales: "<<naxa<<" in range: ["<<axamin<<":"<<axamax<<"]"<<endl;
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







void Error(int rank, int &err, string msg) //raising error with message
{
	err=1;
	if(rank==0)
	{
		cerr<<msg<<endl;
		Writelog(msg);
	}
	
	return;
}






//pointing errors in paramfile
int Error_param(int rank)
{
	int err=0;
	string msg="";
	
	
	if(err_common!=0) //one or more of common parameters are not specified
	{
		Error(rank,err,"[Error]: check "+paramfile+", "+conv(err_common)+" common parameter(s) not specified. How dare you!");
	}
	
	
	if(err_ver!=0) //one or more of common parameters are not specified
	{
		Error(rank,err,"[Error]: check "+paramfile+", "+conv(err_ver)+" VERSION-dependent parameter(s) not specified.\
		Some may not be needed for current VERSION, but better keep them in "+paramfile);
	}
	
	if(err_conv!=0) //numeric parameters are empty or strings
	{
		Error(rank,err,"[Error]: check "+paramfile+", "+conv(err_conv)+" numeric parameter(s) not specified/are strings");
	}
	
	
	
	if(VERSION!="angular" and VERSION!="BOX" and VERSION!="BOX_ellipses" and VERSION!="LC_ellipses")
	{
		Error(rank,err,"[Error]: check "+paramfile+", VERSION should be angular/BOX/BOX_ellipses/LC_ellipses.");
	}
	
	
	if(USE_HDF5!=0 and USE_HDF5!=1)
	{
		Error(rank,err,"[Error]: check "+paramfile+", USE_HDF5 must be either 0 or 1.");
	}
	
	if(nreals<3)
	{
		Error(rank,err,"[Error]: check "+paramfile+", nreals should be >=3 for reliable results");
	}
	
	if(kappa>.5 or kappa<=0)
	{
		Error(rank,err,"[Error]: check "+paramfile+", kappa should be in range (0,0.5] for reliable results");
	}
	
	
	if(Cmin<=0)
	{
		Error(rank,err,"[Error]: check "+paramfile+", Cmin should be positive");
	}
	
	if(Cmin>=Cmax)
	{
		Error(rank,err,"[Error]: check "+paramfile+", Cmin should be < Cmax");
	}
	
	
	if(Cmax>=CMAX)
	{
		Error(rank,err,"[Error]: check "+paramfile+", Cmax should be < "+conv(CMAX));
	}
	
	
	if(Random_provided==1)
	{
		vector<string> frandoms;
		Get_fnames(frandoms,"Randoms");
		int nfrandoms=frandoms.size();
		
		if((Random_file=="*" and nfrandoms<nmodels) or nfrandoms==0)
		{
			Error(rank,err,"[Error]: check "+paramfile+" or Randoms/ directory: some randoms may be missing");
		}
		
		if(Random_file!="*" and Fexist(Random_file)==0)
		{
			Error(rank,err,"[Error]: check "+paramfile+" or Randoms/ directory: Randoms/"+Random_file+" does not exist");
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
					Error(rank,err,"[Error]: check "+paramfile+", file: "+fr+" does not exist.");
				}
				else if(USE_HDF5==0)
				{
					ncols=Fncols(fr);
					if((VERSION=="angular" and ncols!=2) or (VERSION!="angular" and ncols!=3))
					{
						Error(rank,err,"[Error]: file: "+fr+" has wrong number of columns, ");
					}
				}
			}
		}
	}
	
	if(Random_provided!=0 and Random_provided!=1)
	{
		Error(rank,err,"[Error]: check "+paramfile+", Random_provided should be either 0 or 1.");
	}
	
	
	if(Recalc!=0 and Recalc!=1)
	{
		Error(rank,err,"[Error]: check "+paramfile+", Recalc should be either 0 or 1");
	}
	
	
	if(Clean!=0 and Clean!=1)
	{
		Error(rank,err,"[Error]: check "+paramfile+", Clean should be either 0 or 1");
	}
	
	
	if(ErrP!=0 and ErrP!=1)
	{
		Error(rank,err,"[Error]: check "+paramfile+", ErrPoisson should be either 0 or 1");
	}
	
	
	
	if(combine_reals!=0 and combine_reals!=1)
	{
		Error(rank,err,"[Error]: check "+paramfile+", Combine_reals should be either 0 or 1");
	}
	
	
	
	
	
	if(USE_HDF5==0) //only for ASCII: columns beyond data file(s)
	{
		if(VERSION=="angular" and cols_pos.size()!=2) //inappriopriate number of cols specified
		{
			err=1;
			Error(rank,err,"[Error]: check "+paramfile+", cols_pos should have 2 columns in VERSION=angular.");
		}
		
		if(VERSION!="angular" and cols_pos.size()!=3) //inappriopriate number of cols specified
		{
			err=1;
			Error(rank,err,"[Error]: check "+paramfile+", cols_pos should have 3 columns in this VERSION.");
		}
		
		int ncols=Fncols(Fin(Model[0]));
		
		if(cols_pos.size()==3)
		{
			if(stoi(cols_pos[0])<0 or stoi(cols_pos[0])>=ncols
			or stoi(cols_pos[1])<0 or stoi(cols_pos[1])>=ncols
			or stoi(cols_pos[2])<0 or stoi(cols_pos[2])>=ncols)
			{
				err=1;
				if(ncols!=-1) //if ncols==-1, first file doesnt exist, another error shows up
				{
					Error(rank,err,"[Error]: check "+paramfile+", cols_pos exceeding the file - should be in range [0,"+conv(ncols-1)+"].");
				}
			}
		}
		
		if(cols_pos.size()==2)
		{
			if(stoi(cols_pos[0])<0 or stoi(cols_pos[0])>=ncols
			or stoi(cols_pos[1])<0 or stoi(cols_pos[1])>=ncols)
			{
				err=1;
				if(ncols!=-1) //if ncols==-1, first file doesnt exist, another error shows up
				{
					Error(rank,err,"[Error]: check "+paramfile+", cols_pos exceeding the file - should be in range [0,"+conv(ncols-1)+"].");
				}
			}
		}
			
	}
	
	
	if(nmodels==0) //Data/ directory empty
	{
		Error(rank,err,"[Error]: check "+paramfile+", no Datafiles specified or Data/ empty.");
	}
	
	
	string modelsymbol=Get_parameter(paramfile,"Datafiles",err)[0];
	
	if(modelsymbol!="*") //file names written directly in paramfile and some doesnt exist
	{
		for(int i=0;i<nmodels;++i)
		{
			if(Fexist(Fin(Model[i]))==0)
			{
				Error(rank,err,"[Error]: check "+paramfile+", file: "+Fin(Model[i])+" does not exist.");
			}
		}
	}
	
	
	
	
	if(combine_reals==1)
	{
		if(real_template.size()==0)
		{
			Error(rank,err,"[Error]: check "+paramfile+", Combine_reals=1 is used, but Real_template is empty");
		}
		else
		{
			for(int i=0;i<nmodels;++i)
			{
				if(Model[i].find(Replace_string(real_template,"*","")) == string::npos)
				{
					Error(rank,err,"[Error]: check "+paramfile+", Real_template not found in datafile names");
				}
			}
		}
	
	}
	
	
	
	
	
	if(VERSION=="BOX" or VERSION=="BOX_ellipses")
	{
		if(boxsize<=0)
		{
			Error(rank,err,"[Error]: check "+paramfile+", Boxsize should be >0.");
		}
	}
	
	
	if(VERSION=="angular" or VERSION=="BOX")
	{
		if(Rmin<=0 or Rmax<=0 or nR<=0)
		{
			Error(rank,err,"[Error]: check "+paramfile+", Rmin, Rmax and nR should be >0.");
		}
		
		
		if(Rmin>=Rmax)
		{
			Error(rank,err,"[Error]: check "+paramfile+", Rmax should be >= Rmin.");
		}
		
		if(Rmax>=boxsize)
		{
			Error(rank,err,"[Error]: check "+paramfile+", Rmax should be < boxsize");
		}
		
		if(nR>NRMAX)
		{
			Error(rank,err,"[Error]: check "+paramfile+", too big nR; should be <="+conv(NRMAX)+".");
		}
	}
	
	
	if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses")
	{
		if(axamin<=0 or axamax<=0 or axbmin<=0 or axbmax<=0 or naxa<=0 or naxb<=0)
		{
			Error(rank,err,"[Error]: check "+paramfile+", axamin/max, axbmin/max, naxa/b should be >0.");
		}
		
		if(axamin>=axamax or axbmin>=axbmax)
		{
			Error(rank,err,"[Error]: check "+paramfile+", axamax should be >= axamin, same with axb.");
		}
		
		if(axamax>=boxsize or axbmax>=boxsize)
		{
			Error(rank,err,"[Error]: check "+paramfile+", axamax and axbmax should be < boxsize");
		}
		
		if(naxa*naxb>NRMAX)
		{
			Error(rank,err,"[Error]: check "+paramfile+", too big naxa*naxb; should be <="+conv(NRMAX)+".");
		}
	}
	
	
	if(VERSION=="LC_ellipses")
	{
		if(DCMIN<=0 or DCMAX<=0)
		{
			Error(rank,err,"[Error]: check "+paramfile+", DCMIN and DCMAX should be >0.");
		}
		
		if(DCMIN>=DCMAX)
		{
			Error(rank,err,"[Error]: check "+paramfile+", DCMAX should be >= DCMIN.");
		}
	}
	
	
	if(VERSION=="angular")
	{
		if(Areaf<0 and Areaf!=-1)
		{
			Error(rank,err,"[Error]: check "+paramfile+", Areaf should be either >0 or -1.");
		}
		
		if(Areaf>FULLSKYDEG2)
		{
			Error(rank,err,"[Error]: check "+paramfile+", Areaf can not be larger than entire sky: "+conv(FULLSKYDEG2)+" deg2.");
		}
		
		if(ramin<-360. or ramin>360. or ramax<-360. or ramax>360.)
		{
			Error(rank,err,"[Error]: check "+paramfile+", ramin,ramax should not exceed [-360,360] deg");
		}
		
		if(decmin<-90. or decmin>90. or decmax<-90. or decmax>90.)
		{
			Error(rank,err,"[Error]: check "+paramfile+", decmin,decmax should not exceed [-90,90] deg");
		}
	}
	
	
	return err;
}