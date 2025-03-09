#include "Parsfnames.h"


string Fin(string model) //input
{
    return "Data/"+model+".txt";
}




string Mout_onerad(string model, int nnsize) //moment output for one radius
{
    string fin=Fin(model);
    string mout=Replace_string(fin,".txt","_moments_nnR_"+conv(nnsize)+".txt");
    return Replace_string(mout,"Data/","Results/");
}




string Merrout(string model) //output for one model and all moments
{
    string f1=Replace_string(Fin(model),"Data/","Results/");
    return Replace_string(f1,".txt","_merrout.txt");
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
        lfile<<"Models ["<<nmodels<<"]:"<<endl;
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