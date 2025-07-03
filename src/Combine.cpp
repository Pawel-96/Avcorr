#include "Combine.h"



int Merge()
{
	string f,foutXi,foutS,line; //files with Xi _Jand S_J
	int NR; //number of radii
	
	if(VERSION=="angular" or VERSION=="BOX")
	{
		NR=nR;
	}
	else {NR=n_allax;}
	
	for(int i=0;i<nmodels;++i)
	{
		foutXi=Merrout(Model[i]);
		foutS=Replace_string(foutXi,"_merrout.txt","_merrout_Sn.txt");
		ofstream fXi(foutXi.c_str());
		ofstream fS(foutS.c_str());
		Writelog("Merging -> "+foutXi);
		
		for(int j=0;j<NR;++j) //every radius/scale
		{
			f=Mout_onerad(Model[i],j); //file with j-th radii result (XI_J)
			if(Fexist(f)==0){return 1;} //no output (previously reported error in data reading)
			
			ifstream ffxi(f.c_str());
			getline(ffxi,line);
			ffxi.close();
			fXi<<line<<endl;
			if(Clean==1){remove(f.c_str());} //cleaning
			
			f=Replace_string(Mout_onerad(Model[i],j),".txt","_Sn.txt"); //file with j-th radii result (S_J)
			if(Fexist(f)==0){return 1;} //no output (previously reported error in data reading)
			
			ifstream ffs(f.c_str());
			getline(ffs,line);
			ffs.close();
			fS<<line<<endl;
			if(Clean==1){remove(f.c_str());} //cleaning
		}
		
		fXi.close();
		fS.close();
	}

	return 0;
}






//converting 1D arguments back to 2D grid in results (one file)
void Args2grid_onefile(string fname)
{
	if(Fexist(fname)==0){return;} //no output (previously reported error in data reading)
	
	string **data=Empty<string>(ncols_outfile,n_allax); //to avoid nan problems
	string ftemp=Random_str(15);
    double axa,axb;
	
	Fread<string>(fname,data); //reading from file
	ofstream ft(ftemp.c_str());

	for(int j=0;j<n_allax;++j) //writing each line -> changing 1D args -> 2D grid
	{
		Find_ellipses_axes(j,axa,axb); //lines are sorted in Merge()
			
		ft<<axa<<" "<<axb<<" ";
		for(int k=1;k<ncols_outfile;++k) //writing next lines
		{
			ft<<data[k][j]<<" ";
		}
		ft<<endl;
	}
	ft.close();
	rename(ftemp.c_str(),fname.c_str());

	return;
}






void Args2grid() //converting 1D arguments back to 2D grid in results
{
	string fout,fout_sn;
	for(int i=0;i<nmodels;++i)
	{
		fout=Merrout(Model[i]);
		fout_sn=Replace_string(fout,".txt","_Sn.txt");
		Args2grid_onefile(fout);
		Args2grid_onefile(fout_sn);
	}
	return;
}




//extracting realisation number from modelname based on real_template
int Realisation(string modelname)
{
	int npos_templ=real_template.find("*"); //position of * in real template
	int agree=0,position=0,lastcheck;
	int templsize=real_template.size();
	string nrealisation;
	
	while(agree<npos_templ) //position of starting realisation number
	{
		if(agree>0 and modelname[position]!=real_template[agree])
		{
			agree=0; //not total match
		}
		if(modelname[position]==real_template[agree])
		{
			agree+=1;
		}
		
		position+=1;
	}
	
	if(npos_templ==templsize-1){lastcheck=modelname.size();}
	else{lastcheck=modelname.size()-1;} // "*" not at end
	for(int i=position;i<lastcheck;++i)
	{
		nrealisation+=modelname[i];
		if(npos_templ<templsize-1){if(modelname[i+1]==real_template[npos_templ+1]){break;}} //sign after * in template
	}
	return conv(nrealisation);
	
}





//getting independent models if combine_reals=1 (files are models with different reals)
vector<string> Get_independent_models(vector<string> models)
{
	vector<string> IndpModels;
	string modelname,partname,common=Replace_string(real_template,"*","N"); //partname - fragment with info about realisation
	int nnreal,nmodels_all=models.size();
	for(int i=0;i<nmodels_all;++i)
	{
		nnreal=Realisation(models[i]); //realisation of this file
		partname=Replace_string(real_template,"*",conv(nnreal)); //name part with realisation info
		modelname=Replace_string(models[i],partname,common); //removing info about realisation
		IndpModels.push_back(modelname); //putting name without realisation info
	}
	Unique(IndpModels); //remove repeating
	cout<<"***Independent models***"<<endl;
	int nindpmodels=IndpModels.size();
	for(int i=0;i<nindpmodels;++i){cout<<i<<" "<<IndpModels[i]<<endl;}
	
	return IndpModels;
}







//finding number of realisations (checking if file exists, if not, stopping)
int Maxreal(string modelname) //[modelname as common, without realisation info]
{
	string fname,common=Replace_string(real_template,"*","N");
	int ii=1;
	while(true)
	{
		fname=Merrout(Replace_string(modelname,common,Replace_string(real_template,"*",conv(ii+1))));
		if(Fexist(fname)==1){ii+=1;}
		else break;
	}

	return ii;
}






//combining different realisations of model modelname (with real_template: *->N)
void Combine_reals(string modelname)
{
	cout<<"Combining realisations: "<<modelname<<endl;
	vector<vector<float> > WJ_onereal,SJ_onereal; //values for one realization
	float WJ,SJ,WJ_err,SJ_err; //final values
	string common,fname; //modelnames (with different realisations for this model name), template
	common=Replace_string(real_template,"*","N");
	int norders=8,nscales; //number of orders and scales
	int NREALS=Maxreal(modelname); //found number of realizations
	
	//arguments for different version types:
	vector<double> arg1,arg2; //argument/grid (dependent on VERSION)
	
	fname=Merrout(Replace_string(modelname,common,Replace_string(real_template,"*","1")));
	
	if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses") //2D argument grid
	{
		Fread(fname,{&arg1,&arg2},{0,1});
		nscales=n_allax;
	}
	else //1D arguments
	{
		Fread(fname,{&arg1},{0});
		nscales=nR;
	}
	vector<vector<vector<float>>> WJ_all(norders, vector<vector<float>>(nscales,vector<float>(NREALS)));
	vector<vector<vector<float>>> SJ_all(norders, vector<vector<float>>(nscales,vector<float>(NREALS)));
	
	for(int i=0;i<NREALS;++i) //loop over realisations - collecting
	{
		fname=Merrout(Replace_string(modelname,common,Replace_string(real_template,"*",conv(i+1))));
		cout<<"Combining file: "<<fname<<endl;
		
		if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses")
		{
			Fread(fname,WJ_onereal,{2,4,6,8,10,12,14,16}); //collecting WJ_onereal values
			Fread(Replace_string(fname,"merrout.txt","merrout_Sn.txt"),SJ_onereal,{2,4,6,8,10,12,14,16}); //collecting SJ_onereal values
		}
		else
		{
			Fread(fname,WJ_onereal,{1,3,5,7,9,11,13,15}); //collecting WJ_onereal values
			Fread(Replace_string(fname,"merrout.txt","merrout_Sn.txt"),SJ_onereal,{1,3,5,7,9,11,13,15}); //collecting SJ_onereal values
		}
		
		
		for(int j=0;j<norders;++j)
		{
			for(int k=0;k<nscales;++k)
			{
				WJ_all[j][k][i]=WJ_onereal[j][k];
				SJ_all[j][k][i]=SJ_onereal[j][k];
			}
		}
		WJ_onereal.clear();
		SJ_onereal.clear();
	}

	ofstream WJfile(Merrout(modelname).c_str());
	ofstream SJfile(Replace_string(Merrout(modelname),"_merrout.txt","_merrout_Sn.txt").c_str());

	//saving output
	for(int i=0;i<nscales;++i)
	{
		if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses")
		{
			WJfile<<arg1[i]<<" "<<arg2[i]<<" ";
			SJfile<<arg1[i]<<" "<<arg2[i]<<" 1 0 ";
		}
		else
		{
			WJfile<<arg1[i]<<" ";
			SJfile<<arg1[i]<<" 1 0 ";
		}
		
		for(int j=0;j<norders;++j)
		{
			Clean_data<float>(WJ_all[j][i],"erase<=0");
			Clean_data<float>(SJ_all[j][i],"erase<=0");
			WJ=Sum(WJ_all[j][i])/(1.*WJ_all[j][i].size());
			WJ_err=Sigma<float>(WJ_all[j][i]);
			
			if(j>0)
			{
				SJ=Sum(SJ_all[j][i])/(1.*SJ_all[j][i].size());
				SJ_err=Sigma<float>(SJ_all[j][i]);
			}
			
			WJfile<<WJ<<" "<<WJ_err<<" ";
			if(j>0){SJfile<<SJ<<" "<<SJ_err<<" ";}
		}
		WJfile<<endl;
		SJfile<<endl;
	}
	
	WJfile.close();
	SJfile.close();
	return;
	
}
