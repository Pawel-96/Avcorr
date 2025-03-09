#include "Moments.h"




void Connect(double *M) //connected moments (array 2-9)
{
	double *pom=new double[10];
	for(int i=moment_min;i<=moment_max;++i){pom[i]=M[i];}
	
	M[4]=pom[4]-3.0*pow(pom[2],2);
    M[5]=pom[5]-10.0*pom[3]*pom[2];
    M[6]=pom[6]-15.0*pom[4]*pom[2]-10.0*pow(pom[3],2)+30.0*pow(pom[2],3);
    M[7]=pom[7]-21.0*pom[5]*pom[2]-35.0*pom[4]*pom[3]+210.0*pom[3]*pow(pom[2],2);
    M[8]=pom[8]-28.0*pom[6]*pom[2]-56.0*pom[5]*pom[3]-35.0*pow(pom[4],2)+\
    420.0*pom[4]*pow(pom[2],2)+560.0*pow(pom[3],2)*pom[2]-630.0*pow(pom[2],4);
	
    M[9]=pom[9]-36.0*pom[7]*pom[2]-84.0*pom[6]*pom[3]-126.0*pom[5]*pom[4]+\
    756.0*pom[5]*pow(pom[2],2)+2520.0*pom[4]*pom[3]*pom[2]+560.0*pow(pom[3],3)-\
    7560.0*pow(pom[3],2)*pow(pom[2],3);
	return;
}




void Correct_average(double *M, double mean)
{
	double *pom=new double[10];
	for(int i=moment_min;i<=moment_max;++i){pom[i]=M[i];}
	
	//correcting:
	M[2]=pom[2] -mean;
	M[3]=pom[3]-3.0*M[2] -mean;
	M[4]=pom[4]-7.0*M[2]-6.0*M[3] -mean;
	M[5]=pom[5]-15.0*M[2]-25.0*M[3]-10.0*M[4] - mean;
	M[6]=pom[6]-31.0*M[2]-90.0*M[3]-65.0*M[4]-15.0*M[5] - mean;
	M[7]=pom[7]-63.0*M[2]-301.0*M[3]-350.0*M[4]-140.0*M[5]-21.0*M[6] - mean;
	M[8]=pom[8]-127.0*M[2]-966.0*M[3]-1701.0*M[4]-1050.0*M[5]-266.0*M[6]-\
	28.0*M[7] - mean;
	
	M[9]=pom[9]-255.0*M[2]-3025.0*M[3]-7770.0*M[4]-6951.0*M[5]-2646.0*M[6]-\
	462.0*M[7]-36.0*M[8] - mean;
	
	//averaging
	for(int i=2;i<=9;++i){M[i]*=1.*pow(mean+1e-10,-1.*i);}
	return;
}




double Get_moments(double *arr, int C, double *M) //calculating moments from array, writing to M
{
	double mean=Sum(arr,C)/(1.*C),moment_value;
	for(int j=0;j<C;++j){arr[j]-=mean;} //centering

    for(int k=moment_min;k<moment_max+1;++k)
    {
        moment_value=0.0;
        for(int j=0;j<C;++j)
        {
            moment_value+=pow(arr[j],k);
        }
        moment_value/=1.*C;
        M[k]=moment_value;
    }
    return mean;
}







double *Poisson_errs(vector<int> &nc, vector<int> &count, vector<int> &countCDF) //calculating Poisson errors
{
	int CR=countCDF[countCDF.size()-1];
	double *M=new double[10]; //storing moments
	double *Mvar=new double[10]; 
	double *err_moments=new double[10];
	double *err=new double[10];
	double *centered=new double[CR];
	vector<int> list=Hist2list<int>(nc,count);
	int nlist=list.size();
	double *listed=new double[nlist]; for(int i=0;i<nlist;++i){listed[i]=1.*list[i];}
	
	double Mean=Get_moments(listed,CR,M);
	double sum1,sum2; //contribution from counts and mean for small error propagation
	int pow1,pow2;
	double *par=new double[2];
	double **vals=new double *[10];for(int i=0;i<10;++i){vals[i]=new double [Mvar_Poisson];}
	
	for(int i=0;i<CR;++i){centered[i]=1.*listed[i]-Mean;}
	for(int k=moment_min;k<=moment_max;++k)
	{
		sum1=0.;
		sum2=0.;
		pow1=2*(k-1);
		pow2=k-1;
		for(int i=0;i<CR;++i)
		{
			if(i%100==0){cout<<"\rErrP summing: "<<i+1<<"/"<<CR<<" k="<<k;}
			sum1+=pow(centered[i],pow1);
			sum2+=pow(centered[i],pow2);
		}
		err_moments[k]=k*pow(Mean*(sum1 + pow(sum2,2)/(1.*CR)),0.5)/(1.*CR);
		if(isinf(err_moments[k])==1){err_moments[k]=-1000000;}
	}
	//contribution from connection and averaging (no correction!)
	for(int i=0;i<Mvar_Poisson;++i) //probing
	{
		if(i%10==0){cout<<"\rPoisson sample: "<<i+1<<"/"<<Mvar_Poisson;}
		for(int j=moment_min;j<=moment_max;++j) //creating random scatter
		{
			par[0]=M[j];
			par[1]=err_moments[j];
			if(err_moments[j]==-1000000){Mvar[j]=M[j];}
			else Mvar[j]=Random_distr<double>(Gauss_max1,par,M[j]-5.*err_moments[j],M[j]+5.*err_moments[j]);
		}
		Connect(Mvar);
		for(int j=moment_min;j<=moment_max;++j){vals[j][i]=Mvar[j]*pow(Mean+1e-10,-1.*j);}
	}
	cout<<endl;
	
	for(int j=moment_min;j<=moment_max;++j){err[j]=Sigma(vals[j],Mvar_Poisson);}
	return err;	
}






void Get_subsample(vector<int> &nc, vector<int> &countCDF, double *subsample, int C_in, int C_out)
{
	int loc,nCDF=countCDF.size();
	for(int i=0;i<C_out;++i)
	{
		loc=Random(0,C_in);
		for(int j=0;j<nCDF;++j) //searching in CDF which element it is
		{
			if(loc<countCDF[j])
			{
				subsample[i]=nc[j];
				break;
			}
		}
	}
	return;
}











void Calculate_moments(string output, vector<int> &nc, vector<int> &count, double R) //calculating moments from data in *found array
{
	
	string output_SJ=Replace_string(output,".txt","_Sn.txt");
	if((Fexist(output)==1 and Fexist(output_SJ)==1) and Recalc==0){return;} //ending if done previously
	
	int CR=Sum(count); //all circles
	int C_subs=CR/nreals; //subsample size
	double *subsample=new double[C_subs]; //array storing one subsample
	double *M=new double[10]; //storing raw moments
	vector<vector<double>> WJ_all(10,vector<double>(nreals)), SJ_all(10,vector<double>(nreals)); //storing results from all realizations
	vector<double> WJ(10), SJ(10), WJ_err(10), SJ_err(10); //final results with errors
	double mean,errPoisson_Sn;
	Writelog("Sub-cutting -> "+output);
	
	vector<int> countCDF=Cumulative(count); //CDF for subsampling

	for(int i=0;i<nreals;++i) //subsampling
	{
		cout<<"\rSub-cutting: "<<i+1<<"/"<<nreals;
		Get_subsample(nc,countCDF,subsample,CR,C_subs); //writing subsample into subsample array
		mean=Get_moments(subsample,C_subs,M); //calculating moments on this subsample
		Connect(M); //connected moments
		Correct_average(M,mean); //correcting and averaging
		
		WJ_all[0][i]=0.;WJ_all[1][i]=0.;
		SJ_all[0][i]=0.;SJ_all[1][i]=0.;SJ_all[2][i]=1.;
		for(int j=2;j<10;++j)//collecting to big array
		{
			WJ_all[j][i]=M[j]; //Wj
			if(j>2){SJ_all[j][i]=M[j]*pow(M[2],1-j);} //Sj
		}
	}
	cout<<endl<<"Sub-cutted -> calculating dispersions for err"<<endl;
	Combine_results<double>(WJ_all,WJ,WJ_err,"Gauss","erase<=0");
	Combine_results<double>(SJ_all,SJ,SJ_err,"Gauss","erase<=0");
	cout<<"Dispersions calculated, writing output for current model & scale"<<endl;
	
    ofstream foutWJ(output.c_str()); //storing values for Wj
	ofstream foutSJ(output_SJ.c_str()); //storing values for Sj
    foutWJ<<R<<" ";
	foutSJ<<R<<" 1 0 ";
    for(int k=moment_min;k<=moment_max;++k)
    {
        foutWJ<<setprecision(15)<<WJ[k]<<" "<<WJ_err[k]<<" ";
		if(k>2){foutSJ<<setprecision(15)<<SJ[k]<<" "<<SJ_err[k]<<" ";}
    }
	
	if(ErrP==1) //computing Poisson errors?
	{
		cout<<"Getting Poisson err"<<endl;
		double *ErrPoisson;
		ErrPoisson=Poisson_errs(nc,count,countCDF);
		foutSJ<<"0 ";
		
		for(int k=moment_min;k<=moment_max;++k)
		{
			foutWJ<<ErrPoisson[k]<<" ";
			if(k>2)
			{
				if(ErrP==1){errPoisson_Sn=pow( pow(ErrPoisson[k]*pow(WJ[2],1.-k),2)+pow(ErrPoisson[2]*WJ[k]*(1.-k)*pow(WJ[2],-1.*k),2),0.5 );}
				if(ErrP==0){errPoisson_Sn=0.;}
				foutSJ<<setprecision(15)<<errPoisson_Sn<<" ";
			}
		}
	}
	
    foutWJ<<endl;
	foutSJ<<endl;
    foutWJ.close();
	foutSJ.close();
	cout<<"Output written: "<<output<<endl;
	Writelog("Done: -> "+output);
	cout<<endl;
	return;
}