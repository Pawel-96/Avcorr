//****************************************************************************************//
//                                                                                        //
// Copyright 2025 Paweł Drozda                                                            //
//                                                                                        //
// This file is part of Avcorr.                                                           //
//                                                                                        //
//To the extent possible under law, the author has dedicated this software	          //
// to the public domain under the Creative Commons CC0 1.0 Universal License.             //
//                                                                                        //
// Avcorr is distributed in the hope that it will be useful, but                          //
// WITHOUT ANY WARRANTY; without even the implied warranty of                             //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                   //
//                                                                                        //
// You should have received a copy of the CC0-1.0 License                                 //
// along with Avcorr.  If not, see <ttps://creativecommons.org/publicdomain/zero/1.0/>.   //
//                                                                                        //
//                                                                                        //
//****************************************************************************************//

#include "lib/RWMathStat.h"
#include "mpi.h" //parallel computing
#include "src/CIC.h"
#include "src/Combine.h"
#include "src/Moments.h"
#include "src/Parsfnames.h"
#include "lib/Files_HDF5.h"



//optimal pixelization
int Npix(vector<double> &Nav_circ, vector<double> &Nav_pix, int N, double R, double axa, double axb)
{
    double this_nav_circ=0.; //for optimization
	int loc_opt,npix=0;
	if(VERSION=="angular")
	{
		this_nav_circ=N*2.*M_PI*(1.-cos(R*deg2rad))/(Areaf/str2deg2); //average number of objects in circle
		if(this_nav_circ<Nav_circ[0]){npix=sqrt(1.*N*Area_spherical(ramin,ramax,decmin,decmax)*str2deg2/(Nav_pix[0]*Areaf));}
        if(this_nav_circ>Nav_circ[Nav_circ.size()-1]){npix=sqrt(1.*N*Area_spherical(ramin,ramax,decmin,decmax)*str2deg2/(Nav_pix[Nav_circ.size()-1]*Areaf));}
        if(this_nav_circ>=Nav_circ[0] and this_nav_circ<=Nav_circ[Nav_circ.size()-1])
        {
            loc_opt=1.*Nav_circ.size()*(log10(this_nav_circ)-log10(Nav_circ[0])) / (log10(Nav_circ[Nav_circ.size()-1])-log10(Nav_circ[0])); //NAV_CIRC IS IN LOGSPACE
            npix=pow(1.*N*Area_spherical(ramin,ramax,decmin,decmax)*str2deg2/(Nav_pix[loc_opt]*Areaf),0.5);
        }
		
	}
	else
	{
		if(VERSION=="BOX"){this_nav_circ=N*(4.0/3.0)*M_PI*pow(R,3)/pow(1.0*boxsize,3);}
		if(VERSION=="BOX_ellipses"){this_nav_circ=N* 4./3. *M_PI*pow(axb,2)*axa/pow(1.0*boxsize,3);}
		if(VERSION=="LC_ellipses"){this_nav_circ=N*pow(axb,2)*axa/(pow(1.*DCMAX,3)-pow(1.*DCMIN,3));}
		
		if(this_nav_circ<Nav_circ[0]){npix=pow(1.*N/Nav_pix[0],1.0/3.0);}
        if(this_nav_circ>Nav_circ[Nav_circ.size()-1]){npix=pow(1.*N/Nav_pix[Nav_circ.size()-1],1.0/3.0);}
        if(this_nav_circ>=Nav_circ[0] and this_nav_circ<=Nav_circ[Nav_circ.size()-1])
        {
            loc_opt=1.*Nav_circ.size()*(log10(this_nav_circ)-log10(Nav_circ[0])) / (log10(Nav_circ[Nav_circ.size()-1])-log10(Nav_circ[0])); //NAV_CIRC IS IN LOGSPACE
            npix=pow(1.*N/Nav_pix[loc_opt],1.0/3.0);
        }
	}
	
	npix=max(npix,NPIX_min);
	return npix;
}





int Results_onemodel_oneradius(string model, int nnsize, vector<double> &Nav_circ, vector<double> &Nav_pix,
vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand) //calculating moments, reducing, correcting
{
	int msg=0;
    vector<double> x,y,z; //for angular: x,y=ra,dec
    vector<vector<vector<double> > > pix; //for pixelization
    string fin,moments_output,ffound,ffoundhst;
    int npix=0.,N;
    double R=0,axa,axb; //radii/ellipse sizes ;for different cases

    pix.clear();x.clear();y.clear();z.clear();
    fin=Fin(model);
    moments_output=Mout_onerad(model,nnsize); //moment output name
    ffound=Replace_string(moments_output,".txt","_found.fnd");
	ffoundhst=Replace_string(ffound,".fnd",".fndhst");

    cout<<"Catalog: "<<fin<<endl;
    if(Fexist(ffound)==0 and Fexist(ffoundhst)==0) //pixelizing only if counts not done yet
    {
		if(USE_HDF5==0) //ASCII format
		{
			if(VERSION=="angular"){Fread<double>(fin,{&x,&y},{stoi(cols_pos[0]),stoi(cols_pos[1])});}
			else {Fread<double>(fin,{&x,&y,&z},{stoi(cols_pos[0]),stoi(cols_pos[1]),stoi(cols_pos[2])});}
		}
		if(USE_HDF5==1) //HDF5
		{
			auto dcoords=Read_HDF5_dataset(fin,POS_DSET,2,msg);
			
			if(msg!=0)
			{
				Writelog("Check "+paramfile+", dataset "+POS_DSET+" does not exist in "+fin+", stopping:/");
				return 1;
			}
			vector<vector<double>> coords=get<vector<vector<double>> > (dcoords);
			
			//data needs to be transposed
			if((VERSION=="angular" and coords[0].size()==2) or (VERSION!="angular" and coords[0].size()==3))
			{
				Transpose<double>(coords);
			}
			
			x=coords[0];
			y=coords[1];
			if(VERSION!="angular"){z=coords[2];}
			
		}
		
        N=x.size();
		cout<<"Pixelizing "<<endl;

        if(VERSION=="angular" or VERSION=="BOX"){R=Rmin*pow(1.0*qR,nnsize);}
		if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses"){Find_ellipses_axes(nnsize,axa,axb);}
		
		npix=Npix(Nav_circ,Nav_pix,N,R,axa,axb); //optimal pixelization
		
		if(VERSION=="angular"){Pixelize<double>({ramin,ramax,decmin,decmax},{npix,npix},{x,y},pix,"2D_sky"); }
		if(VERSION=="BOX" or VERSION=="BOX_ellipses"){Pixelize<double>({0,boxsize,0,boxsize,0,boxsize},{npix,npix,npix},{x,y,z},pix,"3D"); }
		if(VERSION=="LC_ellipses"){Pixelize<double>({-1.*DCMAX,DCMAX,-1.*DCMAX,DCMAX,-1.*DCMAX,DCMAX},{npix,npix,npix},{x,y,z},pix,"3D");}
        
    }
    cout<<"Calculating moments: "<<endl;
	if(VERSION=="angular"){CIC_angular(pix,moments_output,model,npix,nnsize,Xcn,Ycn,Zcn,nrand);}
	if(VERSION=="BOX"){CIC_BOX(pix,moments_output,model,npix,nnsize,Xcn,Ycn,Zcn,nrand);}
	if(VERSION=="BOX_ellipses"){CIC_BOX_ellipses(pix,moments_output,model,npix,nnsize,Xcn,Ycn,Zcn,nrand);}
	if(VERSION=="LC_ellipses"){CIC_LC_ellipses(pix,moments_output,model,npix,nnsize,Xcn,Ycn,Zcn,nrand);}
    
	return msg;
}







//running ii-th work (out of iimax)
void Run(int ii, int iimax, vector<double> &Nav_circ, vector<double> &Nav_pix, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int nrand)
{
	int nnR=0,nnallax=0,dset_err=0;
	int model_no=ii/(iimax/nmodels); //number of model
	if(VERSION=="angular" or VERSION=="BOX"){nnR=ii%nR;} //number of circle radius
	if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses"){nnallax=ii%n_allax;} //current argument grid point
        

	cout<<"**********Model********** "<<model_no+1<<"/"<<Model.size()<<endl;
	if(VERSION=="angular" or VERSION=="BOX"){dset_err=Results_onemodel_oneradius(Model[model_no],nnR,Nav_circ,Nav_pix,Xcn,Ycn,Zcn,nrand);}
	if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses"){dset_err=Results_onemodel_oneradius(Model[model_no],nnallax,Nav_circ,Nav_pix,Xcn,Ycn,Zcn,nrand);}
        
	if(dset_err==1)
	{
		Error(0,dset_err,"Data reading error(s) occurred, stopping:/"); //rank=0 set just to print
		Writelog("Proceeding with another work, but this will raise error at postprocessing");
		return;
	}
	
	double progress=Progress(iimax);
	cout<<"<----Overall progress: "<<progress<<"----->"<<endl;
	Writelog("*** Overall progress: "+conv(progress)+" ***");

	return ;
}







int main(int argc, char *argv[])
{
	int rank, comm_size;//,iimax=0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	
	int err=Error_param(rank); //error in parameters
	if(err==1)
	{
		Error(rank,err,"Errors occurred, stopping:/");
		return 0;
	}
	
	int nrand;
	string fnav=VERSION=="angular"?"Optimal_pix_angular.txt":"Optimal_pix_BOX.txt"; //file for optimization in pixelizing
	int iimax=VERSION=="angular" || VERSION=="BOX"? nR*nmodels : n_allax*nmodels;
    vector<double> Nav_circ,Nav_pix,Xcn,Ycn,Zcn; //X/Y/Zcn - for angular: xcn,ycn=racn,deccn, for BOX/3D: x,y,z centers
	
	if(rank==0)
	{
		Writelog("","intro"); //writing main information about the run
		Writelog("Running with "+conv(comm_size)+" threads");
		Writelog("Preparing environment");		
	}
    srand(time(NULL)+rank);

	Fread<double>(fnav,{&Nav_circ,&Nav_pix},{0,1}); //reading the data
	Read_randoms(Xcn,Ycn,Zcn,"",nrand); //reading randoms based on option(s)
	
	//----------------------------------------------------------------MPI routine-------------------------------------------------
	int msg_recv=0,msg_send; //messages for proces communication
	int STOP=-1,TAG=1;
	
	if(comm_size<3) //small number of threads, simple work division
	{
		for(int i=0;i<iimax;++i) //computation for all moments and radii
		{
			if(i%comm_size!=rank){continue;} //every comm_size-th assigned to this rank
			Run(i,iimax,Nav_circ,Nav_pix,Xcn,Ycn,Zcn,nrand);
			//in case of data reading error, just getting to next work, but error is raised and code considers that in postprocessing
		}
	}
	else //optimized work division
	{
		if(rank==0) //master thread, assigning job to other ones
		{
			for(int i=1;i<comm_size;++i) //firstly assigning the work to all workers
			{
				msg_send=i-1; //work ID to be done (master just assigns)
				MPI_Send(&msg_send, 1, MPI_INT, i, TAG, MPI_COMM_WORLD); //informing i-th worker worker which work to do
			}
			int first_notstarted=comm_size-1; //lowest no. of not-yet started work
			
			int sender;
			while(first_notstarted<iimax) //repeating until everything is done
			{
				//receiving information that one worker (any_source) finished assigned job
				MPI_Recv(&msg_recv, 1, MPI_INT, MPI_ANY_SOURCE, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				sender = msg_recv; //which worker sent that
				msg_send=first_notstarted; //message informing which work to do now
				
				//sending which work this worker should do now (first unstarted)
				MPI_Send(&msg_send, 1, MPI_INT, sender, TAG, MPI_COMM_WORLD) ;
				first_notstarted+=1; //now first unstarted job is the next one
			}
			
			//everything is done, sending message to stop:
			for(int i=1;i<comm_size;++i)
			{
				MPI_Send(&STOP, 1, MPI_INT, i, TAG, MPI_COMM_WORLD);
			}
		}
		
		if(rank!=0) //worker
		{
			msg_send=rank;
			while(true) //repeating until everything is done (received work ID exceed limit)
			{
				 //receiving information from master (rank=0) about which work to do (msg_recv)
				MPI_Recv(&msg_recv, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				if (msg_recv==STOP) {break;}  //everything is done
				
				Run(msg_recv,iimax,Nav_circ,Nav_pix,Xcn,Ycn,Zcn,nrand); //doing assigned work
				//in case of data reading error, just getting to next work, but error is raised and code considers that in postprocessing
				
				MPI_Send(&msg_send, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD); //informing master thread that work is done (sending rank)
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

    MPI_Finalize(); //ending MPI routine

    if(rank==0) //combining results, done only in one thread
    {
		err=Merge();
		if(err==1){return 0;} //some outputs does not exist (error raised before, in data reading)
		if(VERSION=="BOX_ellipses" or VERSION=="LC_ellipses"){Args2grid();} //converting 1D arguments back to 2D grid in results
		
		if(combine_reals==1) //combining realisation into one results
		{
			Writelog("Combining realisations");
			cout<<"Combining realisations"<<endl;
			vector<string> IndpModels=Get_independent_models(Model);
			int n_indpmodels=IndpModels.size();
			
			for(int i=0;i<n_indpmodels;++i){Combine_reals(IndpModels[i]);}
		}
		Writelog("","outro");
    }

	return 0;
}
