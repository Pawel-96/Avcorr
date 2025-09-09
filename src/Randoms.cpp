#include "Randoms.h"




//creating temporary randoms in continuous field
void Continuous_randoms(string frand_temp, int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn)
{
	double rarand,decrand;
	double R=Rmin*pow(1.0*qR,nnR);
	//number of random points:
	int C=kappa*nreals*Area_spherical(Random_fullr[0],Random_fullr[1],Random_fullr[2],Random_fullr[3])/(2.0*M_PI*(1.0-cos(R*deg2rad)));
	
	Xcn.resize(C);
	Ycn.resize(C);
	Zcn.resize(C);
	
	if(Fexist(frand_temp)==0)
	{
		ofstream frandt(frand_temp.c_str());
		for(int i=0;i<C;++i)
		{
			Random_position<double>({Random_fullr[0],Random_fullr[1],Random_fullr[2],Random_fullr[3],R},rarand,decrand,"sky_border");
			frandt<<rarand<<" "<<decrand<<endl;
			Xcn[i]=rarand;
			Ycn[i]=decrand;
		}

		frandt.close();
	}
	else //file already exist
	{
		Xcn.clear();Ycn.clear();Zcn.clear();
		Fread<double>(frand_temp,{&Xcn,&Ycn},{0,1});
	}
	
	return;
}






//processing randoms for one scale
void Process_randoms_onerad(int rank, int comm_size,int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn,
vector<double> &Nav_circ, vector<double> &Nav_pix, int &nrand)
{
	//temporary randoms within continuous field defined by Random_fullr
	string frand_temp="Randoms/continuous_temp_"+conv(nnR)+".txt";
	string frand="Randoms/Random_proc_"+conv(nnR)+".txt";
	if(Fexist(frand)==1){return;}
	
	vector<double> x,y,z; //for data - here Random_file
	double av,sigma;
	int npix;
	vector<vector<vector<double> > > pix; //for pixelization
	
	if(VERSION=="angular")
	{
		Continuous_randoms(frand_temp,nnR,Xcn,Ycn,Zcn); //creating randoms on continuous field
		
		Read_randoms(rank,comm_size,x,y,z,"",nrand,nnR,Nav_circ,Nav_pix,-1); //reading Random_file into x,y
		
		npix=Npix(Nav_circ,Nav_pix,x.size(),Rmin*pow(1.0*qR,nnR),0,0); //optimal pixelization
		Pixelize<double>({ramin,ramax,decmin,decmax},{npix,npix},{x,y},pix,"2D_sky");
		
		//get average and sigma from continuous field:
		CIC_angular_rnd(pix,"Randoms/temporary",npix,nnR,Xcn,Ycn,Zcn,nrand,0,av,sigma); //counting - mode=0
			
		Xcn.clear();
		Ycn.clear();
		Zcn.clear();
		
		//selecting data from Random_file where counts> average - Random_sigmac (CIC with mode=3): [data already pixelized]
		CIC_angular_rnd(pix,frand,npix,nnR,x,y,z,nrand,1,av,sigma); //counting - mode=1;  position of circles the same as Random_file

		remove(frand_temp.c_str()); //removing
	}
	
	
	return;
}








void Create_processed_randoms(int rank, int comm_size, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn,
vector<double> &Nav_circ, vector<double> &Nav_pix)
{
	int nrand=0;
	if(rank==0)
	{
		cout<<"---Processing randoms for irregular mask---"<<endl;
		Writelog("---Processing randoms for irregular mask---");
	}
	
	//temporary solution------------------------------------------------
	for(int i=0;i<nR;++i) //for all scales
	{
		if(i%comm_size!=rank){continue;} //every comm_size-th assigned to this rank
		Process_randoms_onerad(rank,comm_size,i,Xcn,Ycn,Zcn,Nav_circ,Nav_pix,nrand);
	}
	return;
}










//angular CiC, but version especially for random processing
void CIC_angular_rnd(vector<vector<vector<double> > > &pix, string output, 
int npix_edge, int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand, int mode, double &av, double &sigma)
{
	//mode=0 - reading Random_file in continous region to get av,sigma
	//mode=1 - reading Random_file in entire range to make processed randoms
	
	if(mode==1 and Fexist(output)){return;}
	
    double R,RA_center,DEC_center;//for circle
    double dra=1.0*(ramax-ramin)/npix_edge, dsindec=(sin(decmax*deg2rad)-sin(decmin*deg2rad))/npix_edge;//pix size in ra, sin(pixsize) in dec
    double dist_tocenter[4]; // distance of border points to circle center
    int nnpix=pow(npix_edge,2); //all pixels
    double **pix_border=new double *[4];for(int i=0;i<4;++i){pix_border[i]=new double[nnpix];} //pix:ramin,ramax,decmin,decmax,
    
    R=Rmin*pow(1.0*qR,nnR); //size of circles on sphere [deg]; radius of sphere=1
	double S1=0.,S2=0.;

	int CR=Xcn.size();
	cout<<"Mode="<<mode<<"/1, circles: "<<CR<<endl;

    int galpix,found_thiscircle; //objects found inside this pixel, circle

    double Dmax,alpha_max,sinthetazero,costhetazero;
    int pixra[2],pixdec[2]; //pixel numbers ranges for each circle
	int pixra_taken,pixdec_taken,nnpix_taken; //analyzed pixels

	Writelog("Preparing for CiC: "+output);

	Set_pixborders(pix_border,nnpix,dra,dsindec,npix_edge);

	ofstream out(output.c_str());
	
	
    for(int j=0;j<CR;++j) //for every circle
    {
        if((j+1)%1000==0){cout<<"\rRadius_probe: "<<nnR+1<<"/"<<nR<<" circle: "<<j+1<<"/"<<CR;}
        found_thiscircle=0;

		RA_center=Xcn[j]; //avoiding the same selections for different scales not necessary in this case
		DEC_center=Ycn[j];

        sinthetazero=sin(.5*M_PI-DEC_center*deg2rad);
		costhetazero=cos(.5*M_PI-DEC_center*deg2rad);
		alpha_max=acos(pow(1.0-pow(costhetazero,2)-pow(sin(R*deg2rad),2),0.5)/sinthetazero)*rad2deg; //maximum right-ascension difference for radius R
		pixra[0]=max(0.,floor((RA_center-ramin)/dra)-ceil(alpha_max/dra));
		pixra[1]=min(1.*npix_edge-1.,floor((RA_center-ramin)/dra)+ceil(alpha_max/dra));
		pixdec[0]=max(0.,floor((sin((DEC_center-R)*deg2rad)-sin(decmin*deg2rad))/dsindec));
		pixdec[1]=min(1.*npix_edge-1.,floor((sin((DEC_center+R)*deg2rad)-sin(decmin*deg2rad))/dsindec));

		pixra_taken=pixra[0];
		pixdec_taken=pixdec[0];

        while(true) //loop over relevant pixels
		{
			nnpix_taken=pixdec_taken*npix_edge+pixra_taken;
			galpix=pix[nnpix_taken].size();

            dist_tocenter[0]=Dist(RA_center,DEC_center,pix_border[0][nnpix_taken],pix_border[2][nnpix_taken],"spherical");
            dist_tocenter[1]=Dist(RA_center,DEC_center,pix_border[0][nnpix_taken],pix_border[3][nnpix_taken],"spherical");
            dist_tocenter[2]=Dist(RA_center,DEC_center,pix_border[1][nnpix_taken],pix_border[2][nnpix_taken],"spherical");
            dist_tocenter[3]=Dist(RA_center,DEC_center,pix_border[1][nnpix_taken],pix_border[3][nnpix_taken],"spherical");

            Dmax=max(dist_tocenter[0],max(dist_tocenter[1],max(dist_tocenter[2],dist_tocenter[3])));

            if(Dmax<=R and pix[nnpix_taken][0][0]!=-1000000){found_thiscircle+=galpix;} //entire pixel inside the circle
            else
            {
                for(int m=0;m<galpix;++m) //objects inside this pixel (pixel crossing this circle partially)
                {
                    if(Dist(RA_center,DEC_center,pix[nnpix_taken][m][0],pix[nnpix_taken][m][1],"spherical")<=R){++found_thiscircle;}
                }
            }

			pixra_taken+=1;
			if(pixra_taken>pixra[1]){pixdec_taken+=1;pixra_taken=pixra[0];}
			if(pixdec_taken>pixdec[1]){break;} //out of range, everything analyzed
		}
		
		
		if(mode==0)
		{
			S1+=found_thiscircle;
			S2+=found_thiscircle*found_thiscircle;
		}
		
		if(mode==1)
		{
			if(found_thiscircle>=av-sigma*Random_sigmac) //counts not diverging by more than chosen nsigma from average
			{
				out<<Xcn[j]<<" "<<Ycn[j]<<endl;
			}
		}
		

    }
	out.close();
	
	
	if(mode==0)
	{
		av=S1/CR;
		sigma=sqrt(S2/(CR-1.) - (S1*S1)/((CR-1.)*CR));
		remove(output.c_str());//nothing there
		cout<<"Randoms for nnR="<<nnR<<": <N>="<<av<<", sigma="<<sigma<<endl;
	}
	cout<<endl;
    cout<<"Randoms CiC mode "<<mode<<"/1 done"<<endl;

	Writelog("Randoms CiC mode "+conv(mode)+"/1 done");

	
	return;
}

