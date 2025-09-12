#include "Randoms.h"




//creating temporary randoms in continuous field
//angular: [Random_fullr: ramin,ramax,decmin,decmax]
//LC_ellipses: [Random_fullr: Rmin,Rmax,ramin,ramax,decmin,decmax]
void Continuous_randoms(string frand_temp, int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn)
{
	double Rrand,rarand,decrand; //R,ra,dec of random position
	double R=0,axa,axb; //angular radius (VERSION=angular) and ellipses axes (VERSION=LC_ellipses)
	int C=0; //number of random points
	
	
	if(VERSION=="angular")
	{
		R=Rmin*pow(1.0*qR,nnR);
		C=kappa*nreals*Area_spherical(Random_fullr[0],Random_fullr[1],Random_fullr[2],Random_fullr[3])/(2.0*M_PI*(1.0-cos(R*deg2rad)));
	}
	if(VERSION=="LC_ellipses")
	{
		Find_ellipses_axes(nnR,axa,axb); //finding axes [axa,axb] for given nnR (=nnallax here)
		//kappa*nreals * 4/3 pi (Rmax^3 - Rmin^3) * (Area_spherical(fullr)/4pi) /ellipse_volume; ellipse_vol=4/3 pi axa* axb^2
		C=kappa*nreals*(Random_fullr[1]*Random_fullr[1]*Random_fullr[1] - Random_fullr[0]*Random_fullr[0]*Random_fullr[0])*\
		Area_spherical(Random_fullr[2],Random_fullr[3],Random_fullr[4],Random_fullr[5])/(4.*M_PI*axa*axb*axb);
	}
	
	Xcn.resize(C);
	Ycn.resize(C);
	Zcn.resize(C);
	
	if(Fexist(frand_temp)==0)
	{
		ofstream frandt(frand_temp.c_str());
		
		if(VERSION=="angular")
		{
			for(int i=0;i<C;++i)
			{
				Random_position<double>({Random_fullr[0],Random_fullr[1],Random_fullr[2],Random_fullr[3],R},rarand,decrand,"sky_border");
				frandt<<rarand<<" "<<decrand<<endl;
				Xcn[i]=rarand;
				Ycn[i]=decrand;
			}
		}
		
		if(VERSION=="LC_ellipses")
		{
			//for random distance:
			double powr=1./3.;
			double minnr3=(Random_fullr[0]+1.001*axa)*(Random_fullr[0]+1.001*axa)*(Random_fullr[0]+1.001*axa); //min_r^3 + avoiding tan(pi/2)
			double maxxr3=(Random_fullr[1]-axa)*(Random_fullr[1]-axa)*(Random_fullr[1]-axa); //max_r^3
			
			double alpha; //angular minimal distance from borders at chosen Rrand
			double alpha_tol=0.25*Dist(Random_fullr[2],Random_fullr[4],Random_fullr[3],Random_fullr[5],"astro"); //upper tolerance
			
			for(int i=0;i<C;++i)
			{
				while(true)
				{
					Rrand=pow(Random(minnr3,maxxr3),powr); //radial position (to have PDF \propto r^2 - constant 3D density)
					
					//angular part:
					alpha=atan(axb/sqrt(Rrand*Rrand - axa*axa));
					
					//setting angular random positions:
					Random_position<double>({Random_fullr[2],Random_fullr[3],Random_fullr[4],Random_fullr[5],alpha},rarand,decrand,"sky_border");
					
					if(alpha<alpha_tol){break;} //ellipsoid can be within [ramin,ramax,decmin,decmax] cone for that Rrand
				}
				
				
				Xcn[i]=Rrand*cos(rarand*deg2rad)*cos(decrand*deg2rad); //back to cartesian
				Ycn[i]=Rrand*sin(rarand*deg2rad)*cos(decrand*deg2rad);
				Zcn[i]=Rrand*sin(decrand*deg2rad);
				
				frandt<<Xcn[i]<<" "<<Ycn[i]<<" "<<Zcn[i]<<endl;
				
			}
		}
		

		frandt.close();
	}
	else //file already exist
	{
		Xcn.clear();Ycn.clear();Zcn.clear();
		if(VERSION=="angular"){Fread<double>(frand_temp,{&Xcn,&Ycn},{0,1});}
		if(VERSION=="LC_ellipses"){Fread<double>(frand_temp,{&Xcn,&Ycn,&Zcn},{0,1,2});}
		
	}
	
	return;
}






//processing randoms for one scale
void Process_randoms_onerad(int rank, int comm_size,int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn,
vector<double> &Nav_circ, vector<double> &Nav_pix, int &nrand)
{
	//temporary randoms within continuous field defined by Random_fullr
	string ftemp="Randoms/temp_"+Random_str(10);
	string frand_temp="Randoms/continuous_temp_"+conv(nnR)+".txt";
	string frand="Randoms/Random_proc_"+conv(nnR)+".txt";
	if(Fexist(frand)==1){return;}
	
	vector<double> x,y,z; //for data - here Random_file
	double av,sigma;
	int npix;
	vector<vector<vector<double> > > pix; //for pixelization
	
	Continuous_randoms(frand_temp,nnR,Xcn,Ycn,Zcn); //creating randoms on continuous field
	Read_randoms(rank,comm_size,x,y,z,"",nrand,nnR,Nav_circ,Nav_pix,-1); //reading Random_file into x,y,(z) depending on VERSION
	
	
	if(VERSION=="angular")
	{
		npix=Npix(Nav_circ,Nav_pix,x.size(),Rmin*pow(1.0*qR,nnR),0,0); //optimal pixelization
		Pixelize<double>({ramin,ramax,decmin,decmax},{npix,npix},{x,y},pix,"2D_sky");
		
		//get average and sigma from continuous field:
		CIC_angular_rnd(pix,ftemp,npix,nnR,Xcn,Ycn,Zcn,nrand,0,av,sigma); //counting - mode=0
			
		Xcn.clear();
		Ycn.clear();
		
		//selecting data from Random_file where counts> average - Random_sigmac (CIC with mode=3): [data already pixelized]
		CIC_angular_rnd(pix,frand,npix,nnR,x,y,z,nrand,1,av,sigma); //counting - mode=1;  position of circles the same as Random_file

	}
	
	if(VERSION=="LC_ellipses")
	{
		double axa,axb;
		Find_ellipses_axes(nnR,axa,axb);
		npix=Npix(Nav_circ,Nav_pix,x.size(),0,axa,axb); //optimal pixelization
		
		Pixelize<double>({-1.*DCMAX,DCMAX,-1.*DCMAX,DCMAX,-1.*DCMAX,DCMAX},{npix,npix,npix},{x,y,z},pix,"3D");
		
		//get average and sigma from continuous field:
		CIC_LC_ellipses_rnd(pix,ftemp,npix,nnR,Xcn,Ycn,Zcn,nrand,0,av,sigma); //counting - mode=0
			
		Xcn.clear();
		Ycn.clear();
		Zcn.clear();
		
		//selecting data from Random_file where counts> average - Random_sigmac (CIC with mode=3): [data already pixelized]
		CIC_LC_ellipses_rnd(pix,frand,npix,nnR,x,y,z,nrand,1,av,sigma); //counting - mode=1;  position of circles the same as Random_file
		
	}
	
	remove(frand_temp.c_str()); //removing
	
	
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
	
	int njobs=VERSION=="angular"?nR:n_allax; //number of all jobs
	
	//temporary solution------------------------------------------------
	for(int i=0;i<njobs;++i) //for all scales
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
	
	if(mode==1 and Fexist(output)==1){return;}
	
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

            dist_tocenter[0]=Dist(RA_center,DEC_center,pix_border[0][nnpix_taken],pix_border[2][nnpix_taken],"astro");
            dist_tocenter[1]=Dist(RA_center,DEC_center,pix_border[0][nnpix_taken],pix_border[3][nnpix_taken],"astro");
            dist_tocenter[2]=Dist(RA_center,DEC_center,pix_border[1][nnpix_taken],pix_border[2][nnpix_taken],"astro");
            dist_tocenter[3]=Dist(RA_center,DEC_center,pix_border[1][nnpix_taken],pix_border[3][nnpix_taken],"astro");

            Dmax=max(dist_tocenter[0],max(dist_tocenter[1],max(dist_tocenter[2],dist_tocenter[3])));

            if(Dmax<=R and pix[nnpix_taken][0][0]!=-1000000){found_thiscircle+=galpix;} //entire pixel inside the circle
            else
            {
                for(int m=0;m<galpix;++m) //objects inside this pixel (pixel crossing this circle partially)
                {
                    if(Dist(RA_center,DEC_center,pix[nnpix_taken][m][0],pix[nnpix_taken][m][1],"astro")<=R){++found_thiscircle;}
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








//CiC for LC_ellipses, but version especially for random processing
void CIC_LC_ellipses_rnd(vector<vector<vector<double> > > &pix, string output,
int npix_edge, int nnallax, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand, int mode, double &av, double &sigma)
{
	//mode=0 - reading Random_file in continous region to get av,sigma
	//mode=1 - reading Random_file in entire range to make processed randoms
	
	if(mode==1 and Fexist(output)==1){return;}
	
	double *CENTER=new double[3]; //positions of ellipse centers
    double dpix=2.0*DCMAX/npix_edge;//pix edge size
    double *dist_tocenter=new double[8]; // distance of border points to circle center
    double axa,axb,minax,maxax;
	double F1[3],F2[3],CFOCI,alpha; //ellipse foci and CFOCI(half foci distance),alpha - direction rescaling
	double *direction=new double[3]; //for finding focal points
	double *w=new double[3]; //w: alpha*dir+w=galaxy_position
	double S1=0.,S2=0.;
	
	int CR=Xcn.size();
	cout<<"Mode="<<mode<<"/1, ellipses: "<<CR<<endl;
	
	Find_ellipses_axes(nnallax,axa,axb);
	minax=min(axa,axb);
	maxax=max(axa,axb);
	if(axa>=axb){CFOCI=pow(pow(axa,2)-pow(axb,2),0.5);}
	else{CFOCI=pow(pow(axb,2)-pow(axa,2),0.5);}

    int galpix,found_thiscircle; //objects found inside this pixel, circle
	double Dmax;
	int pixx[2],pixy[2],pixz[2]; //pixel numbers ranges for each circle
	int pixx_taken,pixy_taken,pixz_taken,nnpix_taken; //analyzed pixels
	
	Writelog("Preparing for CiC: "+output);
	
	ofstream out(output.c_str());

    for(int j=0;j<CR;++j) //for every ellipse
    {
		if((j+1)%100==0){cout<<"\rProbe: "<<nnallax+1<<"/"<<n_allax<<" ellipsoid: "<<j+1<<"/"<<CR;}
        found_thiscircle=0;
		
		CENTER[0]=Xcn[j]; //avoiding the same selections for different scales not necessary in this case
		CENTER[1]=Ycn[j];
		CENTER[2]=Zcn[j];
		
		for(int k=0;k<3;++k){direction[k]=CENTER[k];}
		Normalize(direction);
		
		//setting focal points: [the axa<axb case is for each object inside m loop]
		if(axa>=axb)
		{
			for(int k=0;k<3;++k)
			{
				F1[k]=CENTER[k] - direction[k]*CFOCI;
				F2[k]=CENTER[k] + direction[k]*CFOCI;
			}
		}

		
		//relevant pixel ranges [for this catalog topology different than previously]
		pixx[0]=max(0.,floor((CENTER[0]+DCMAX)/dpix)-ceil(maxax/dpix)); 
		pixx[1]=min(1.*npix_edge-1.,floor((CENTER[0]+DCMAX)/dpix)+ceil(maxax/dpix));
		pixy[0]=max(0.,floor((CENTER[1]+DCMAX)/dpix)-ceil(maxax/dpix));
		pixy[1]=min(1.*npix_edge-1.,floor((CENTER[1]+DCMAX)/dpix)+ceil(maxax/dpix));
		pixz[0]=max(0.,floor((CENTER[2]+DCMAX)/dpix)-ceil(maxax/dpix));
		pixz[1]=min(1.*npix_edge-1.,floor((CENTER[2]+DCMAX)/dpix)+ceil(maxax/dpix));
		
		pixx_taken=pixx[0];
		pixy_taken=pixy[0];
		pixz_taken=pixz[0];

        while(true) //loop over relevant pixels
		{
			nnpix_taken=pixz_taken*pow(npix_edge,2) +pixy_taken*npix_edge +pixx_taken;
			galpix=pix[nnpix_taken].size();
			Dmax=Set_3Dvertex_distances(dist_tocenter,CENTER[0],CENTER[1],CENTER[2],pixx_taken,pixy_taken,pixz_taken,dpix);

            if(Dmax<=minax and pix[nnpix_taken][0][0]!=-1000000){found_thiscircle+=galpix;} //entire pixel inside the ellipse
            else
            {
                for(int m=0;m<galpix;++m) //objects inside this pixel (pixel crossing this circle partially)
                {
					//from ellipse definition: sum of distances to focii is =2a on ellipse
					if(axa<axb) //sophisticated case, must compute focii for each object (on ring)
					{
						alpha=(pix[nnpix_taken][m][0]*direction[0]+pix[nnpix_taken][m][1]*direction[1]+pix[nnpix_taken][m][2]*direction[2]); //direction already normalized
						for(int nn=0;nn<3;++nn){w[nn]=pix[nnpix_taken][m][nn]-alpha*direction[nn];} //perpendicular vec connecting alpha*dir and galaxy
						Normalize(w);
						
						for(int nn=0;nn<3;++nn) //setting focal points for sophisticated case
						{
							F1[nn]=CENTER[nn] - CFOCI*w[nn];
							F2[nn]=CENTER[nn] + CFOCI*w[nn];
						}
					}
					
                    if(Dist(pix[nnpix_taken][m][0],pix[nnpix_taken][m][1],pix[nnpix_taken][m][2],F1[0],F1[1],F1[2])+
					   Dist(pix[nnpix_taken][m][0],pix[nnpix_taken][m][1],pix[nnpix_taken][m][2],F2[0],F2[1],F2[2])<=2.*maxax)
					{
						++found_thiscircle;
					}
                }
            }

			pixx_taken+=1;
			if(pixx_taken>pixx[1])
			{
				pixy_taken+=1;
				if(pixy_taken>pixy[1])
				{
					pixz_taken+=1;
					pixy_taken=pixy[0];
				}
				pixx_taken=pixx[0];
			}
			if(pixz_taken>pixz[1]){break;} //out of range, everything analyzed
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
				out<<Xcn[j]<<" "<<Ycn[j]<<" "<<Zcn[j]<<endl;
			}
		}
		
		
    }
	out.close();
	
	if(mode==0)
	{
		av=S1/CR;
		sigma=sqrt(S2/(CR-1.) - (S1*S1)/((CR-1.)*CR));
		remove(output.c_str()); //nothing there
		cout<<"Randoms for nnallax="<<nnallax<<": <N>="<<av<<", sigma="<<sigma<<endl;
	}
	cout<<endl;
    cout<<"Randoms CiC mode "<<mode<<"/1 done"<<endl;

	Writelog("Randoms CiC mode "+conv(mode)+"/1 done");
	
	return;
}

