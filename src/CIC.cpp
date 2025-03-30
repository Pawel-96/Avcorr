#include "CIC.h"



int CRR(string ffound,string ffoundhst,vector<double> &Xcn, string model, double val1, double val2, double val3, double val4) //number of circles/ellipses drawn
{	//val1/2/3/4: for angular and BOX: R/0/0/0,  for BOX_ellipses: axa/axb/0/0, for LC_ellipses: axa/axb/DCMIN/DCMAX
	int C=0,CR;
	if(VERSION=="angular"){C=kappa*(Areaf/str2deg2)/(2.0*M_PI*(1.0-cos(val1*deg2rad)));}
	if(VERSION=="BOX"){C=kappa*pow(boxsize/val1,3)/(4.0/3.0 *M_PI);}
	if(VERSION=="BOX_ellipses"){ C=kappa*pow(boxsize,3)/(val1*pow(val2,2)*4./3. *M_PI);}
	if(VERSION=="LC_ellipses"){C=kappa*(pow(val4,3)-pow(val3,3))/(val1*pow(val2,2));}
	
	CR=min(max(C*nreals,Cmin),Cmax);
	
	if(Fexist(ffound)==1){CR=min(CR,Fnlines(ffound));} //3 tree steps: random_provided [0/1],random_file=="*"[y/n],ffound exists?[y/n]
	if(Fexist(ffoundhst)==1)
	{
		vector<double> nhist;
		Fread<double>(ffoundhst,{&nhist},{1});
		int n=Sum(nhist); 
		CR=n;//min(CR,n);
	}
	if(Random_provided==1 and Fexist(ffound)==0 and Fexist(ffoundhst)==0)
	{
		if(Random_file=="*"){CR=Fnlines("Randoms/Random_"+model+EXT);}//min(CR,(long int)(Lines("Randoms/Random_"+model+".txt")));}
		else{CR=Xcn.size();}//min(CR,(long int)(Xcn.size()));}
	}
	
	return CR;
}




void Set_pixborders(double **pix_border, int npix, double dra, double dsindec, int npix_edge) //borders of pixels: ramin/max, decmin/max
{
    for(int k=0;k<npix;++k)
    {
        pix_border[0][k]=ramin+ dra*(k%npix_edge);//pix_ramin
        pix_border[1][k]=pix_border[0][k] +dra;//pix_ramax
        pix_border[2][k]=asin(sin(decmin*deg2rad)+(1.0*(int)(k/npix_edge))*dsindec)*rad2deg;//pix_decmin
        pix_border[3][k]=asin(sin(decmin*deg2rad)+(1.0 + 1.0*(int)(k/npix_edge))*dsindec)*rad2deg;//pix_decmax
    }
    return;
}






double Set_3Dvertex_distances(double *dist_tocenter, double X, double Y, double Z, int pixx_taken, int pixy_taken, int pixz_taken, double dpix)
{
	double Dmax;
	if(VERSION=="BOX" or VERSION=="BOX_ellipses")
	{
		dist_tocenter[0]=Dist(X,Y,Z,  pixx_taken*dpix,        pixy_taken*dpix,        pixz_taken*dpix);
		dist_tocenter[1]=Dist(X,Y,Z,  pixx_taken*dpix,        pixy_taken*dpix,        (pixz_taken+1.0)*dpix);
		dist_tocenter[2]=Dist(X,Y,Z,  pixx_taken*dpix,        (pixy_taken+1.0)*dpix,  pixz_taken*dpix);
		dist_tocenter[3]=Dist(X,Y,Z,  (pixx_taken+1.0)*dpix,  pixy_taken*dpix,        pixz_taken*dpix);
		dist_tocenter[4]=Dist(X,Y,Z,  (pixx_taken+1.0)*dpix,  (pixy_taken+1.0)*dpix,  pixz_taken*dpix);
		dist_tocenter[5]=Dist(X,Y,Z,  (pixx_taken+1.0)*dpix,  pixy_taken*dpix,        (pixz_taken+1.0)*dpix);
		dist_tocenter[6]=Dist(X,Y,Z,  pixx_taken*dpix,        (pixy_taken+1.0)*dpix,  (pixz_taken+1.0)*dpix);
		dist_tocenter[7]=Dist(X,Y,Z,  (pixx_taken+1.0)*dpix,  (pixy_taken+1.0)*dpix,  (pixz_taken+1.0)*dpix);
	}
	
	if(VERSION=="LC_ellipses")
	{
		dist_tocenter[0]=Dist(X,Y,Z,  -1.*DCMAX+pixx_taken*dpix,        -1.*DCMAX+pixy_taken*dpix,        -1.*DCMAX+pixz_taken*dpix);
		dist_tocenter[1]=Dist(X,Y,Z,  -1.*DCMAX+pixx_taken*dpix,        -1.*DCMAX+pixy_taken*dpix,        -1.*DCMAX+(pixz_taken+1.0)*dpix);
		dist_tocenter[2]=Dist(X,Y,Z, -1.*DCMAX+pixx_taken*dpix,        -1.*DCMAX+(pixy_taken+1.0)*dpix,  -1.*DCMAX+pixz_taken*dpix);
		dist_tocenter[3]=Dist(X,Y,Z, -1.*DCMAX+(pixx_taken+1.0)*dpix,  -1.*DCMAX+pixy_taken*dpix,        -1.*DCMAX+pixz_taken*dpix);
		dist_tocenter[4]=Dist(X,Y,Z,  -1.*DCMAX+(pixx_taken+1.0)*dpix,  -1.*DCMAX+(pixy_taken+1.0)*dpix,  -1.*DCMAX+pixz_taken*dpix);
		dist_tocenter[5]=Dist(X,Y,Z,  -1.*DCMAX+(pixx_taken+1.0)*dpix,  -1.*DCMAX+pixy_taken*dpix,        -1.*DCMAX+(pixz_taken+1.0)*dpix);
		dist_tocenter[6]=Dist(X,Y,Z,  -1.*DCMAX+pixx_taken*dpix,        -1.*DCMAX+(pixy_taken+1.0)*dpix,  -1.*DCMAX+(pixz_taken+1.0)*dpix);
		dist_tocenter[7]=Dist(X,Y,Z,  -1.*DCMAX+(pixx_taken+1.0)*dpix,  -1.*DCMAX+(pixy_taken+1.0)*dpix,  -1.*DCMAX+(pixz_taken+1.0)*dpix);	
	}
	
	Dmax=Minmax<double>(dist_tocenter,8,"max");
	return Dmax;
}





void Find_ellipses_axes(int nnallax, double &axa, double &axb) //finding axes values for this nnallax
{
	int nnaxa,nnaxb;
	nnaxa=nnallax/naxb; //number of semi-major axis grid point
	nnaxb=nnallax-nnaxa*naxb; //number of semi-minor axis grid point
	axa=axamin*pow(1.*qaxa,nnaxa); //current semi-major axis value
	axb=axbmin*pow(1.*qaxb,nnaxb); //current semi-minor axis value
	return;
}






//adding count to count array for data storage
void Add_count(vector<int> &count, int &countmin, int &countmax, int found_thiscircle)
{
	if(found_thiscircle>=countmin and found_thiscircle<=countmax)
		{
			count[found_thiscircle-countmin]+=1;
			return;
		}
		
		if(countmin==-1) //first count
		{
			count.push_back(1);
			countmin=found_thiscircle;
			countmax=found_thiscircle;
			return;
		}
		
		if(found_thiscircle<countmin) //count smaller than actual array range
		{
			count.insert(count.begin(),(int)(countmin-found_thiscircle),0); //adding zeros to count array (on the beginnning)
			count[0]+=1;
			countmin=found_thiscircle;
			return;
		}
		
		if(found_thiscircle>countmax) //count bigger than actual array range
		{
			count.insert(count.end(),(int)(found_thiscircle-countmax),0); //adding zeros to count array (on ending)
			count[found_thiscircle-countmin]+=1;
			countmax=found_thiscircle;
			return;
		}
}







//removing counts from histogram where count[i]=0
void Shift_zeros(vector<int> &nc, vector<int> &count)
{
	vector<int> nc_new,count_new;
	int nnc=nc.size();
	for(int i=0;i<nnc;++i)
	{
		if(count[i]>0)
		{
			nc_new.push_back(nc[i]);
			count_new.push_back(count[i]);
		}
	}
	nc.clear();
	count.clear();
	
	nc=nc_new;
	count=count_new;
	
	return;
}







//checking whether count exist
int Check_count_existence(string ffound, string ffoundhst, string output, vector<int> &nc, vector<int> &count,  double R)
{
	int exst=0;
	if(Fexist(ffound)==1 or Fexist(ffoundhst)==1) //f[] collected earlier
    {
		Writelog("CiC done -> collecting counts [exists] <- "+ffound);
        if(Fexist(ffound)==1)
		{
			vector<int> found;
			Fread<int>(ffound,{&found},{0});
			Hist(found,nc,count,"discrete"); //creating nc,count histogram from read data
		}
		if(Fexist(ffoundhst)==1){Fread<int>(ffoundhst,{&nc,&count},{0,1});}
        Calculate_moments(output,nc,count,R);
        exst=1;
    }
	return exst;
}






void CIC_angular(vector<vector<vector<double> > > &pix, string output, string model, int npix_edge, int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn) //collecting numbers of objects  ->calculating moments for nnR-th radius
{
    double R,RA_center,DEC_center;//for circle
    double dra=1.0*(ramax-ramin)/npix_edge, dsindec=(sin(decmax*deg2rad)-sin(decmin*deg2rad))/npix_edge;//pix size in ra, sin(pixsize) in dec
    double dist_tocenter[4]; // distance of border points to circle center
    int nnpix=pow(npix_edge,2); //all pixels
    double **pix_border=new double *[4];for(int i=0;i<4;++i){pix_border[i]=new double[nnpix];} //pix:ramin,ramax,decmin,decmax,
    
    R=Rmin*pow(1.0*qR,nnR); //size of circles on sphere [deg]; radius of sphere=1

	string ffound=Replace_string(output,".txt","_found.fnd"); //file with counts
	string ffoundhst=Replace_string(ffound,".fnd",".fndhst"); //count in compact format
	int CR=CRR(ffound,ffoundhst,Xcn,model,R,0.,0.,0.);
	cout<<"Circles: "<<CR<<endl;

    vector<int> nc, count; //storing counts histogram
	int countmin=-1,countmax=-1; //minimal count (for adding to histogram)
    int galpix,found_thiscircle; //objects found inside this pixel, circle

    double Dmax,alpha_max,sinthetazero,costhetazero;
    int pixra[2],pixdec[2]; //pixel numbers ranges for each circle
	int pixra_taken,pixdec_taken,nnpix_taken; //analyzed pixels

	Writelog("Preparing for CiC: "+output);
	int exst=Check_count_existence(ffound,ffoundhst,output,nc,count,R);
    if(exst==1){return;}

	Set_pixborders(pix_border,nnpix,dra,dsindec,npix_edge);

	if(Random_provided==1 and Random_file=="*")
	{
			Xcn.clear();Ycn.clear();
			Fread<double>("Randoms/Random_"+model+EXT,{&Xcn,&Ycn},{0,1});
	}

    for(int j=0;j<CR;++j) //for every circle
    {
        if((j+1)%100000==0){cout<<"\rRadius_probe: "<<nnR+1<<"/"<<nR<<" circle: "<<j+1<<"/"<<CR;}
        found_thiscircle=0;
		
		if(Random_provided==1)
		{
			RA_center=Xcn[j];
			DEC_center=Ycn[j];		
		}
        else{Random_position({ramin,ramax,decmin,decmax,R},RA_center,DEC_center,"sky_border");}

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
		Add_count(count,countmin,countmax,found_thiscircle); //adding found value to histogram (better data storage)

    }
	cout<<endl;
    cout<<"CiC done, then -> "<<output<<endl;

	Writelog("CiC done -> saving counts -> "+ffoundhst);
	Arange(nc,countmin,countmax,1);
	Shift_zeros(nc,count); //removing zeros from hist to save space
	Fwrite<int>(ffoundhst,{&nc,&count});
    Calculate_moments(output,nc,count,R);
	
	return;
}





//periodic pixel (assumption: pix_taken not exceeding [-npix_1D,2npix_1D] range)
int Periodic_pixel(int pix_taken, int npix_1D)
{
	if(pix_taken>=0 and pix_taken<npix_1D){return pix_taken;}
	
	if(pix_taken<0){return pix_taken+npix_1D;}
	else if(pix_taken>=npix_1D){return pix_taken-npix_1D;}
	
	return -1;
}




//distance in 3D periodic box
double Dist_periodic(double x1, double y1, double z1, double x2, double y2, double z2, double Lbox)
{
	double Lhalf=.5*Lbox;
	double dx=fabs(x1-x2),dy=fabs(y1-y2),dz=fabs(z1-z2);
	if(dx>Lhalf){dx=Lbox-dx;}
	if(dy>Lhalf){dy=Lbox-dy;}
	if(dz>Lhalf){dz=Lbox-dz;}
	
	return sqrt(dx*dx + dy*dy + dz*dz);
}









void CIC_BOX(vector<vector<vector<double> > > &pix, string output, string model, int npix_edge, int nnR,vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn) //collecting numbers of objects  ->calculating moments for nnR-th radius
{
    double R,X_center,Y_center,Z_center;//for sphere
    double dpix=1.0*boxsize/npix_edge;//pix edge size
    double *dist_tocenter=new double[8]; // distance of border points to circle center
    
    R=Rmin*pow(1.0*qR,nnR); //size of sphere radius

	string ffound=Replace_string(output,".txt","_found.fnd"); //file with found[] array
	string ffoundhst=Replace_string(ffound,".fnd",".fndhst"); //count in compact format
	int CR=CRR(ffound,ffoundhst,Xcn,model,R,0.,0.,0.);
	cout<<"Spheres: "<<CR<<endl;

    vector<int> nc, count; //storing counts histogram
	int countmin=-1,countmax=-1; //minimal count (for adding to histogram)
    int galpix,found_thiscircle; //objects found inside this pixel, circle
	
	double Dmax;
	int pixx[2],pixy[2],pixz[2]; //pixel numbers ranges for each circle
	int pixx_taken,pixy_taken,pixz_taken,nnpix_taken; //analyzed pixels
	
	Writelog("Preparing for CiC: "+output);
    int exst=Check_count_existence(ffound,ffoundhst,output,nc,count,R);
    if(exst==1){return;}
	
	if(Random_provided==1 and Random_file=="*")
	{
			Xcn.clear();Ycn.clear();Zcn.clear();
			Fread<double>("Randoms/Random_"+model+EXT,{&Xcn,&Ycn,&Zcn},{0,1,2});
	}

    for(int j=0;j<CR;++j) //for every circle
    {
        if((j+1)%100000==0){cout<<"\rRadius_probe: "<<nnR+1<<"/"<<nR<<" sphere: "<<j+1<<"/"<<CR;}
        found_thiscircle=0;
		if(Random_provided==1)
		{
			X_center=Xcn[j];
			Y_center=Ycn[j];
			Z_center=Zcn[j];			
		}
		else {
			X_center=Random(0.,boxsize);
			Y_center=Random(0.,boxsize);
			Z_center=Random(0.,boxsize);			
		}

		//extreme relevant pixels (pixx/y/z[0/1]); might be out of the box; then periodicity
		pixx[0]=floor(X_center/dpix)-ceil(R/dpix);
		pixx[1]=floor(X_center/dpix)+ceil(R/dpix);
		pixy[0]=floor(Y_center/dpix)-ceil(R/dpix);
		pixy[1]=floor(Y_center/dpix)+ceil(R/dpix);
		pixz[0]=floor(Z_center/dpix)-ceil(R/dpix);
		pixz[1]=floor(Z_center/dpix)+ceil(R/dpix);
		
		pixx_taken=pixx[0];
		pixy_taken=pixy[0];
		pixz_taken=pixz[0];

        while(true) //loop over relevant pixels
		{
			nnpix_taken=Periodic_pixel(pixz_taken,npix_edge)*npix_edge*npix_edge+
							Periodic_pixel(pixy_taken,npix_edge)*npix_edge+
							Periodic_pixel(pixx_taken,npix_edge); //converting back to 1D listing; assuming periodicity
			galpix=pix[nnpix_taken].size();
			Dmax=Set_3Dvertex_distances(dist_tocenter,X_center,Y_center,Z_center,pixx_taken,pixy_taken,pixz_taken,dpix);
			
            if(Dmax<=R and pix[nnpix_taken][0][0]!=-1000000){found_thiscircle+=galpix;} //entire pixel inside the circle
            else
            {
                for(int m=0;m<galpix;++m) //objects inside this pixel (pixel crossing this circle partially)
                {
                    if(Dist_periodic(X_center,Y_center,Z_center,
					   pix[nnpix_taken][m][0],pix[nnpix_taken][m][1],pix[nnpix_taken][m][2],boxsize)<=R){++found_thiscircle;}
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
		Add_count(count,countmin,countmax,found_thiscircle); //adding found value to histogram (better data storage)
		
    }
	cout<<endl;
    cout<<"CiC done, then -> "<<output<<endl;
	
	Writelog("CiC done -> saving counts -> "+ffoundhst);
	Arange(nc,countmin,countmax,1);
	Shift_zeros(nc,count); //removing zeros from hist to save space
	Fwrite<int>(ffoundhst,{&nc,&count});
    Calculate_moments(output,nc,count,R);
	
	return;
}









void CIC_BOX_ellipses(vector<vector<vector<double> > > &pix, string output, string model, int npix_edge, int nnallax,vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn) //collecting numbers of objects  ->calculating moments for nnR-th radius
{
    double X_center,Y_center,Z_center;//for sphere
    double dpix=1.0*boxsize/npix_edge;//pix edge size
    double *dist_tocenter=new double[8]; // distance of border points to circle center
	double Lhalf=.5*boxsize;
    double axa,axb,minax; //axes and smaller axis
	double dx,dy,dz,inellips; //useful for ellipsoid equation: (dx/axa)^2+(dy/axb)^2+(dz/axb)^2<1
	//double CFOCI,frac,Dfoc; //CFOCI(half foci distance), frac (for rescaling), Dfoc -distance to focal ring x2 (axb>axa case)
    
	Find_ellipses_axes(nnallax,axa,axb);
	minax=min(axa,axb);
	//if(axa>=axb){CFOCI=pow(pow(axa,2)-pow(axb,2),0.5);}
	//else{CFOCI=pow(pow(axb,2)-pow(axa,2),0.5);}
	
    string ffound=Replace_string(output,".txt","_found.fnd"); //file with found[] array
	string ffoundhst=Replace_string(ffound,".fnd",".fndhst"); //count in compact format
	int CR=CRR(ffound,ffoundhst,Xcn,model,axa,axb,0.,0.);
	cout<<"Ellipses: "<<CR<<endl;

    vector<int> nc, count; //storing counts histogram
	int countmin=-1,countmax=-1; //minimal count (for adding to histogram)
    int galpix,found_thiscircle; //objects found inside this pixel, ellipsoid
	
	double Dmax;
	int pixx[2],pixy[2],pixz[2]; //pixel numbers ranges for each circle
	int pixx_taken,pixy_taken,pixz_taken,nnpix_taken; //analyzed pixels
	
	Writelog("Preparing for CiC: "+output);
    int exst=Check_count_existence(ffound,ffoundhst,output,nc,count,nnallax);
    if(exst==1){return;}
	
	if(Random_provided==1 and Random_file=="*")
	{
			Xcn.clear();Ycn.clear();Zcn.clear();
			Fread<double>("Randoms/Random_"+model+EXT,{&Xcn,&Ycn,&Zcn},{0,1,2});
	}

    for(int j=0;j<CR;++j) //for every ellipse
    {
        if((j+1)%1000==0){cout<<"\rRadius_probe: "<<nnallax+1<<"/"<<n_allax<<" ellipsoid: "<<j+1<<"/"<<CR;}
        found_thiscircle=0;
		if(Random_provided==1)
		{
			X_center=Xcn[j];
			Y_center=Ycn[j];
			Z_center=Zcn[j];			
		}
		else {
			X_center=Random(0.,boxsize);
			Y_center=Random(0.,boxsize);
			Z_center=Random(0.,boxsize);			
		}

		
		pixx[0]=floor(X_center/dpix)-ceil(axa/dpix);
		pixx[1]=floor(X_center/dpix)+ceil(axa/dpix);
		pixy[0]=floor(Y_center/dpix)-ceil(axb/dpix);
		pixy[1]=floor(Y_center/dpix)+ceil(axb/dpix);
		pixz[0]=floor(Z_center/dpix)-ceil(axb/dpix);
		pixz[1]=floor(Z_center/dpix)+ceil(axb/dpix);
		
		pixx_taken=pixx[0];
		pixy_taken=pixy[0];
		pixz_taken=pixz[0];

        while(true) //loop over relevant pixels
		{
			nnpix_taken=Periodic_pixel(pixz_taken,npix_edge)*npix_edge*npix_edge+
							Periodic_pixel(pixy_taken,npix_edge)*npix_edge+
							Periodic_pixel(pixx_taken,npix_edge); //converting back to 1D listing; assuming periodicity
			galpix=pix[nnpix_taken].size();
			Dmax=Set_3Dvertex_distances(dist_tocenter,X_center,Y_center,Z_center,pixx_taken,pixy_taken,pixz_taken,dpix);

            if(Dmax<=minax and pix[nnpix_taken][0][0]!=-1000000){found_thiscircle+=galpix;} //entire pixel inside the ellipsoid
            else
            {
                for(int m=0;m<galpix;++m) //objects inside this pixel (pixel crossing this circle partially)
                {
					dx=fabs(pix[nnpix_taken][m][0]-X_center);
					dy=fabs(pix[nnpix_taken][m][1]-Y_center);
					dz=fabs(pix[nnpix_taken][m][2]-Z_center);
					
					if(dx>Lhalf){dx=boxsize-dx;}
					if(dy>Lhalf){dy=boxsize-dy;}
					if(dz>Lhalf){dz=boxsize-dz;}
					
					inellips=dx*dx/(axa*axa) + (dy*dy + dz*dz)/(axb*axb); //from ellipse equation
					if(inellips<1.){++found_thiscircle;} //from ellipsoid definition    
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
        Add_count(count,countmin,countmax,found_thiscircle); //adding found value to histogram (better data storage)
    }
	cout<<endl;
    cout<<"CiC done, then -> "<<output<<endl;
	
	Writelog("CiC done -> saving counts -> "+ffoundhst);
	Arange(nc,countmin,countmax,1);
	Shift_zeros(nc,count); //removing zeros from hist to save space
	Fwrite<int>(ffoundhst,{&nc,&count});
    Calculate_moments(output,nc,count,nnallax);
	
	return;
}






void CIC_LC_ellipses(vector<vector<vector<double> > > &pix, string output, string model, int npix_edge, int nnallax,vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn) //collecting numbers of objects  ->calculating moments for nnR-th radius
{
	double *CENTER=new double[3]; //positions of ellipse centers
    double dpix=2.0*DCMAX/npix_edge;//pix edge size
    double *dist_tocenter=new double[8]; // distance of border points to circle center
    double axa,axb,minax,maxax;
	double F1[3],F2[3],CFOCI,alpha; //ellipse foci and CFOCI(half foci distance),alpha - direction rescaling
	double *direction=new double[3]; //for finding focal points
	double *w=new double[3]; //w: alpha*dir+w=galaxy_position
	
	Find_ellipses_axes(nnallax,axa,axb);
	minax=min(axa,axb);
	maxax=max(axa,axb);
	if(axa>=axb){CFOCI=pow(pow(axa,2)-pow(axb,2),0.5);}
	else{CFOCI=pow(pow(axb,2)-pow(axa,2),0.5);}
		
    string ffound=Replace_string(output,".txt","_found.fnd"); //file with found[] array
	string ffoundhst=Replace_string(ffound,".fnd",".fndhst"); //count in compact format
	int CR=CRR(ffound,ffoundhst,Xcn,model,axa,axb,DCMIN,DCMAX);
	cout<<"Ellipses: "<<CR<<endl;

    vector<int> nc, count; //storing counts histogram
	int countmin=-1,countmax=-1; //minimal count (for adding to histogram)
    int galpix,found_thiscircle; //objects found inside this pixel, circle
	double Dmax;
	int pixx[2],pixy[2],pixz[2]; //pixel numbers ranges for each circle
	int pixx_taken,pixy_taken,pixz_taken,nnpix_taken; //analyzed pixels
	
	Writelog("Preparing for CiC: "+output);
    int exst=Check_count_existence(ffound,ffoundhst,output,nc,count,nnallax);
    if(exst==1){return;}
	
	if(Random_provided==1 and Random_file=="*")
	{
			Xcn.clear();Ycn.clear();Zcn.clear();
			Fread<double>("Randoms/Random_"+model+EXT,{&Xcn,&Ycn,&Zcn},{0,1,2});
	}
	

    for(int j=0;j<CR;++j) //for every ellipse
    {
		if(Random_provided==1)
		{
			CENTER[0]=Xcn[j];
			CENTER[1]=Ycn[j];
			CENTER[2]=Zcn[j];			
		}
        else{Random_position({DCMIN+axa,DCMAX-axa},CENTER[0],CENTER[1],CENTER[2],"shell");}
		
		if((j+1)%100==0){cout<<"\rProbe: "<<nnallax+1<<"/"<<n_allax<<" ellipsoid: "<<j+1<<"/"<<CR;}
        found_thiscircle=0;
		
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
        Add_count(count,countmin,countmax,found_thiscircle); //adding found value to histogram (better data storage)
    }
	cout<<endl;
    cout<<"CiC done, then -> "<<output<<endl;
	
	Writelog("CiC done -> saving counts -> "+ffoundhst);
	Arange(nc,countmin,countmax,1);
	Shift_zeros(nc,count); //removing zeros from hist to save space
	Fwrite<int>(ffoundhst,{&nc,&count});
    Calculate_moments(output,nc,count,nnallax);
	
	return;
}