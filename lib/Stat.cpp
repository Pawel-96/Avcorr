#include "Stat.h"









string Random_str(int n) //random string of length n
{
	string signs="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
	string res;
	int nn=signs.size();
	for(int i=0;i<n;++i)
	{
		res+=signs[rand() %nn];
	}
	return res;
}





int Chose_probability(double p) //returns 1 if chosen with probability p, otherwise 0
{
	double Max=1000.0;
    double a=Random(0.0,Max);
    if(a<p*Max){return 1;}
    else return 0;
}







//finding which pixel is nnpix_{x,y,z}-th pixel in each axis - 2D/3D case based on last argument existence
int Find_pixel(initializer_list<int> npix, int nnpix_x, int nnpix_y, int nnpix_z)
{
	vector<int>NPIX(npix.begin(),npix.end());
	
	if(nnpix_z!=-1000000) return NPIX[0]*NPIX[1]*nnpix_z+NPIX[0]*nnpix_y+nnpix_x;
	else return NPIX[0]*nnpix_y+nnpix_x;
}



//finding which pixel in each axis nnpix_{x,y,z} is pixel thispix - 3D case
void Find_pixel(initializer_list<int> npix, int thispix, int &nnpix_x, int &nnpix_y, int &nnpix_z)
{
	vector<int>NPIX(npix.begin(),npix.end());
    nnpix_z=thispix/(NPIX[0]*NPIX[1]);
    nnpix_y=(thispix-nnpix_z*NPIX[0]*NPIX[1])/NPIX[0];
    nnpix_x=thispix-nnpix_z*NPIX[0]*NPIX[1]-nnpix_y*NPIX[0];
	return;
}




//finding which pixel in each axis nnpix_{x,y} is pixel thispix - D case
void Find_pixel(initializer_list<int> npix, int thispix, int &nnpix_x, int &nnpix_y)
{
    vector<int>NPIX(npix.begin(),npix.end());
	nnpix_y=thispix/NPIX[0];
	nnpix_x=thispix - nnpix_y*NPIX[0];
	return;
}








//least square method
double LSM(vector<double> &x, vector<double> &y, double &a, double &ua, double &b, double &ub)
{
	double S=0.,Sx=0.,Sy=0.,Sxy=0.,Sxx=0,Syy=0,xi2=0.,delta;
	double er;
	int n=x.size();
	
	for(int i=0;i<n;++i)
	{
		Sx+=x[i];
		Sy+=y[i];
		Sxy+=x[i]*y[i];
		Sxx+=pow(x[i],2);
		Syy+=pow(y[i],2);
	}
	S=n;
	delta=S*Sxx-pow(Sx,2);
	a=(S*Sxy-Sx*Sy)/delta;
	b=(Sxx*Sy-Sx*Sxy)/delta;
	ua=pow(S/(S-2) * (Syy-a*Sxy-b*Sy)/delta,0.5);
	ub=pow(pow(ua,2)*Sxx/S,0.5);
	
	for(int i=0;i<n;++i)
	{
		er=y[i]-a*x[i]-b;
		xi2+=er*er;
	}
	xi2/=1.*n;
	
	return xi2;
}






