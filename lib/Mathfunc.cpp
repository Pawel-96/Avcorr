#include "Mathfunc.h"


double Dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
}




double Dist(double x1, double y1, double x2, double y2, string option)
{
	if(option=="cartesian")
	{
		return sqrt(pow(x1-x2,2)+pow(y1-y2,2));
	}
	
	if(option=="astro") //astronomical: ra,dec [deg], result in deg too
	{
		double dec1r=y1*deg2rad,dec2r=y1*deg2rad;
		return acos(sin(dec1r)*sin(dec2r)+cos(dec1r)*cos(dec2r)*cos((x1-x2)*deg2rad))*rad2deg;
	}
	
	if(option=="spherical") //spherical: ptheta,phi [rad], result in [rad]
	{
		return acos(cos(x1)*cos(x2)+sin(x1)*sin(x2)*cos(y1-y2));
	}
	
	return -1;
}





double Meanlog(double x1, double x2) //logarithmic mean
{
	double power1=log10(x1);
	double power2=log10(x2);
	return pow(10.0,0.5*(power1+power2));
}




int SGN(double a) //signum function
{
	if(a>0) return 1;
	if(a<0) return -1;
	if(a==0) return 0;
	else return 0;
}







//angle [rad] between 2 vectors
double Vector_angle(double p1, double p2, double p3, double q1, double q2, double q3) 
{
	double norm1=sqrt(p1*p1+p2*p2+p3*p3);
    double norm2=sqrt(q1*q1+q2*q2+q3*q3);
	double cosvec=(p1*q1+p2*q2+p3*q3)/(norm1*norm2+1e-15);
	double rotx,roty,rotz;
	rotx=p2*q3-p3*q2;
	roty=q1*p3-p1*q3;
	rotz=p1*q2-p2*q1;

	double sinvec=sqrt(rotx*rotx+roty*roty+rotz*rotz)/(norm1*norm2+1e-15);
	double angle=acos(cosvec);
	if(sinvec<0.){angle=M_PI-angle;}
	return angle;	
}




//length of vector p projected along q
double Projected_vector(double p1, double p2, double p3, double q1, double q2, double q3) 
{
	double pp=sqrt(p1*p1+p2*p2+p3*p3);
	double cosvec=(p1*q1+p2*q2+p3*q3)/(sqrt(q1*q1+q2*q2+q3*q3)*pp+1e-15); //cos of angle between p and q
	return pp*cosvec;
}






double Phi(double x, double y) //polar angle on x,y plane [rad]
{
	if (x>0 and y>=0) return atan(1.0*y/(1.0*x));
	if (x>0 and y<0)   return atan(1.0*y/(1.0*x))+2.0*M_PI;
	if (x<0)           return atan(1.0*y/(1.0*x))+1.0*M_PI;
	if (x==0 and y>0)  return 0.5*M_PI;
	if (x==0 and y<0)  return 1.5*M_PI;
	if (x==0 and y==0) return 0.0;
	else return 0.0;
}







void Convert_coords(double x, double y, double z, double &a, double &b, double &c, string option) //converting coordinates: 3D->3D
{
	if(option=="cartesian2spherical") //[x,y,z] -> [r,theta,phi]
	{
		a=sqrt(pow(x,2)+pow(y,2)+pow(z,2)); //r
		b=acos(z/a); //theta
		c=Phi(x,y); //phi
		return;
	}
	
	if(option=="spherical2cartesian") //[r,theta,phi] -> [x,y,z]
	{
		a=x*sin(y)*cos(z); //x
		b=x*sin(y)*sin(z); //y
		c=x*cos(y); //z
		return;
	}
	
	if(option=="cartesian2cylindrical") //[x,y,z] -> [rho,phi,z]
	{
		a=sqrt(pow(x,2)+pow(y,2)); //rho
		b=Phi(x,y); //phi
		c=z; //z
		return;
	}
	
	if(option=="cylindrical2cartesian") //[rho,phi,z] -> [x,y,z]
	{
		a=x*cos(y); //x
		b=x*sin(y); //y
		c=z; //z
		return;
	}
	
	return;
}









void Convert_coords(double x, double y, double z, double &a, double &b, string option) //converting coordinates: 3D->2D
{
	if(option=="cartesian2astro") //[x,y,z] -> [ra,dec] [deg]
	{
		b=90. - acos(z/a)*rad2deg; //dec
		a=Phi(x,y)*rad2deg; //ra
		return;
	}
	
	return;
}





void Convert_coords(double x, double y, double &a, double &b, double &c, string option) //converting coordinates: 2D->3D
{
	if(option=="astro2cartesian") //[ra,dec] [deg] -> [x,y,z]
	{
		a=cos(y*deg2rad)*cos(x*deg2rad); //x
		b=cos(y*deg2rad)*cos(x*deg2rad); //y
		c=sin(y*deg2rad); //z
		return;
	}
	
	return;
}






void Convert_coords(double x, double y, double &a, double &b, string option) //converting coordinates: 2D->2D
{
	if(option=="astro2spherical") //[ra,dec] [deg] -> [theta,phi] [rad]
	{
		a=(90.-x)*deg2rad;
		b=x*deg2rad;
		return;
	}
	
	if(option=="spherical2astro") // [theta,phi] [rad] -> [ra,dec] [deg]
	{
		a=y*rad2deg;
		b=90. -x*rad2deg;
		return;
	}
	
	return;
}






double Area_spherical(double ramin, double ramax, double decmin, double decmax) //[srd]returns area of sphere fragment (radius=1)
{
    double thetamin=0.5*M_PI-decmax*deg2rad,thetamax=0.5*M_PI-decmin*deg2rad;
    double phimin=ramin*deg2rad,phimax=ramax*deg2rad;
    return (phimax-phimin)*(cos(thetamin)-cos(thetamax));
}








// Helper to make static_assert work with types
template <typename>
struct always_false : false_type {};



/*
***********Check if returning -1000000 in T functions doesnt create problems if T!=int!----------------
*/
