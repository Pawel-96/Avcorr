#ifndef MATHFUNC_H
#define MATHFUNC_H

#include <iostream>
#include <vector>
#include <math.h>

#include "Const.h"

using namespace std;


double Dist(double x1, double y1, double z1, double x2, double y2, double z2); //Distance in 3D cartesian coordinates
double Dist(double x1, double y1, double x2, double y2, string option="cartesian"); //Distance in 2D; option: cartesian/astro/spherical

double Meanlog(double x1,double x2); //logarithmic mean

int SGN(double a); //signum function


template <typename T>
T Round_cut(T val, int n) //rounding-cut for n-digits: Round_cut(5.378,2)=5.37
{
    T res=floor(val*pow(10,n));
    res/=pow(10,n);
    return res;
}

template<typename T> //vector product pxq - vector elements as params
void Vector_product(T p1, T p2, T p3, T q1, T q2, T q3, T &vec)
{
	vec[0]=p2*q3-p3*q2;
	vec[1]=q1*p3-p1*q3;
	vec[2]=p1*q2-p2*q1;
	return;
}




template<typename T> //vector product - vectors or arrays as params
void Vector_product(T &p, T &q, T &vec)
{
	vec[0]=p[1]*q[2]-p[2]*q[1];
	vec[1]=q[0]*p[2]-p[0]*q[2];
	vec[2]=p[0]*q[1]-p[1]*q[0];
	return;
}

//angle [rad] between 2 vectors
double Vector_angle(double p1, double p2, double p3, double q1, double q2, double q3);

//length of vector p projected along q
double Projected_vector(double p1, double p2, double p3, double q1, double q2, double q3);

double Phi(double x, double y); //polar angle on x,y plane [rad]

void Convert_coords(double x, double y, double z, double &a, double &b, double &c, string option="cartesian2spherical"); //converting coordinates: 3D->3D
void Convert_coords(double x, double y, double z, double &a, double &b, string option); //converting coordinates: 3D->2D
void Convert_coords(double x, double y, double &a, double &b, double &c, string option); //converting coordinates: 2D->3D
void Convert_coords(double x, double y, double &a, double &b, string option); //converting coordinates: 2D->2D






template <typename T> //rotating in 3D: axis=X/Y/Z, [ang]=deg
void Rotate(T &x, T &y, T &z, T ang, string axis) 
{
    T theta=ang*deg2rad;
    T pomx=x,pomy=y,pomz=z;

    if(axis=="X")
    {
        y=cos(theta)*pomy-sin(theta)*pomz;
        z=sin(theta)*pomy+cos(theta)*pomz;
    }

    if(axis=="Y")
    {
        x=cos(theta)*pomx+sin(theta)*pomz;
        z=-1.0*sin(theta)*pomx+cos(theta)*pomz;
    }

    if(axis=="Z")
    {
        x=cos(theta)*pomx-sin(theta)*pomy;
        y=sin(theta)*pomx+cos(theta)*pomy;
    }
    return;
}


double Area_spherical(double ramin, double ramax, double decmin, double decmax); //[srd]returns area of sphere fragment (radius=1)


template <typename T> //integrating 1D function
T Integrate(T (*f)(T, T *), T *par, T a, T b,T eps)
{
	T I0,h0=b-a;
	T h1=h0/2;
	I0=h0/2 *(f(a,par)+f(b,par));
	T I1=h1/2 *(f(a,par)+2*f(0.5*(a+b),par)+f(b,par));
	T I;
	if(fabs((I0-I1)/(I0+1e-10))<eps)
	{return I1;} else
	I=Integrate(f,par,a,(a+b)/2,eps)+Integrate(f,par,(a+b)/2,b,eps);
	return I;
}



template <typename T> //integrating 2D function
T Integrate_2D(T (*f)(T,T,T *), T *par, T xmin, T xmax, T ymin, T ymax, T eps)
{
	T h0x=xmax-xmin,h1x=h0x/2;
	T h0y=ymax-ymin,h1y=h0y/2;
	T I0,I1,I;
	I0=h0x*h0y/4. *(f(xmin,ymin,par)+f(xmin,ymax,par)+f(xmax,ymin,par)+f(xmax,ymax,par));
	I1=h1x*h1y/4. *(  f(xmin,ymin,par)+f(xmin,ymax,par)+f(xmax,ymin,par)+f(xmax,ymax,par)+
							  2.*(f(0.5*(xmin+xmax),ymin,par)+f(0.5*(xmin+xmax),ymax,par)+f(xmin,0.5*(ymin+ymax),par)+f(xmax,0.5*(ymin+ymax),par))+
							  4.*f(0.5*(xmin+xmax),0.5*(ymin+ymax),par)  );

	if(fabs((I0-I1)/(I0+1e-010))<eps) return I1;
	else
	{
		I=Integrate_2D(f,par,xmin,0.5*(xmin+xmax),ymin,0.5*(ymin+ymax),eps)+
		Integrate_2D(f,par,0.5*(xmin+xmax),xmax,ymin,0.5*(ymin+ymax),eps)+
		Integrate_2D(f,par,xmin,0.5*(xmin+xmax),0.5*(ymin+ymax),ymax,eps)+
		Integrate_2D(f,par,0.5*(xmin+xmax),xmax,0.5*(ymin+ymax),ymax,eps);
	}
	
	return I;
}






template <typename T> //interpolating - 1D case [arrays]
double Interpolate(T &x, T &f, int n, double x0, string option="linear")
{
	if(option=="linear")
	{
		double x1,x2,y1,y2;

		for(int i=0;i<n-1;++i)
		{
			if(x0>x[i] and x0<=x[i+1])
			{
				x1=x[i];y1=f[i];
				x2=x[i+1];y2=f[i+1];
				break;
			}
		}

		if(x0==x1) return y1;
		else if(x0!=x1) return 1.0*y1+1.0*(((y2-y1)/(x2-x1))*(x0-x1));
		else return f[0];
	}
	
	if(option=="lagrange")
	{
		double sum=0.0,add=1.0;

		for(int i=0;i<n-1;++i)
		{   
			add=f[i];
			for(int j=0;j<n-1;++j)
			{
				if(j==i) continue;
				add*=(x0-x[j])/(x[i]-x[j]+1e-10);
			}
			sum+=add;
		}
		return sum;
	}
	
	if(option=="knownstep") //regular data, linear interpolation
	{
		if(x0<x[0]){return f[0];}
		if(x0>x[n-1]){return f[n-1];}
		double x1,x2,y1,y2,dx=x[1]-x[0];

		x1=x[0]+floor((x0-x[0])/dx)*dx;
		x2=x1+dx;
		int P1=(x1-x[0])/dx; //locations in arrays
		y1=f[P1];
		y2=f[P1+1];

		if(x0==x1) return y1;
		else if(x0!=x1) return 1.0*y1+1.0*(((y2-y1)/(x2-x1))*(x0-x1));
		else return f[0];
	}
	
	if(option=="c")
	{
		double sum=0.,term;
		for(int i=0;i<n;++i)
		{
			term=f[i];
			for(int j=0;j<n;++j)
			{
				if(j!=i){term*=(x0-x[j])/(x[i]-x[j]);}
			}
			sum+=term;
		}
		return sum;
	}
	return -1000000;
}






template <typename T> //interpolating - 1D case [vectors]
double Interpolate(vector<T> &x, vector<T> &f, double x0, string option="linear")
{
	int n=x.size();
	if(option=="linear")
	{
		double x1,x2,y1,y2;

		for(int i=0;i<n-1;++i)
		{
			if(x0>x[i] and x0<=x[i+1])
			{
				x1=x[i];y1=f[i];
				x2=x[i+1];y2=f[i+1];
				break;
			}
		}

		if(x0==x1) return y1;
		else if(x0!=x1) return 1.0*y1+1.0*(((y2-y1)/(x2-x1))*(x0-x1));
		else return f[0];
	}
	
	if(option=="lagrange")
	{
		double sum=0.0,add=1.0;

		for(int i=0;i<n-1;++i)
		{   
			add=f[i];
			for(int j=0;j<n-1;++j)
			{
				if(j==i) continue;
				add*=(x0-x[j])/(x[i]-x[j]+1e-10);
			}
			sum+=add;
		}
		return sum;
	}
	
	if(option=="knownstep") //regular data, linear interpolation
	{
		if(x0<x[0]){return f[0];}
		if(x0>x[n-1]){return f[n-1];}
		double x1,x2,y1,y2,dx=x[1]-x[0];

		x1=x[0]+floor((x0-x[0])/dx)*dx;
		x2=x1+dx;
		int P1=(x1-x[0])/dx; //locations in arrays
		y1=f[P1];
		y2=f[P1+1];

		if(x0==x1) return y1;
		else if(x0!=x1) return 1.0*y1+1.0*(((y2-y1)/(x2-x1))*(x0-x1));
		else return f[0];
	}
	
	if(option=="c")
	{
		double sum=0.,term;
		for(int i=0;i<n;++i)
		{
			term=f[i];
			for(int j=0;j<n;++j)
			{
				if(j!=i){term*=(x0-x[j])/(x[i]-x[j]);}
			}
			sum+=term;
		}
		return sum;
	}
	return -1000000;
}



template <typename T> //linear interpolating - 2D case [arays]
double Interpolate2D(T &x, T &y, T &f, int n, double x0, double y0)
{
	double x1,x2,y1,y2,FQ11,FQ12,FQ21,FQ22,l1,l2;
	for(int i=0;i<n-1;++i)
	{
		if(x[i]<=x0 and x[i+1]>x0 and x[i]!=x[i+1]){x1=x[i];x2=x[i+1];} //not positions yet ->it isn't sorted like 1D data
		if(y[i]<=y0 and y[i+1]>y0 and y[i]!=y[i+1]){y1=y[i];y2=y[i+1];}
	}

	for(int i=0;i<n-1;++i)
	{
		if(x[i]==x1 and y[i]==y1){FQ11=f[i];}
		if(x[i]==x2 and y[i]==y1){FQ21=f[i];}
		if(x[i]==x1 and y[i]==y2){FQ12=f[i];}
		if(x[i]==x2 and y[i]==y2){FQ22=f[i];}
	}
	l1=FQ11*(y2-y0)+FQ12*(y0-y1),l2=FQ21*(y2-y0)+FQ22*(y0-y1);
	return ((x2-x0)*l1+(x0-x1)*l2)/((x2-x1)*(y2-y1)+1e-10);
}





template <typename T> //linear interpolating - 2D case [vectors]
double Interpolate2D(vector<T> &x, vector<T> &y, vector<T> &f, double x0, double y0)
{
	int n=x.size();
	double x1,x2,y1,y2,FQ11,FQ12,FQ21,FQ22,l1,l2;
	for(int i=0;i<n-1;++i)
	{
		if(x[i]<=x0 and x[i+1]>x0 and x[i]!=x[i+1]){x1=x[i];x2=x[i+1];} //not positions yet ->it isn't sorted like 1D data
		if(y[i]<=y0 and y[i+1]>y0 and y[i]!=y[i+1]){y1=y[i];y2=y[i+1];}
	}

	for(int i=0;i<n-1;++i)
	{
		if(x[i]==x1 and y[i]==y1){FQ11=f[i];}
		if(x[i]==x2 and y[i]==y1){FQ21=f[i];}
		if(x[i]==x1 and y[i]==y2){FQ12=f[i];}
		if(x[i]==x2 and y[i]==y2){FQ22=f[i];}
	}
	l1=FQ11*(y2-y0)+FQ12*(y0-y1),l2=FQ21*(y2-y0)+FQ22*(y0-y1);
	return ((x2-x0)*l1+(x0-x1)*l2)/((x2-x1)*(y2-y1)+1e-10);
}






template <typename T> //inversing the function func(x,par)=y; [a,b] - argument ranges for f; y- new argument
T Inverse_function(T (*func)(T,T*), T *par, T a, T b, T eps, T y) //method: bisection
{
	T fa,fb,fc,c;
    int counter=0;
    while(true)
    {
        c=0.5*(a+b);
        fa=func(a,par)-y;
		fb=func(b,par)-y;
		fc=func(c,par)-y;
        if (fa*fc<0) {b=c;}
        if (fc*fb<0) {a=c;}
        if (fabs(fc)<eps) {return c;}
        ++counter;
        if(counter>1000){return -1000000;}
	}
	
}





template <typename T>//inversing the function func(x,par)=y; start - starting point; y- new argument
T Inverse_function(T (*func)(T,T*), T *par, T start, T eps, T y) //method: Newton
{
    T arg=start;
	T fval,deriv;
	int counter=0;
	while(true)
	{
		fval=func(arg,par)-y;
		deriv=(func(arg+eps,par)-func(arg,par))/eps;
		
		arg-=fval/deriv;
		if(counter>1000){return -1000000;}
		if(fval/deriv<eps){break;}
		++counter;
	}
	return arg;
}





template <typename T> //bisection method: finding argument where func=0
T Bisection(T (*f)(T,T*), T *par, T a, T b, T eps) 
{
    T fa,fb,fc,c;
    int counter=0;
    while(true)
    {
        c=0.5*(a+b);
        fa=f(a,par);fb=f(b,par);fc=f(c,par);
        if (fa*fc<0) {b=c;}
        if (fc*fb<0) {a=c;}
        if (fabs(fc)<eps) {return c;}
        ++counter;
        if(counter>1000){return -1000000;}
    }
}





#endif