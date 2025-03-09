#ifndef ARRVEC_H
#define ARRVEC_H

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <algorithm>

#include "Const.h"

using namespace std;


double Random(double a, double b); //random number [a,b] range [srand(time(NULL)) required before]


template <typename T>
void Print(vector<T> &vec) //printing vector
{
	int s=vec.size();
	for(int i=0;i<s;++i)
	{
		cout<<vec[i]<<" ";
	}
	return;
}




template <typename T>
void Print(T *arr, int s) //printing array
{
	for(int i=0;i<s;++i)
	{
		cout<<arr[i]<<" ";
	}
	return;
}

//template <typename T>
//double Sum(T &tab, int w=0, int h=0); //sum of elements, works for 1D or 2D array and vector




//summing 1D vector
template <typename T>
double Sum(vector<T> &tab) //sum of elements, works for 1D or 2D array and vector
{
	double sum=0.;
	int s=tab.size();
	for(int i=0;i<s;++i)
	{
		sum+=tab[i];
	}
	
	return sum;
}





//summing 1D array
template <typename T>
double Sum(T *tab, int w)
{
	double sum=0.;
	for(int i=0;i<w;++i)
	{
		sum+=tab[i];
	}
	
	return sum;
}





//summing 2D array
template <typename T>
double Sum(T **tab, int w, int h)
{
	double sum=0.;
	for(int i=0;i<w;++i)
	{
		for(int j=0;j<h;++j)
		{
			sum+=tab[i][j];
		}
	}
	
	return sum;
}







template <typename T> //spacing for vectors
void Spacing(vector<T> &tab, double a, double b, int n, string option="lin") //allows either array or vector
{
	double dsindec=(sin(b*deg2rad)-sin(a*deg2rad))/n;
	tab.resize(n);
	for(int i=0;i<n;++i)
	{
		if(option=="lin") tab[i]=a+1.0*i*(b-a)/n;
		if(option=="log") tab[i]=pow(10,log10(a)+1.0*i*log10(b/a)/n);
		if(option=="fullranges") tab[i]=a+1.0*i*(b-a)/(n-1.);
		if(option=="dsin") tab[i]=asin(sin(a*deg2rad)+i*dsindec)*rad2deg;
	}
	return;
}





template <typename T> //spacing for array
void Spacing(T *tab, double a, double b, int n, string option="lin") //allows either array or vector
{
	double dsindec=dsindec=(sin(b*deg2rad)-sin(a*deg2rad))/n;

	for(int i=0;i<n;++i)
	{
		if(option=="lin") tab[i]=a+1.0*i*(b-a)/n;
		if(option=="log") tab[i]=pow(10,log10(a)+1.0*i*log10(b/a)/n);
		if(option=="fullranges") tab[i]=a+1.0*i*(b-a)/(n-1.);
		if(option=="dsin") tab[i]=asin(sin(a*deg2rad)+i*dsindec)*rad2deg;
	}
	return;
}





template <typename T>
vector<T> Spacing(double a, double b, int n, string option="lin") //n-elements in space [a,b] (lin-linear,log-logarithmic)
{
    vector<T> tab(n);
	double dsindec=dsindec=(sin(b*deg2rad)-sin(a*deg2rad))/n;
	
	for(int i=0;i<n;++i)
	{
        if(option=="lin") tab[i]=a+1.0*i*(b-a)/n;
		if(option=="log") tab[i]=pow(10,log10(a)+1.0*i*log10(b/a)/n);
		if(option=="fullranges") tab[i]=a+1.0*i*(b-a)/(n-1.);
		if(option=="dsin") tab[i]=asin(sin(a*deg2rad)+i*dsindec)*rad2deg;
	}
	return tab;
}







template <typename T>
void Arange(T *tab, double a, double b, double delta) //spacing from a to b, with step=delta (arrays)
{
	int n=(b-a)/delta +1;
	for(int i=0;i<n;++i)
	{
		tab[i]=a+i*delta;
	}
	return;
}





template <typename T>
void Arange(vector<T> &tab, double a, double b, double delta) //spacing from a to b, with step=delta (vectors)
{
	int n=(b-a)/delta +1;
	tab.resize(n);
	for(int i=0;i<n;++i)
	{
		tab[i]=a+i*delta;
	}
	return;
}







template <typename T>
vector<T> Arange(double a, double b, double delta) //spacing from a to b, with step=delta
{
	int n=(b-a)/delta +1;
    vector<T> tab(n);
	tab.resize(n);
	for(int i=0;i<n;++i)
	{
		tab[i]=a+i*delta;
	}
	return tab;
}






//Min or Max of  array 1D/2D/vector; overloaded function
template <typename T>
T Minmax(T *tab, int s, string option)
{
	T res=tab[0];
    if(option=="min")
    {
        for(int i=1;i<s;++i){if(tab[i]<res){res=tab[i];}}
    }
    if(option=="max")
    {
        for(int i=1;i<s;++i){if(tab[i]>res){res=tab[i];}}
    }
    return res;
}




//case for 2D array
template <typename T>
T Minmax(T **tab, int w, int h, string option)
{
	auto res=tab[0][0];
    if(option=="min")
    {
		for(int i=0;i<w;++i)
		{
			for(int j=0;j<h;++j)
			{
				if(tab[i]<res){res=tab[i];}
			}
		}
    }
    if(option=="max")
    {
		for(int i=0;i<w;++i)
		{
			for(int j=0;j<h;++j)
			{
				if(tab[i]>res){res=tab[i];}
			}
		}
    }
    return res;
}




template <typename T>
T Minmax(vector<T> &tab,  string option)
{
	T res=tab[0];
	int ntab=tab.size();
    if(option=="min")
    {
        for(int i=1;i<ntab;++i){if(tab[i]<res){res=tab[i];}}
    }
    if(option=="max")
    {
        for(int i=1;i<ntab;++i){if(tab[i]>res){res=tab[i];}}
    }
    return res;
}





template <typename T>
void Unique(vector<T> &tab)
{
    int n=tab.size();
	int check;
	vector<T> pom=tab;
	
	tab.clear();
	tab.push_back(pom[0]);

    for(int i=1;i<n;++i)
    {
        check=0;
        for(int j=0;j<i;++j)
        {
            if(pom[j]==pom[i]){++check;break;}
        }
        if(check==0){tab.push_back(pom[i]);}
    }
    return;
}





template <typename T>
void Change_index(T &tab, int i, int j) //change index of array or vector - 1D case
{
	auto pom=tab[i];
	tab[i]=tab[j];
	tab[j]=pom;
	return;
}




template <typename T>
void Change_index(T &tab, int i, int j, int k, int l) //change index of array or vector - 2D case
{
	auto pom=tab[i][j];
	tab[i][j]=tab[k][l];
	tab[k][l]=pom;
	return;
}




template <typename T>
T **Empty(int w, int h) //empty array
{
	T **tab = new T*[w];
    for (int i=0;i<w;++i)
	{
        tab[i]=new T[h];
    }
	return tab;
}





template <typename T>
T **Zeros(int w, int h)
{
	T **tab = new T*[w];
    for (int i=0;i<w;++i)
	{
        tab[i]=new T[h];
		for(int j=0;j<h;++j)
		{
			tab[i][j]=0.;
		}
    }
	return tab;
}



template <typename T>
T *Zeros(int s)
{
	T *tab=new T[s];
	for(int i=0;i<s;++i) tab[i]=0;
	return tab;
}






template <typename T> //transposing vector
void Transpose(vector<vector<T>> &vec)
{
	int h=vec.size();
	int w=vec[0].size();
	vector<vector<T>> pom(w,vector<T>(h));
	for(int i=0;i<w;++i)
	{
		for(int j=0;j<h;++j) pom[i][j]=vec[j][i];
	}
	vec.clear();
	vec=pom;
	return;
}



template <typename T> //transposing vector
T **Transpose(T **tab, int h, int w) //transposing array
{
	T **transposed=Zeros<T>(w,h);
	for(int i=0;i<w;++i)
	{
		for(int j=0;j<h;++j) transposed[i][j]=tab[j][i];
	}
	return transposed;
}





template <typename T> //randomizing 1D data - arrays
void Randomize(T &data, int n)
{
	for(int i=0;i<n;++i) Change_index(data,Random(0,n),Random(0,n));
	return;
}


template <typename T> //randomizing 1D data - vectors
void Randomize(T &data)
{
	int n=data.size();
	for(int i=0;i<n;++i) Change_index(data,Random(0,n),Random(0,n));
	return;
}



template <typename T, typename U> //adding value - works both for vectors and arrays
void Add(T &data, int s, U value)
{
	for(int i=0;i<s;++i) data[i]+=value;
	return;
}


template <typename T> //adding array/vec
void Add(T &data, int s, T &values)
{
	for(int i=0;i<s;++i) data[i]+=values[i];
	return;
}



template <typename T, typename U> //multiplying by value - works both for vectors and arrays
void Multiply(T &data, int s, U value)
{
	for(int i=0;i<s;++i) data[i]*=value;
	return;
}


template <typename T> //multiplying array/vec
void Multiply(T &data, int s, T &values)
{
	for(int i=0;i<s;++i) data[i]*=values[i];
	return;
}



template <typename T, typename U> //multiplying by value - works both for vectors and arrays
void Divide(T &data, int s, U value)
{
	for(int i=0;i<s;++i) data[i]/=value;
	return;
}


template <typename T> //multiplying array/vec
void Divide(T &data, int s, T &values)
{
	for(int i=0;i<s;++i) data[i]/=values[i];
	return;
}




template <typename T> //normalize vector [array/vec]
void Normalize(T &data, string dim="3D")
{
	if(dim=="3D")
	{
		auto norm=sqrt(pow(data[0],2)+pow(data[1],2)+pow(data[2],2));
		data[0]/=norm;data[1]/=norm;data[2]/=norm;
	}
	if(dim=="2D")
	{
		auto norm=sqrt(pow(data[0],2)+pow(data[1],2));
		data[0]/=norm;data[1]/=norm;
	}
	return;
}




template <typename T, typename Func> //function on array/vector
void Func_arrvec(T &data, int s, Func func)
{
	for(int i=0;i<s;++i) data[i]=func(data[i]);
	return;
}




template <typename T> //rearranging irregular data: x[i] will be equally spaced [x must be sorted]
void Rearrange_data(vector<T> &x, vector<T> &y, string interpolation="linear") 
{
    int n=x.size();
    auto Minx=x[0], Maxx=x[n-1];
	auto stepx=1.0*(Maxx-Minx)/n;
    vector<T> p,q;
    p=x;
    q=y;
    
    for(int i=0;i<n;++i)
    {
        x[i]=Minx+1.*i*stepx;
        y[i]=Interpolate(p,q,n,x[i],interpolation);
    }
    return;
}




template <typename T> //sorting data vector based on ncol-th column; data[col][row]
void Sort(vector<vector<T> >&data, int ncol)
{
	int nrows=data[0].size(),ncols=data.size();
	vector<int> indices=Arange<int>(0,nrows-1,1);
	vector<T> col_tosort=data[ncol]; //column based on which sorting occurs
	
	sort(indices.begin(),indices.end(), [&col_tosort](int i, int j)
	{return col_tosort[i]<col_tosort[j];}); //sorting indices based on data[ncol] vector
	
	vector<vector<T> > pom=data; //copying data
	
	for(int i=0;i<ncols;++i)
	{
		for(int j=0;j<nrows;++j) data[i][j]=pom[i][indices[j]]; //sorting based on indices
	}
	pom.clear();
	return;
}




template <typename T> //sorting data array based on ncol-th column; data[col][row]
void Sort(T **data, int ncols, int nrows, int ncol)
{
	int *indices=new int[nrows];
	Arange(indices,0,nrows-1,1);
	
	T *col_tosort=new T[nrows]; //column based on which sorting occurs
	for(int i=0;i<nrows;++i) col_tosort[i]=data[ncol][i];
	
	sort(indices,indices+nrows, [col_tosort](int i, int j)
	{return col_tosort[i]<col_tosort[j];}); //sorting indices based on data[ncol] vector
	
	T **pom=Zeros<T>(ncols,nrows);
	for(int i=0;i<ncols;++i)
	{
		for(int j=0;j<nrows;++j) pom[i][j]=data[i][j];
	}
	
	for(int i=0;i<ncols;++i)
	{
		for(int j=0;j<nrows;++j) data[i][j]=pom[i][indices[j]]; //sorting based on indices
	}
	return;
}



#endif