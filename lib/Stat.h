#ifndef STAT_H
#define STAT_H

#include <numeric>

#include "Mathfunc.h"
#include "Arrvec.h"



template <typename T>
double Corr(vector<T> &x, vector<T> &y) //correlation coefficient
{
    int n=x.size();
    double avx=Sum<T>(x)/(1.0*n);
    double avy=Sum<T>(y)/(1.0*n);

    double sumxy=0.0,sumxx=0.0,sumyy=0.0,res;
    for(int i=0;i<n;++i)
    {
        sumxy+=(x[i]-avx)*(y[i]-avy);
        sumxx+=pow(x[i]-avx,2);
        sumyy+=pow(y[i]-avy,2);
    }
    res=sumxy/sqrt(sumxx*sumyy);
    return res;
}






template <typename T>
double Sigma(T *tab, int s) //Gaussian deviation on dataset (array)
{
    double sigma=0., av=Sum<T>(tab,s)/(1.*s);
    for(int i=0;i<s;++i)
    {
        sigma+=pow(tab[i]-av,2);
    }
    sigma/=(double)(1.0*s-1.0);
    return sqrt(sigma);
}





template <typename T>
double Sigma(vector<T> &tab) //Gaussian deviation on dataset (vector)
{
	int s=tab.size();

    double sigma=0., av=Sum<T>(tab)/(1.*s);
    for(int i=0;i<s;++i)
    {
        sigma+=pow(tab[i]-av,2);
    }
    sigma/=(double)(1.0*s-1.0);
    return sqrt(sigma);
}






template <typename T> 
T Gauss(T x, T *par, string option) //Gaussian distribution value
{
    T miu=par[0],sigma=par[1];
	if(option=="normal") return exp(-0.5*pow((x-miu)/sigma,2))/(sigma*sqrt(2.0*M_PI));
	if(option=="max1") return exp(-0.5*pow((x-miu)/sigma,2));
	else return -1;
}





template <typename T> 
T Gauss_max1(T x, T *par) //Gaussian distribution, but max1
{
    T miu=par[0],sigma=par[1];
	return exp(-0.5*pow((x-miu)/sigma,2))/(sigma*sqrt(2.0*M_PI));

}

string Random_str(int n); //random string of length n

int Chose_probability(double p); //returns 1 if chosen with probability p, otherwise 0


template <typename T> //random position - 2D
void Random_position(initializer_list<T> params, T &ra, T &dec, string option)
{
	vector<T> pars(params.begin(),params.end());
	double p;
	
	if(option=="sky") //params: {ramin,ramax,decmin,decmax} [deg]
	{
		ra=Random(pars[0],pars[1]); //[ramin,ramax]
		for(;;)
		{
			dec=Random(pars[2],pars[3]); //[decmin,decmax]
			p=cos(dec*deg2rad);
			if(Chose_probability(p)==1){break;}
		}
		return;
	}
	
	if(option=="sky_border") //params: {ramin,ramax,decmin,decmax,delta} [deg]
	{
		double dec_raminclose,dec_ramaxclose; //distance from circle center to furthest point in pixel
		for(;;)
		{
			ra=Random(pars[0]+pars[4],pars[1]-pars[4]);
			for(int k=0;;++k)
			{
				dec=Random(pars[2]+pars[4],pars[3]-pars[4]);
				p=cos(dec*deg2rad);
				if(Chose_probability(p)==1){break;}
			}

			if(sin(dec*deg2rad)==0.0){dec_raminclose=0.0;dec_ramaxclose=0.0;}
			else
			{
				//declination of point on ra=ramin line which is closests to center of drawn circle
				dec_raminclose=atan(tan(0.5*M_PI-dec*deg2rad)*cos((pars[0]-ra)*deg2rad))*rad2deg;
				//declination of point on ra=ramax line which is closests to center of drawn circle
				dec_ramaxclose=atan(tan(0.5*M_PI-dec*deg2rad)*cos((pars[1]-ra)*deg2rad))*rad2deg;
			}
			if(Dist(ra,dec,pars[0],dec_raminclose,"spherical")>pars[4] and Dist(ra,dec,pars[1],dec_ramaxclose)>pars[4]){break;} //(for DEC is ok)
		}
	return;
	}
	return;
}





template <typename T> //random position - 3D
void Random_position(initializer_list<T> params, T &X, T &Y, T &Z, string option)
{
	vector<T> pars(params.begin(),params.end());
	
	if(option=="shell") //params: {rmin,rmax}
	{
		T r=Random(pars[0],pars[1]); //dist from center [0,0,0]
		double p,theta,phi=Random(0.,2.*M_PI);
		
		for(;;)
		{
			theta=Random(0.,M_PI);
			p=sin(theta);
			if(Chose_probability(p)==1){break;}
		}
		X=r*sin(theta)*cos(phi);
		Y=r*sin(theta)*sin(phi);
		Z=r*cos(theta);
		return;
	}
	
	if(option=="box") //params: {boxsize}
	{
		X=Random(0.,pars[0]);
		Y=Random(0.,pars[0]);
		Z=Random(0.,pars[0]);
		return;
	}
	
	return;
}





template <typename T> //random value from a 1D custom distribution ( f must be max1)
T Random_distr(T(*f)(T,T*),T *par, T a, T b)
{
    T x,p;
    int chosen;
    while(true)
    {
        x=Random(a,b);
        p=f(x,par);
        chosen=Chose_probability(p);
        if(chosen==1){return x;}
    }
	return -1000000;
}





template <typename T> //random value from a 1D custom distribution defined by vector (max1)
T Random_distr(vector<T> &arg, vector<T> &hist, int n, T a, T b, string interpolate_option="linear")
{
    T x,p;
    int chosen;
    while(true)
    {
        x=Random(a,b);
        p=Interpolate(arg,hist,n,x,interpolate_option);
        chosen=Chose_probability(p);
        if(chosen==1){return x;}
    }
	return -1000000;
}






template <typename T> //setting [xres,yres] to random custom 2D distribution defined by function f
void Random_distr(T (*f)(T,T,T*), T *par, T xmin, T xmax, T ymin, T ymax, T &xres, T &yres)
{
    double p,is_ok; //probability of chosing, is ok?
    while(true)
    {
        xres=Random(xmin,xmax);
        yres=Random(ymin,ymax);

        p=f(xres,yres,par)/f(0.0,0.0,par);
        is_ok=Chose_probability(p);
        if(is_ok==1){return;}
    }
    return;
}





template <typename T, typename L>
int Hist(vector<T> &data, vector<T> &arg, vector<L> &hist, string option, double a=0, double b=0, int n=-1) //histogram
{
	int dsize=data.size();
	if constexpr (is_same_v<decay_t<T>, vector<string>>) //data made of strings; output: string,int/float/double
	{
		arg=data;
		Unique(arg);
		hist.clear();
		hist.resize(arg.size(),0); //filling with zeros
		for(int i=0;i<dsize;++i)
		{
			for(int j=0;j<arg.size();++j)
			{
				if(data[i]==arg[j])
				{
					hist[j]+=1;
					break;
				}
			}
		}
		return 0;
    }
	
	if(a==0 and b==0) //ranges not set, finding
	{
		a=Minmax<T>(data,"min");
		b=Minmax<T>(data,"max");
	}
	
	int loc,outliers=0;
	if(option=="discrete")
	{
		Arange(arg,a,b,1);
		int nnarg=arg.size();
		hist.clear();
		hist.resize(nnarg,0); //filling with zeros
		
		for(int i=0;i<dsize;++i)
		{
			loc=data[i]-a;
			if(loc<0 or loc+1>nnarg) {outliers+=1; continue;} //out of range
			hist[loc]+=1.;
		}
	}
	
	else
	{
		double delta=(b-a)/(1.*n);
		Spacing<T>(arg,a,b,n);
		hist.clear();
		hist.resize(arg.size(),0); //filling with zeros
		
		for(int i=0;i<dsize;++i) //counting
		{
			loc=(data[i]-a)/delta;
			if(loc<0 or loc+1>n) {outliers+=1; continue;} //out of range
			hist[loc]+=1.;
		}
		
		//for option=="hist "- not normalizing, for other - normalization
		if(option!="hist")
		{
			double norm=1.;
			if(option=="PDF") norm=Sum<L>(hist)*delta; //normalizing such that Sum(hist_i delta_i) =1
			if(option=="max1") norm=Minmax(hist,"max"); //achieving 1 at most
			
			for(int i=0;i<n;++i){hist[i]/=norm;} //normalizing
		}
	}

	return outliers;
}





template <typename T> //cumulative distribution
vector<T> Cumulative(vector<T> &tab)
{
	int n=tab.size();
	vector<T> CDF(n);
	
	CDF[0]=tab[0];
	for(int i=1;i<n;++i) CDF[i]=CDF[i-1]+tab[i];
	
	return CDF;
}






template <typename T> //converting discrete histogram into 1D list (randomized) - array
T *Hist2list(T *arg, T *count, int s) //args should be int; if other type, int(val)
{
	int nn=Sum<T>(count,s), i_all=0;
	T *list=new T[nn];
	
	for(int i=0;i<s;++i) //firstly adding blocks of the same counts from within one "bin"
	{
		for(int j=0;j<count[i];++j) //number of counts within that bin
		{
			list[i_all]=1.*arg[i];
			i_all+=1;
		}	
	}
	
	Randomize(list,nn);
	return list;
}




//func overload for faster executing
template <typename T> //converting discrete histogram into 1D list (randomized) - vector
vector<T> Hist2list(vector<T> arg, vector<T> count) //args should be int; if other type, int(val)
{
	int s=arg.size();
	int nn=Sum<T>(count), i_all=0;
	vector<T> list(nn);

	for(int i=0;i<s;++i) //firstly adding blocks of the same counts from within one "bin"
	{
		for(int j=0;j<count[i];++j) //number of counts within that bin
		{
			list[i_all]=1.*arg[i];
			i_all+=1;
		}	
	}
	
	Randomize(list);
	return list;
}






/*pixelizing objects; for 2D case (3D analogical) pixels distribution is:
12  13  14  15
8   9   10  11
4   5   6   7
0   1   2   3
*/
template <typename T>
void Pixelize(initializer_list<T> ranges_list, initializer_list<int> npix_axis, initializer_list< vector<T> > coords_list, vector<vector<vector<T> > > &pix, string option, string msg="show")
{	//pix[i][j][k]; i-number of pixel,j-number of object, k-coordinate and other data
	//example for 3D case pix[5][7][1]=10means that 7th object in 5th pixel has coordinate y=10 (last index: 0,1,2 ->x,y,z)
    
	pix.clear(); //clearing before writing data there
	vector<int> nnpix(npix_axis.begin(), npix_axis.end()); //collecting into vector
	vector<T> ranges(ranges_list.begin(),ranges_list.end()); //collecting ranges to vector
	
	vector<vector<T>> coord(coords_list); //coordinates
	for (const auto& vec:coords_list) coord.push_back(vec); //obtaining form list

	int ndim=nnpix.size(); //dimension
	int nobj=coord[0].size(); //number of objects to pixelize
	int nval=coord.size(); //number of values: coordinates + optional (no of optional=nval-ndim)
	int npix=accumulate(nnpix.begin(),nnpix.end(), 1.0, multiplies<int>()); //number of all pixels 
	int nax[ndim],foundpix; //no. of pixel in each axis, sum of objects in this pixel
	double dpix[ndim]; //pixel size in each dimension (for sky - delta_ra, delta_sin(dec))
	
	if(msg=="show")
	{
		cout<<"Pixelizing: ";
		for(int i=0;i<ndim;++i) cout<<"["<<ranges[2*i]<<":"<<ranges[2*i+1]<<"] ";
		cout<<"npix: ";
		for(int i=0;i<ndim;++i) cout<<nnpix[i]<<" ";
		cout<<endl;
	}
	
	vector<T> obj(nval,-1000000); //vector for one object
	vector<vector<T> >onepix{obj}; //vector for one pixel - filled with default object "empty"
	pix.assign(npix,onepix); //adding "empty" pixels
	
	if(option=="3D") //3D case: [x,y,z]
	{
		dpix[0]=(ranges[1]-ranges[0])/(1.*nnpix[0]); //pixel size in x axis
		dpix[1]=(ranges[3]-ranges[2])/(1.*nnpix[1]); //pixel size in y axis
		dpix[2]=(ranges[5]-ranges[4])/(1.*nnpix[2]);//pixel size in z axis
		
		for(int i=0;i<nobj;++i)
		{
			if(coord[0][i]<ranges[0] || coord[0][i]>ranges[1] ||
			   coord[1][i]<ranges[2] || coord[1][i]>ranges[3] ||
			   coord[2][i]<ranges[4] || coord[2][i]>ranges[5]) continue; //object out of ranges
			
			nax[0]=(coord[0][i]-ranges[0])/dpix[0]; //number of pixel in x axis
			nax[1]=(coord[1][i]-ranges[2])/dpix[1]; //number of pixel in y axis
			nax[2]=(coord[2][i]-ranges[4])/dpix[2]; //number of pixel in z axis
			foundpix=nax[2]*nnpix[1]*nnpix[0] +nax[1]*nnpix[0] +nax[0]; //number of pixel which this object will belong to
			
			if(pix[foundpix][0][0]==-1000000) //pixel doesn't have any member yet
			{
				for(int j=0;j<nval;++j) pix[foundpix][0][j]=coord[j][i]; //filling "empty" pixel
			}
			else //pixel has some data already
			{
				for(int j=0;j<nval;++j) obj[j]=coord[j][i]; //filling object vector with this object data
				pix[foundpix].push_back(obj); //adding object to pixel
			}
		}
	}
	
	if(option=="2D") //2D case: [x,y]
	{
		dpix[0]=(ranges[1]-ranges[0])/(1.*nnpix[0]); //pixel size in x axis
		dpix[1]=(ranges[3]-ranges[2])/(1.*nnpix[1]); //pixel size in y axis
		
		for(int i=0;i<nobj;++i)
		{
			if(coord[0][i]<ranges[0] || coord[0][i]>ranges[1] ||
			   coord[1][i]<ranges[2] || coord[1][i]>ranges[3]) continue; //object out of ranges
			
			nax[0]=(coord[0][i]-ranges[0])/dpix[0]; //number of pixel in x axis
			nax[1]=(coord[1][i]-ranges[2])/dpix[1]; //number of pixel in y axis
			foundpix=nax[1]*nnpix[0] +nax[0]; //number of pixel which this object will belong to
			
			if(pix[foundpix][0][0]==-1000000) //pixel doesn't have any member yet
			{
				for(int j=0;j<nval;++j) pix[foundpix][0][j]=coord[j][i]; //filling "empty" pixel
			}
			else //pixel has some data already
			{
				for(int j=0;j<nval;++j) obj[j]=coord[j][i]; //filling object vector with this object data
				pix[foundpix].push_back(obj); //adding object to pixel
			}
		}
	}
	
	if(option=="2D_sky") //sky case: [ra,dec] [deg]
	{
		dpix[0]=(ranges[1]-ranges[0])/(1.*nnpix[0]); //pixel size in ra
		dpix[1]=(sin(ranges[3]*deg2rad)-sin(ranges[2]*deg2rad))/(1.*nnpix[1]); //pixel size in sin(dec)
		
		for(int i=0;i<nobj;++i)
		{
			if(coord[0][i]<ranges[0] || coord[0][i]>ranges[1] ||
			   coord[1][i]<ranges[2] || coord[1][i]>ranges[3]) continue; //object out of ranges
			
			nax[0]=(coord[0][i]-ranges[0])/dpix[0]; //number of pixel in ra
			nax[1]=(sin(coord[1][i]*deg2rad)-sin(ranges[2]*deg2rad))/dpix[1]; //number of pixel in dec
			foundpix=nax[1]*nnpix[0] +nax[0]; //number of pixel which this object will belong to
			
			if(pix[foundpix][0][0]==-1000000) //pixel doesn't have any member yet
			{
				for(int j=0;j<nval;++j) pix[foundpix][0][j]=coord[j][i]; //filling "empty" pixel
			}
			else //pixel has some data already
			{
				for(int j=0;j<nval;++j) obj[j]=coord[j][i]; //filling object vector with this object data
				pix[foundpix].push_back(obj); //adding object to pixel
			}
		}
	}

	if(msg=="show") cout<<"Pixelized succesfully"<<endl;
    return;
}






//finding which pixel is nnpix_{x,y,z}-th pixel in each axis - 2D/3D case based on last argument existence
int Find_pixel(initializer_list<int> npix, int nnpix_x, int nnpix_y, int nnpix_z=-1000000);

//finding which pixel in each axis nnpix_{x,y,z} is pixel thispix - 3D case
void Find_pixel(initializer_list<int> npix, int thispix, int &nnpix_x, int &nnpix_y, int &nnpix_z);

//finding which pixel in each axis nnpix_{x,y} is pixel thispix - D case
void Find_pixel(initializer_list<int> npix, int thispix, int &nnpix_x, int &nnpix_y);


//***add more options*******************************************************************************
template <typename T> //cleaning vector of data based on option
void Clean_data(vector<T> &vec, string option)
{
	vector<T> pom;
	int s=vec.size();
	if(option=="erase<=0")
	{
		for(int i=0;i<s;++i)
		{
			if(vec[i]>0.) pom.push_back(vec[i]);
		}
	}
	vec.clear();
	vec=pom;
	return;
}






//combining results into estimate +/- err; data[npoints][nreals] (errs computed for sample of nreals measurements)
template <typename T> 
void Combine_results(vector<vector<T> > &data, vector<T> &res, vector<T> &err, string estimation, string option,
int nsubsamples=1, int subsample_size=1)
{
	int npoints=data.size(); //number of independent measurements (not for err calculation)
	int nreals; //number of realizations - measurement repetitions
	vector<T> datapoint; //storing here different realizations of the same measurement
	
	res.resize(npoints);
	err.resize(npoints);
	
	for(int i=0;i<npoints;++i)
	{
		datapoint=data[i];
		if(option!="") Clean_data(datapoint,option);
		nreals=datapoint.size(); //valid realizations after cleaning
		
		if(estimation=="Gauss")
		{
			res[i]=Sum<T>(datapoint)/(1.*nreals);
			err[i]=Sigma(datapoint);
		}
		
		if(estimation=="jackknife")
		{
			vector<T> means(nreals); //Jackknife means
			for(int j=0;j<nreals;++j) means[j]=(Sum<T>(datapoint)-datapoint[j])/(nreals-1.);
			res[i]=Sum(means)/(1.*nreals);
			err[i]=Sigma(means)*nreals;
			
		}
		
		if(estimation=="bootstrap")
		{
			vector<T> means(nsubsamples,0.);
			for(int j=0;j<nsubsamples;++j) //every subsample
			{
				//drawing random subsample with returning:
				for(int k=0;k<subsample_size;++k) means[j]+=datapoint[(int)(Random(0,nreals))];
				means[j]/=1.*subsample_size;
			}
			res[i]=Sum<T>(means)/(1.*nsubsamples);
			err[i]=Sigma(means);
			
		}
	}
	return;
}


//least square method
double LSM(vector<double> &x, vector<double> &y, double &a, double &ua, double &b, double &ub);




#endif