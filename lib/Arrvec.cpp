#include "Arrvec.h"


double Random(double a, double b) //random number [a,b] range [srand(time(NULL)) required before]
{
	double res=a+(1.0*rand()/RAND_MAX)*(b-a);
	return res;
}





/*
template <typename T>
double Sum(T &tab, int w, int h) //sum of elements, works for 1D or 2D array and vector
{
	double sum=0.;
	if constexpr (is_pointer_v<T>)
    {
        if (constexpr (is_pointer_v<remove_pointer_t<T>>)) //2D array
        {
            for(int i=0;i<w;++i)
			{
				for(int j=0;j<h;++j)
				{
					sum+=tab[i][j];
				}
			}
        }
        else //1D array
        {
            for(int i=0;i<w;++i)
			{
				sum+=tab[i];
			}
        }
    }
    else if (constexpr (is_same_v<decay_t<T>, vector<double>> ||
              is_same_v<decay_t<T>, vector<float>> ||
              is_same_v<decay_t<T>, vector<int>>)) //1D vector case of type int/float/double
    {
        for(int i=0;i<tab.size();++i)
		{
			sum+=tab[i];
		}
    }
    else if (constexpr (is_same_v<decay_t<T>, vector<vector<double>>>)) //2D vector
    {
        for(int i=0;i<tab.size();++i)
		{
			for(int j=0;j<tab[i].size();++j)
			{
				sum+=tab[i][j];
			}
		}
    }
	else
	{
		cout<<"Invalid data type :/"<<endl;
		return -1000000;
	}
	
	return sum;
}
*/



