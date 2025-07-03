#ifndef COMBINE_H
#define COMBINE_H

#include "../lib/RWMathStat.h"
#include "CIC.h"
#include "Parsfnames.h"


int Merge();

//converting 1D arguments back to 2D grid in results (one file)
void Args2grid_onefile(string fname);

void Args2grid(); //converting 1D arguments back to 2D grid in results

//extracting realisation number from modelname based on real_template
int Realisation(string modelname);

//getting independent models if combine_reals=1 (files are models with different reals)
vector<string> Get_independent_models(vector<string> models);

//finding number of realisations (checking if file exists, if not, stopping)
int Maxreal(string modelname); //[modelname as common, without realisation info]

//combining different realisations of model modelname (with real_template: *->N)
void Combine_reals(string modelname);






#endif