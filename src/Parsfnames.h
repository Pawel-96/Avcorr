#ifndef PARSFNAMES_H
#define PARSFNAMES_H

#include "../lib/RWMathStat.h"
#include <chrono>
#include <regex>  //finding pattern in str

const string paramfile="param.txt";
const string logfile="Results/log.txt";
const int NRMAX=100000; //maximum nR or naxa*naxb
const int CMAX=1000000000; //total limit of Cmax
const double FULLSKYDEG2=41252.96125; //full sky area [deg^2]

inline int err_common=0; //error in common params
inline int err_ver=0; //error in params from one of versions
inline int err_conv=0; //error in conversion from string


const string VERSION=Get_parameter(paramfile,"VERSION",err_common)[0];	    //code version: angular/BOX/BOX_ellipses/LC_ellipses

const int USE_HDF5=conv_check(Get_parameter(paramfile,"USE_HDF5",err_common)[0],err_conv); 	//HDF5 [1] or ASCII [0] input data
const string POS_DSET=Get_parameter(paramfile,"POS_DSET",err_common)[0];		//dataset with positions (if USE_HDF5=1)
const vector<string> cols_pos=Get_parameter(paramfile,"cols_pos",err_common);  //columns with positions (if USE_HDF5=0)

//********************common parameters:********************
const int moment_min=2;                                 			//must be 2
const int moment_max=9;                                			    //must be 9
const double kappa=conv_check(Get_parameter(paramfile,"kappa",err_common)[0],err_conv);	    //multiplicity factor for number of circles
const int Cmin=conv_check(Get_parameter(paramfile,"Cmin",err_common)[0],err_conv);       //minimum number of circles drawn
const int Cmax=conv_check(Get_parameter(paramfile,"Cmax",err_common)[0],err_conv);       //maximum number of circles drawn
const int nreals=conv_check(Get_parameter(paramfile,"nreals",err_common)[0],err_conv);	    //number of independent randoms sets (randoms must be at least C*nreals)
const int NPIX_min=50; 									            //minimum pixelization [NEW]
const int Mvar_Poisson=1000;										//probing for Poisson error
const int Recalc=conv_check(Get_parameter(paramfile,"Recalc",err_common)[0],err_conv);        //if =1, calculating moments again
const int Clean=conv_check(Get_parameter(paramfile,"Clean",err_common)[0],err_conv);        //if =1, cleaning in the end
const int ErrP=conv_check(Get_parameter(paramfile,"ErrPoisson",err_common)[0],err_conv);        //compute Poisson errors? [0/1]
const int combine_reals=conv_check(Get_parameter(paramfile,"Combine_reals",err_common)[0],err_conv); //[0/1] are Data/* files models with different reals?
const string real_template=Get_parameter(paramfile,"Real_template",err_common)[0]; //how are realisations marked in names, e.g. _BOX*_

inline string EXT="";
const vector<string> Model=Conditional_modelreading("Data",paramfile,EXT,"Datafiles");
const int nmodels=Model.size();

const int Random_provided=conv_check(Get_parameter(paramfile,"Random_provided",err_common)[0],err_conv);  //random file provided? [0/1]
const string Random_file=Get_parameter(paramfile,"Random_file",err_common)[0]; //random file name

//********************parameters - angular:********************
const double ramin=conv_check(Get_parameter(paramfile,"ramin",err_ver)[0],err_conv); 		//[deg]catalog rightascension lower range
const double ramax=conv_check(Get_parameter(paramfile,"ramax",err_ver)[0],err_conv); 		//[deg]catalog rightascension lower range
const double decmin=conv_check(Get_parameter(paramfile,"decmin",err_ver)[0],err_conv);		//[deg]catalog declination lower range
const double decmax=conv_check(Get_parameter(paramfile,"decmax",err_ver)[0],err_conv); 	//[deg] catalog declination upper range

const double Areaf_read=conv_check(Get_parameter(paramfile,"Areaf",err_ver)[0],err_conv);
const double Areaf=Areaf_read==-1 ? Area_spherical(ramin,ramax,decmin,decmax)*str2deg2 : Areaf_read; //catalog sky area [deg2]

//********************parameters - angular and BOX:********************
const double Rmin=conv_check(Get_parameter(paramfile,"Rmin",err_ver)[0],err_conv); 		//[deg] smallest angular scale considered
const double Rmax=conv_check(Get_parameter(paramfile,"Rmax",err_ver)[0],err_conv); 		//[deg] biggest angular scale considered
const int nR=conv_check(Get_parameter(paramfile,"nR",err_ver)[0],err_conv); 				//number of angular scales considered
const double qR=pow(10.,1.*log10(Rmax/Rmin)/nR);        			//R multiplicity factor

//********************parameters - BOX and BOX_ellipses:********************
const double boxsize=conv_check(Get_parameter(paramfile,"Boxsize",err_ver)[0],err_conv); 	//box size in the same units as random spheres


//********************parameters - BOX_ellipses and LC_ellipses:********************
const double axamin=conv_check(Get_parameter(paramfile,"axamin",err_ver)[0],err_conv); 	//smallest semi-major axis
const double axamax=conv_check(Get_parameter(paramfile,"axamax",err_ver)[0],err_conv); 	//biggest semi-major axis
const double axbmin=conv_check(Get_parameter(paramfile,"axbmin",err_ver)[0],err_conv); 	//smallest semi-minor axis
const double axbmax=conv_check(Get_parameter(paramfile,"axbmax",err_ver)[0],err_conv); 	//biggest semi-minor axis
const int naxa=conv_check(Get_parameter(paramfile,"naxa",err_ver)[0],err_conv); 			//number of semi-major axes
const int naxb=conv_check(Get_parameter(paramfile,"naxb",err_ver)[0],err_conv); 			//number of semi-minor axes
const double qaxa=pow(10.,1.*log10(axamax/axamin)/naxa);        	//semi-major axis multiplicity factor
const double qaxb=pow(10.,1.*log10(axbmax/axbmin)/naxb);        //semi-minor axis multiplicity factor
const int n_allax=naxa*naxb;										//number of grid points

//********************parameters - LC_ellipses:********************
const double DCMIN=conv_check(Get_parameter(paramfile,"DCMIN",err_ver)[0],err_conv);		//minimum distance analyzed
const double DCMAX=conv_check(Get_parameter(paramfile,"DCMAX",err_ver)[0],err_conv);		//maximum distance analyzed

/*-------------------------------------------------------------------------------------------------------------------------*/
const int ncols_outfile=ErrP==0?17:25; //number of columns (check Args2grid)



string Fin(string model); //input

string Mout_onerad(string model, int nnsize); //moment output for one radius

string Merrout(string model); //output for one model and all moments

void Writelog(string text, string option="write"); //log

double Progress(int iimax); //counting overall progress

void Error(int &err, string msg); //raising error with message

int Error_param(); //pointing errors in paramfile

#endif