#ifndef PARSFNAMES_H
#define PARSFNAMES_H

#include "../lib/RWMathStat.h"
#include<chrono>
#include <regex>  //finding pattern in str

const string paramfile="param.txt";
const string logfile="Results/log.txt";
const int NRMAX=100000; //maximum nR or naxa*naxb
const int CMAX=1000000000; //total limit of Cmax
const double FULLSKYDEG2=41252.96125; //full sky area [deg^2]

const string VERSION=Get_parameter(paramfile,"VERSION")[0];	    //code version: angular/BOX/BOX_ellipses/LC_ellipses

const int USE_HDF5=conv(Get_parameter(paramfile,"USE_HDF5")[0]); 	//HDF5 [1] or ASCII [0] input data
const string POS_DSET=Get_parameter(paramfile,"POS_DSET")[0];		//dataset with positions (if USE_HDF5=1)
const vector<string> cols_pos=Get_parameter(paramfile,"cols_pos");  //columns with positions (if USE_HDF5=0)


//********************common parameters:********************
const int moment_min=2;                                 			//must be 2
const int moment_max=9;                                			    //must be 9
const double kappa=conv(Get_parameter(paramfile,"kappa")[0]);	    //multiplicity factor for number of circles
const int Cmin=conv(Get_parameter(paramfile,"Cmin")[0]);       //minimum number of circles drawn
const int Cmax=conv(Get_parameter(paramfile,"Cmax")[0]);       //maximum number of circles drawn
const int nreals=conv(Get_parameter(paramfile,"nreals")[0]);	    //number of independent randoms sets (randoms must be at least C*nreals)
const int NPIX_min=50; 									            //minimum pixelization [NEW]
const int Mvar_Poisson=1000;										//probing for Poisson error
const int Recalc=conv(Get_parameter(paramfile,"Recalc")[0]);        //if =1, calculating moments again
const int Clean=conv(Get_parameter(paramfile,"Clean")[0]);        //if =1, cleaning in the end
const int ErrP=conv(Get_parameter(paramfile,"ErrPoisson")[0]);        //compute Poisson errors? [0/1]
const int combine_reals=conv(Get_parameter(paramfile,"Combine_reals")[0]); //[0/1] are Data/* files models with different reals?
const string real_template=Get_parameter(paramfile,"Real_template")[0]; //how are realisations marked in names, e.g. _BOX*_

inline string EXT="";
const vector<string> Model=Conditional_modelreading("Data",paramfile,EXT,"Datafiles");
const int nmodels=Model.size();

const int Random_provided=conv(Get_parameter(paramfile,"Random_provided")[0]);  //random file provided? [0/1]
const string Random_file=Get_parameter(paramfile,"Random_file")[0]; //random file name

//********************parameters - angular:********************
const double ramin=conv(Get_parameter(paramfile,"ramin")[0]); 		//[deg]catalog rightascension lower range
const double ramax=conv(Get_parameter(paramfile,"ramax")[0]); 		//[deg]catalog rightascension lower range
const double decmin=conv(Get_parameter(paramfile,"decmin")[0]);		//[deg]catalog declination lower range
const double decmax=conv(Get_parameter(paramfile,"decmax")[0]); 	//[deg] catalog declination upper range

const double Areaf_read=conv(Get_parameter(paramfile,"Areaf")[0]);
const double Areaf=Areaf_read==-1 ? Area_spherical(ramin,ramax,decmin,decmax)*str2deg2 : Areaf_read; //catalog sky area [deg2]

//********************parameters - angular and BOX:********************
const double Rmin=conv(Get_parameter(paramfile,"Rmin")[0]); 		//[deg] smallest angular scale considered
const double Rmax=conv(Get_parameter(paramfile,"Rmax")[0]); 		//[deg] biggest angular scale considered
const int nR=conv(Get_parameter(paramfile,"nR")[0]); 				//number of angular scales considered
const double qR=pow(10.,1.*log10(Rmax/Rmin)/nR);        			//R multiplicity factor

//********************parameters - BOX and BOX_ellipses:********************
const double boxsize=conv(Get_parameter(paramfile,"Boxsize")[0]); 	//box size in the same units as random spheres


//********************parameters - BOX_ellipses and LC_ellipses:********************
const double axamin=conv(Get_parameter(paramfile,"axamin")[0]); 	//smallest semi-major axis
const double axamax=conv(Get_parameter(paramfile,"axamax")[0]); 	//biggest semi-major axis
const double axbmin=conv(Get_parameter(paramfile,"axbmin")[0]); 	//smallest semi-minor axis
const double axbmax=conv(Get_parameter(paramfile,"axbmax")[0]); 	//biggest semi-minor axis
const int naxa=conv(Get_parameter(paramfile,"naxa")[0]); 			//number of semi-major axes
const int naxb=conv(Get_parameter(paramfile,"naxb")[0]); 			//number of semi-minor axes
const double qaxa=pow(10.,1.*log10(axamax/axamin)/naxa);        	//semi-major axis multiplicity factor
const double qaxb=pow(10.,1.*log10(axbmax/axbmin)/naxb);        //semi-minor axis multiplicity factor
const int n_allax=naxa*naxb;										//number of grid points

//********************parameters - LC_ellipses:********************
const double DCMIN=conv(Get_parameter(paramfile,"DCMIN")[0]);		//minimum distance analyzed
const double DCMAX=conv(Get_parameter(paramfile,"DCMAX")[0]);		//maximum distance analyzed

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