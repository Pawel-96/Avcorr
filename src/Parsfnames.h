#ifndef PARSFNAMES_H
#define PARSFNAMES_H

#include "../lib/RWMathStat.h"
#include <chrono>
#include <regex>  //finding pattern in str

const string paramfile="param.txt";
const string logfile="Results/log.txt";
const int NRMAX=100000; //maximum nR or naxa*naxb
const int CMAX=2000000000; //total limit of Cmax
const double FULLSKYDEG2=41252.96125; //full sky area [deg^2]

extern int err_common; //error in common params
extern int err_ver; //error in params from one of versions
extern int err_conv; //error in conversion from string


extern const string VERSION;	    //code version: angular/BOX/BOX_ellipses/LC_ellipses

extern const int USE_HDF5; 	//HDF5 [1] or ASCII [0] input data
extern const string POS_DSET;		//dataset with positions (if USE_HDF5=1)
extern const vector<string> cols_pos;  //columns with positions (if USE_HDF5=0)

//********************common parameters:********************
extern const int moment_min;                                 			//must be 2
extern const int moment_max;                                			    //must be 9
extern const double kappa;	    //multiplicity factor for number of circles
extern const int Cmin;       //minimum number of circles drawn
extern const int Cmax;       //maximum number of circles drawn
extern const int nreals;	    //number of independent randoms sets (randoms must be at least C*nreals)
extern const int NPIX_min; 									            //minimum pixelization [NEW]
extern const int Mvar_Poisson;										//probing for Poisson error
extern const int Recalc;        //if =1, calculating moments again
extern const int Clean;        //if =1, cleaning in the end
extern const int ErrP;        //compute Poisson errors? [0/1]
extern const int combine_reals; //[0/1] are Data/* files models with different reals?
extern const string real_template; //how are realisations marked in names, e.g. _BOX*_

extern  string EXT;
extern string PATH;
extern const vector<string> Model;
extern const int nmodels;

extern const int Random_provided;  //random file provided? [0/1]
extern const string Random_file; //random file name

//********************parameters - angular:********************
extern const double ramin; 		//[deg]catalog rightascension lower range
extern const double ramax; 		//[deg]catalog rightascension lower range
extern const double decmin;		//[deg]catalog declination lower range
extern const double decmax; 	//[deg] catalog declination upper range

extern const double Areaf_read;
extern const double Areaf; //catalog sky area [deg2]

//********************parameters - angular and BOX:********************
extern const double Rmin; 		//[deg] smallest angular scale considered
extern const double Rmax; 		//[deg] biggest angular scale considered
extern const int nR; 				//number of angular scales considered
extern const double qR;        			//R multiplicity factor

//********************parameters - BOX and BOX_ellipses:********************
extern const double boxsize; 	//box size in the same units as random spheres


//********************parameters - BOX_ellipses and LC_ellipses:********************
extern const double axamin; 	//smallest semi-major axis
extern const double axamax; 	//biggest semi-major axis
extern const double axbmin; 	//smallest semi-minor axis
extern const double axbmax; 	//biggest semi-minor axis
extern const int naxa; 			//number of semi-major axes
extern const int naxb; 			//number of semi-minor axes
extern const double qaxa;        	//semi-major axis multiplicity factor
extern const double qaxb;        //semi-minor axis multiplicity factor
extern const int n_allax;										//number of grid points

//********************parameters - LC_ellipses:********************
extern const double DCMIN;		//minimum distance analyzed
extern const double DCMAX;		//maximum distance analyzed

/*-------------------------------------------------------------------------------------------------------------------------*/
extern const int ncols_outfile; //number of columns (check Args2grid)




string Fin(string model); //input

string Mout_onerad(string model, int nnsize); //moment output for one radius

string Merrout(string model); //output for one model and all moments

void Writelog(string text, string option="write"); //log

double Progress(int iimax); //counting overall progress

void Error(int rank, int &err, string msg); //raising error with message

int Error_param(int rank); //pointing errors in paramfile

#endif