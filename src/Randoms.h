#ifndef RANDOMS_H
#define RANDOMS_H

#include "../lib/RWMathStat.h"

#include "Parsfnames.h"
#include "CIC.h"


//creating temporary randoms in continuous field
void Continuous_randoms(string frand_temp, double R, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn);


//processing randoms for one scale
void Process_randoms_onerad(int rank, int comm_size, int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn,
vector<double> &Nav_circ, vector<double> &Nav_pix, int &nrand);


void Create_processed_randoms(int rank, int comm_size, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn,
vector<double> &Nav_circ, vector<double> &Nav_pix);


//angular CiC, but version especially for random processing
void CIC_angular_rnd(vector<vector<vector<double> > > &pix, string output, 
int npix_edge, int nnR, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand, int mode, double &av, double &sigma);


//CiC for LC_ellipses, but version especially for random processing
void CIC_LC_ellipses_rnd(vector<vector<vector<double> > > &pix, string output,
int npix_edge, int nnallax, vector<double> &Xcn, vector<double> &Ycn, vector<double> &Zcn, int &nrand, int mode, double &av, double &sigma);

#endif
